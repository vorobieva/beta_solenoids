#! /usr/bin/env python
InfoString = ''' 
TO GENERATE BACKBONES AROUND STARTING STRUCTURE 
    WHICH MAY BE MODIFIED 
BY STRECHING REPEAT REGION DESIGNATED BY NUMBER
'''

# '''
# from repo 
import solenoid_tools

# libraries
from scipy import spatial
import multiprocessing
import numpy as np
import subprocess
import argparse
import glob
import sys
import os
import re

import rosetta
rosetta.init()
# rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
from rosetta.protocols import grafting 

from rosetta.core.pose import initialize_atomid_map
from rosetta.core.id import AtomID, AtomID_Map_Real, AtomID_Map_bool

pymol_link = rosetta.PyMolMover()

sys.argv.extend(['-pdbs', '1EZG.pdb', '-out', './' ])
# '''

def fuse(Pose1, Pose2):


def append_poses(RepeatPose, JunctionPose, Duplications):
  ''' 
    Give two overlapping poses
  '''

  RepetitivePose = rosetta.Pose()
  RepetitivePose.assign(RepeatPose)
  
  # for Duplication in range(Duplications):
  RepeatUnitCopy = rosetta.Pose()
  RepeatUnitCopy.assign(RepeatPose)

  MatchingResidues = solenoid_tools.match_hotspots_and_motifs_to_dock(RepeatPose, JunctionPose, 0.1)

  print 'MatchingResidues', MatchingResidues

  OverhangStart = MatchingResidues[RepeatPose.n_residue()][0] + 1

  OverhangCoordList = [] # for coords of overhang atoms
  RepeatCoordList = [] # for coords of begining of repeat, which will be superimposed onto those of the overhang atoms to translate protein

  for Repeat0index, OverhangPosition in enumerate( range(OverhangStart, JunctionPose.n_residue()+1) ):
    RepeatPosition = Repeat0index + 1

    for AtomName in ['N', 'C', 'CA']:
      OverhangCoordList.append( list(JunctionPose.residue(OverhangPosition).xyz(AtomName)) )
      RepeatCoordList.append( list(RepeatPose.residue(RepeatPosition).xyz(AtomName)) )

    print RepeatPosition, OverhangPosition

  OverHangArray = np.array(OverhangCoordList)
  RepeatArray = np.array(RepeatCoordList)

  RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( OverHangArray, RepeatArray )
  print 'RMSD:', RMSD
  rosetta.Pose.apply_transform_Rx_plus_v(RepeatUnitCopy, rMtx, tVec)


  if Duplications: # Reset and call recursively ?

    MovedRepeatIdenticalCoordList = []
    JunctionIdenticalCoordList = []
    
    for Res in range( RepeatPose.n_residue() ):
      if len(MatchingResidues[Res]):
        assert len(MatchingResidues[Res]) == 1
        
        for AtomName in ['N', 'C', 'CA']:
          MovedRepeatIdenticalCoordList.append( list(RepeatUnitCopy.residue(Res).xyz(AtomName)) )
          MovedRepeatIdenticalCoordList.append( list(JunctionPose.residue(MatchingResidues[Res][0]).xyz(AtomName)) )

    RepeatIdenticalArray = np.array(MovedRepeatIdenticalCoordList)
    JunctionIdenticalArray = np.array(JunctionIdenticalCoordList)

    RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( RepeatIdenticalArray, JunctionIdenticalArray )
    print 'RMSD:', RMSD  
    rosetta.Pose.apply_transform_Rx_plus_v(JunctionPose, rMtx, tVec)

    Duplications -= 1
    append_poses( JunctionPose, Duplications)

  else: 

    # return rMtx, tVec 

def apply_motifs(SubjectPose, MotifPose, MotifsOverSubject):
   ''' takes subject pose, pose of motif, and list length of subject if 'position' pointer to motif pose residue to apply to that index's position '''
   for Position in xrange( 1, SubjectPose.n_residue()+1 ):

      if len(MotifsOverSubject[Position]) == 1:
         # replace_residue( (Pose)arg1, (int)seqpos, (Residue)new_rsd_in, (vector1_pair_string_string)atom_pairs) -> None :
         SubjectPose.replace_residue(Position, MotifPose.residue( MotifsOverSubject[Position][0] ), True)

      elif len(MotifsOverSubject[Position]) > 1:
         print 'Randomly selecting residue for position %d '%Position
         SubjectPose.replace_residue(Position, MotifPose.residue( random.choice( MotifsOverSubject[Position] ) ), True)


def repeat_pose_generator(Pose, RepeatStart, RepeatSpacings):
  '''  '''
  EndShift = RepeatSpacings[0] - 1
  RepeatEnd = RepeatStart + EndShift

  RepeatBase = grafting.return_region(Pose, RepeatStart, RepeatEnd)

  return RepeatBase

def mcmc_mover(Pose, ScoreFunction, kT=1.0, Gen=100):
  ''' threadable backbone breathing and side chain repacking with constraints '''
  print 'Starting work on', Pose

  movemap = rosetta.MoveMap()
  movemap.set_bb(True)
  small_mover = rosetta.SmallMover(movemap, kT, 1)
  shear_mover = rosetta.ShearMover(movemap, kT, 1)

  MontyCarlos = rosetta.MonteCarlo(Pose, ScoreFunction, kT)

  small_mover(Pose)

  for g in range(Gen):
    small_mover(Pose)
    MontyCarlos.boltzmann(Pose)
    shear_mover(Pose)
    MontyCarlos.boltzmann(Pose)


  print MontyCarlos.show_scores()
  print MontyCarlos.show_counters()
  print MontyCarlos.show_state()


class alpha_carbon:
  """ Calpha node for wDAG searches of proteins for repeats """
  def __init__(self, Number, CoordinateArray):
    if type(CoordinateArray) == list:
      CoordinateArray = np.array(CoordinateArray)
    self.CoordinateArray = CoordinateArray
    self.Number = int(Number)
    self.DownstreamNeighbors = {} # keyed with residue number, value is displacement vector to that residue's own CA

class pose_wdag:
  """ 'wDag' is a little generous. Init method takes protein and spawns Calpha instances to populate """
  def __init__( self, Pose, MinRepeats=2, DistanceTarget=4.7, DistanceFlex=0.5, AngleFlex=15.0 ):
    # Pose = Pose
    self.MinRepeats = MinRepeats
    self.DistanceTarget = DistanceTarget
    self.DistanceFlex = DistanceFlex
    self.AngleFlex = np.radians(AngleFlex) # converts to radians

    # make instance of alpha_carbon for each row
    self.CalphaInstances = []
    CalphaCoords = []
    for P in range( 1, Pose.n_residue() + 1):
      print P
      print Pose.residue(P)
      print Pose.residue(P).xyz('CA')
      CoordList = [Value for Value in Pose.residue(P).xyz('CA') ]
      print CoordList
      CoordList = list( Pose.residue(P).xyz('CA') )
      print CoordList
      print '\n*2'
      self.CalphaInstances.append( alpha_carbon(P, CoordList) )
      CalphaCoords.append( CoordList )

    self.CalphaArray = np.array(CalphaCoords)


  def find_repeat_chains(self):
    CheckRedundancyChecker = {}
    Repeats = []
    Count = 0
    for Calpha in self.CalphaInstances:
      for NeighborNumber in Calpha.DownstreamNeighbors:
        '''
        Repeat chains start with each of Calphas downstream neighbors 
        and are recursively extended so long as VecA ~= VecB
        
           VecA      VecB
        CA- - - > CA- - - > CA
        '''
        Count += 1
        RepeatChain = [ Calpha.Number ]
        RecursiveCalpha = Calpha
        DownstreamNumber = NeighborNumber
        
        try:
          CheckRedundancyChecker['%d_%d'%(Calpha.Number, DownstreamNumber)]

        except KeyError:
          ChainExtended = 1

        while ChainExtended:
          RepeatChain.append(DownstreamNumber)
          DownstreamCalpha = self.CalphaDict[DownstreamNumber]          
          VectorA = RecursiveCalpha.DownstreamNeighbors[DownstreamNumber]
          
          ChainExtended = 0 # not extended unless valid extension found
          for NextDownstreamNumber in DownstreamCalpha.DownstreamNeighbors:
            
            # To prevent redundant calculations, upstream downstream pairs already investigated aren't considered         
            try:
              CheckRedundancyChecker['%d_%d'%(DownstreamNumber, NextDownstreamNumber)]
            # only residue pairs not in redundancy checking dictionary are considering
            except KeyError:
              VectorB = DownstreamCalpha.DownstreamNeighbors[NextDownstreamNumber]
            
              if spatial.distance.cosine(VectorA, VectorB) < self.AngleFlex:
                CheckRedundancyChecker['%d_%d'%(DownstreamNumber, NextDownstreamNumber)] = True
                DownstreamNumber = NextDownstreamNumber
                RecursiveCalpha = DownstreamCalpha
                ChainExtended = 1
                break

        if len(RepeatChain) >= self.MinRepeats:
          Repeats.append(RepeatChain)
          # print 'select repeat%d, resi %s'%(Count,'+'.join([str(Number).split('.')[0] for Number in RepeatChain]))

    return Repeats

  def find_downstream_neighbors(self):
    ''' To set up DAG Calphas become nodes'''
    self.CalphaDict = {}
    for i, Instance in enumerate(self.CalphaInstances):
      # To make the graph directed, repeat chains will always be arranged from 
      for Subsequent in self.CalphaInstances[i:]:
        # Displacement from lower numbered residue to higher number
        Displacement = Subsequent.CoordinateArray - Instance.CoordinateArray
        # Magnitute of displacement vector
        Distance = (np.sum([Component**2 for Component in Displacement]))**0.5
        # check if displacement of given Instance (upstream) and Subsequent (downstream) has magnitute within target range
        if np.abs(Distance - self.DistanceTarget) <= self.DistanceFlex:
          # if in range record Subsequent as downstream neighbor 
          Instance.DownstreamNeighbors[Subsequent.Number] = Displacement
      # add residue alpha carbon instance to dictionary
      self.CalphaDict[Instance.Number] = Instance



def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' nc_cst_gen.py arguments ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=list, help=' input pdbs ', required=True)
  ArgParser.add_argument('-start', type=list, help=' repeat start residue ', default=False)
  # ArgParser.add_argument('-end', type=list, help=' input pdbs ', required=True)
  ArgParser.add_argument('-repeat', type=list, help=' number of repeats to make ', default=3)
  ArgParser.add_argument('-dat', type=str, help=' input symm dat file ', default=None)
  ArgParser.add_argument('-cst', type=str, help=' input symm cst file ', default=None)
  ArgParser.add_argument('-out', type=str, help=' output directory ', default='./')
  Args = ArgParser.parse_args()
  if len(Args.pdbs[0]) == 1:
    Args.pdbs = [''.join(Args.pdbs)]
  
  # return Args
  # Super list to hold multithreaded design trajectories from all input pdbs
  PdbJobs = []

  # ScoreFunction = rosetta.create_score_function("standard")
  # ScoreFunction.set_weight(rosetta.atom_pair_constraint, 5.0)

  for Pdb in Args.pdbs:
    
    # load Pdb
    Pose = rosetta.pose_from_pdb(Pdb)

    # Non-rosetta part, together these four lines finds primary sequence repeats to duplicate     
    wDag = pose_wdag(Pose) # default settings are for beta solenoids !
    wDag.find_downstream_neighbors()
    RepeatChains = wDag.find_repeat_chains()
    ConsolidatedRepeatStarts, TandemRepeats = solenoid_tools.consolidate_repeats(RepeatChains)

    RepeatPoses = []
    for RepeatStart in ConsolidatedRepeatStarts:
      
      RepeatSegmentPose = repeat_pose_generator(Pose, RepeatStart, TandemRepeats[RepeatStart])
      
      for Res in range(1, RepeatSegmentPose.n_residue()+1):
        if RepeatSegmentPose.residue(Res).name() == 'CYD':
          rosetta.core.conformation.change_cys_state(Res, 'CYS', RepeatSegmentPose.conformation())

    #   append_poses

      RepeatPoses.append(Pose)
    return RepeatPoses

    # Rosetta part
    # mcmc_mover(Pose)

    # Sub list to hold multithreaded design trajectories from all frames within current itertaions pdb
    # FrameJobs = []
    # if Arg.start:
    #   Arg.start
    # else:
      # for Frame in Frames:
      #   SingleProcess = multiprocessing.Process(target=worker, args=(i,))
      #   FrameJobs.append(SingleProcess)
      #   SingleProcess.start()

    # PdbJobs.append(FrameJobs)

# if __name__ == "__main__":
#   sys.exit(main())