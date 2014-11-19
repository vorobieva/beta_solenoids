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
import itertools
import numpy as np
import subprocess
import argparse
import glob
import sys
import os
import re

import rosetta
# rosetta.init()
rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
from rosetta.protocols import grafting 

from rosetta.core.pose import initialize_atomid_map
# from rosetta.core.id import AtomID, AtomID_Map_Real, AtomID_Map_bool

# pymol_link = rosetta.PyMolMover()

# ''' 
# sys.argv.extend(['-pdbs', '1EZG.pdb', '-out', './' ])


def dump_many_poses(IterablePoses, Tag):
  for i, Pose in enumerate(IterablePoses):
    rosetta.dump_pdb( Pose, '%s_n%d.pdb'%(Tag, (1+i)) )

def pose_repeat_unit_finder(Pose, RepeatChains=False):
  '''  '''
  if RepeatChains == False:
    # Non-rosetta part, together these four lines finds primary sequence repeats to duplicate     
    wDag = solenoid_tools.pose_wdag(Pose) # default settings are for beta solenoids !
    wDag.find_downstream_neighbors()
    RepeatChains = wDag.find_repeat_chains()
  print 'RepeatChains', RepeatChains
    
  ConsolidatedRepeatStarts, TandemRepeats = solenoid_tools.consolidate_repeats(RepeatChains)
  return ConsolidatedRepeatStarts, TandemRepeats

def fuse(Pose1, Pose2, SubsetSize=2):
  # Should continue to fiddle with the hardcoded var below,
  # Originally 0.5, only good for indentical copies,
  # then 1.5, works for close copy
  # trying 2.0
  MatchingResidueHash = solenoid_tools.match_superimposed_pose_residues(Pose1, Pose2, 1.5)
  # checks there is a one to one correspondance for all residue matches
  for MatchRes in MatchingResidueHash:
    assert len(MatchingResidueHash[MatchRes]) <= 1

  # list comprehension through matches, add one for each position with match
  NofMatchRes = sum([ 1 for Match in MatchingResidueHash if len(MatchingResidueHash[Match]) ])
  try:
    assert SubsetSize <= NofMatchRes
  except AssertionError:
    dump_many_poses([Pose1, Pose2], 'FusedFusion')
    # print 'MatchingResidueHash:', MatchingResidueHash
    assert SubsetSize <= NofMatchRes, ' Designated subset length should not exceed that of the overlap between poses. Poses dumped for inspection '
  # contains positions like [ (MatchRes1InPose1, MatchRes1InPose2), (MatchRes2InPose1, MatchRes2InPose2) .. ]
  CorrespondingResidues = []

  # iterates through positions in pose1
  for P1 in range( 1, Pose1.n_residue()+1 ):
    if len(MatchingResidueHash[P1]) == 1:
      P2 = MatchingResidueHash[P1][0]
      CorrespondingResidues.append((P1, P2))

  LengthIncenative = 1.2
  BestRMSD = 999
  BestSubset = []
  BestTransformation = ()

  for i in range( len(CorrespondingResidues) - SubsetSize + 1):
    Pose1Coords = []
    Pose2Coords = []    
    IterationsSubset = CorrespondingResidues[i:i+SubsetSize]
    for ResidueMatch in IterationsSubset:
      P1 = ResidueMatch[0]
      P2 = ResidueMatch[1]
      for AtomName in ['N','C','O','CA']:
        Pose1Coords.append( list(Pose1.residue(P1).xyz(AtomName)) )
        Pose2Coords.append( list(Pose2.residue(P2).xyz(AtomName)) )      
    # makes (subset length)by3 array out of list of lists
 
    Pose1Array = np.array(Pose1Coords)
    Pose2Array = np.array(Pose2Coords)

    RMSD, Rotation, Translation = solenoid_tools.rmsd_2_np_arrays_rosetta(Pose1Array, Pose2Array)

    if RMSD < BestRMSD:
      BestRMSD = RMSD
      BestSubset = IterationsSubset
      BestTransformation = ( Rotation, Translation )

  # Unpack within overlap subset and corresponding transformation vectors 
  Rotation, Translation = BestTransformation
  rosetta.Pose.apply_transform_Rx_plus_v(Pose2, Rotation, Translation)

  # print 'BestSubset', BestSubset
  # make a lot of sense for even overlaps, makes less sense for odd overlaps
  Cutpoint = SubsetSize / 2
  EndOfPose1 = BestSubset[Cutpoint-1][0]
  StartOfPose2 = BestSubset[Cutpoint][1]

  # print 'EndOfPose1', EndOfPose1
  # print 'StartOfPose2', StartOfPose2
  FusionPose = grafting.return_region(Pose1, 1, EndOfPose1)
  # rosetta.dump_pdb(FusionPose, 'FusionPose1.pdb')
  for Pose2Position in range( StartOfPose2, Pose2.n_residue()+1 ):
    Pose2Residue = Pose2.residue(Pose2Position)
    FusionPose.append_residue_by_bond(Pose2Residue)
  
  return FusionPose, BestRMSD, CorrespondingResidues


def get_residue_array(Pose, Residues):
  CoordList = [] 
  for Residue in Residues:
    for AtomName in ['N','C','O','CA']:
      CoordList.append( list(Pose.residue(Residue).xyz(AtomName)) )
  return np.array(CoordList)


def extrapolate_repeat_pose(Repeat1Pose, Repeat2Pose, Duplications):
  '''  MUST give two overlapping, equal length, one-repeat-unit long poses !!!  '''
  assert Repeat1Pose.n_residue() == Repeat2Pose.n_residue(), ' Repeat poses must be same length '
  RepeatLength = Repeat1Pose.n_residue()
  Residues = [P for P in range(1, RepeatLength+1)]
 
  # First fusion is easy since the input poses contain duplicate residues
  FusionPose, RMSD, CorrespondingResidues = fuse(Repeat1Pose, Repeat2Pose)
  
  # Copy and store first subunit as base unit for subsequent transformations/fusions 
  FusionBasePose = rosetta.Pose()
  FusionBasePose.assign(FusionPose)
  # rosetta.dump_pdb(FusionBasePose, 'FusionBasePose.pdb')

  # Umatched region should be N terminal region of repeat pose 1 without matching pose 2 residues
  # Used repeated during loop below to superimpose on to end of growing fusion pose
  UnmatchedPose1Residues = [ Number for Number in range(1, CorrespondingResidues[0][0]) ]

  # N-terminal end of repeat unit pose 1 that extends past repeat unit pose 2
  NterminalUnmatchedArray = get_residue_array(Repeat1Pose, UnmatchedPose1Residues)

  # one iteration is run for each duplication event, 
  # each time add a fusion pose to end of growing 
  for Duplication in range(Duplications):
    # Grabs residue from back of RepeatPose2 for each unmatched residue at begining of Pose1
    EndMatch = []
    # also start making array of coord
    for UnmatchResidue in UnmatchedPose1Residues:
      EndMatch.append(Residues[-1*UnmatchResidue])

    EndMatch.reverse()
    EndMatchArray = get_residue_array(Repeat2Pose, EndMatch)

    # print 'EndMatch', EndMatch
    RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( EndMatchArray, NterminalUnmatchedArray )
    # print 'Corse RMSD:  ', RMSD
    # print 'rMtx:  ', rMtx
    # print 'tVec:  ', tVec

    DuplicatePose = rosetta.Pose()
    DuplicatePose.assign(FusionBasePose)
    rosetta.Pose.apply_transform_Rx_plus_v(DuplicatePose, rMtx, tVec)
    # rosetta.dump_pdb(DuplicatePose, 'Dup%d_DuplicatePose.pdb'%Duplication)

    FusionPose, RMSD, CorrespondingResidues = fuse(FusionPose, DuplicatePose)
    # print 'Refined RMSD', RMSD
    # rosetta.dump_pdb(FusionPose, 'Dup%d_FusionPose.pdb'%Duplication)
    
    NterminalUnmatchedArray = get_residue_array(FusionPose, UnmatchedPose1Residues)
    
    Residues = [ P for P in range(1, FusionPose.n_residue() + 1) ]
    Repeat2Pose = rosetta.Pose()
    Repeat2Pose.assign(FusionPose)
    # rosetta.dump_pdb(Repeat2Pose, 'Dup%d_NewRepeat2Pose.pdb'%Duplication)

  return FusionPose

class constraint_extrapolator:
  """ constraint extrapolator"""
  def __init__(self, CstFilename):
    self.CstFilename = CstFilename
    # contains: ConstraintDict[Res1] = {Name1: (Name2, Res2, ConstraintParameters)} 
    #  and      ConstraintDict[Res2] = {Name2: (Name1, Res1, ConstraintParameters)}
    self.CstDictionary = self.parse_cst(CstFilename)

  def parse_cst(self, CstFilename):
    with open(CstFilename, 'r') as CstFile:
      CstLines = CstFile.readlines()
    ConstraintDict = {}
    for Line in CstLines:
      if Line.startswith('AtomPair'):
        CstType, Name1, Res1, Name2, Res2 = tuple( Line.split()[:5] )
        ConstraintParameters = Line.split()[5:]
        
        Res1 = int(Res1)
        Res2 = int(Res2)

        # these check for residue entry in constraint dict
        try:
          check = ConstraintDict[Res1]
        except KeyError:
          ConstraintDict[Res1] = {}
        try:
          check = ConstraintDict[Res2]
        except KeyError:
          ConstraintDict[Res2] = {}

        # add atom name entries for residue subdict as needed
        try:
          ConstraintDict[Res1][Name1].append( (Name2, Res2, ConstraintParameters) )
        except KeyError:
          ConstraintDict[Res1][Name1] = [ (Name2, Res2, ConstraintParameters) ] 
        try:
          ConstraintDict[Res2][Name2].append( (Name1, Res1, ConstraintParameters) )
        except KeyError:
          ConstraintDict[Res2][Name2] = [ (Name1, Res1, ConstraintParameters) ]
    
    return ConstraintDict

  def reassemble_atompair_cst(self, Name1, Res1, Name2, Res2, ConstraintParameters):
    return 'AtomPair %s %d %s %d %s'%(Name1, Res1, Name2, Res2, ' '.join(ConstraintParameters) )

  def shift(self, Position):
    return Position - self.NewPoseStartShift

  def in_range(self, Position):
    if Position < self.Range[0]:
      return False
    elif Position > self.Range[1]:
      return False
    else:
      return True

  def extrapolate_from_repeat_unit(self, StartArchetype, EndArchetype, RepeatUnitLength, NewPose):
    ''' renumbers based on repeat unit pose '''
    self.Range = (1, NewPose.n_residue())
    self.NewPoseStartShift = StartArchetype - 1 # for 1 indexing
    # Loop through positions in range of archetype
    # To avoid double counting first only add constraints from archetype residues to 
    # more C-terminal residues 
    
    # To hold extrapolated constraints
    NewCsts = []
    NewResNum = NewPose.n_residue()
    UnitShiftMultiples = (NewResNum / RepeatUnitLength) + 1
    UnitShiftList = [ RepeatUnitLength * Multiple for Multiple in range( UnitShiftMultiples+1 ) ] 
    # print 'UnitShiftList', UnitShiftList

    for ArchetypePosition in range(StartArchetype, EndArchetype+1):
      try:
        ArchetypePositionCstDict = self.CstDictionary[ArchetypePosition]
      except:
        continue
      # print        
      # print 'ArchetypePosition: ', ArchetypePosition
      # print 'shifted: ', self.shift(ArchetypePosition)
      # print 'ArchetypePositionCstDict: ', ArchetypePositionCstDict
      
      for AtomName in ArchetypePositionCstDict:
        # print 'AtomName: ', AtomName
        for Constraint in ArchetypePositionCstDict[AtomName]:
          OtherName, OtherPosition, ConstraintParameters = Constraint
          # print 'OtherPosition', OtherPosition
          if OtherPosition > ArchetypePosition:
            for UnitShift in UnitShiftList:
              ShiftedArchetype = self.shift(ArchetypePosition) + UnitShift
              ShiftedOther = self.shift(OtherPosition) + UnitShift
              
              if self.in_range(ShiftedArchetype) and self.in_range(ShiftedOther):
                
                # perhaps instead of removing these constraints, the residue should be changed. This would require cst and pose extrapolation to be intergrated
                if NewPose.residue(ShiftedArchetype).has(AtomName) and NewPose.residue(ShiftedOther).has(OtherName):
                  CST = self.reassemble_atompair_cst(AtomName, ShiftedArchetype, OtherName, ShiftedOther, ConstraintParameters)
                  # print CST
                  NewCsts.append(CST)
    
    # print ' debug exit '
    # sys.exit()         
    return NewCsts

def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' generate_backbones.py ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=str, nargs='+', help=' Input pdbs ', required=True)
  # ArgParser.add_argument('-start', type=int, help=' Repeat start residue ', default=False)
  
  # ArgParser.add_argument('-end', type=str, help=' input pdbs ', required=True)
  ArgParser.add_argument('-repeat', type=int, help=' Number of repeats to make ', default=5)
  ArgParser.add_argument('-max_turns_per_repeat', type=int, help=' Upper bound of number of turns per repeat unit to extrapolate repeat pose from. Lower bound always is 1 ', default=2)
  
  ArgParser.add_argument('-man_repeat', type=str, help=' manual input of repeat chains (i.e. stripes of residues along pdb), _ seperates residues in chain, __ seperates chains ', default=False)

  # ArgParser.add_argument('-dat', type=str, help=' input symm dat file ', default=None)
  ArgParser.add_argument('-csts', type=str, nargs='+', help=' Input symm cst files (MUST be provided in same order as input PDBS)', default=None)
  ArgParser.add_argument('-out', type=str, help=' Output directory ', default='./')

  Args = ArgParser.parse_args()
  # print 'Args.pdbs', Args.pdbs
  if len(Args.pdbs) == 1:
    Args.pdbs = [''.join(Args.pdbs)]
  if Args.csts:
    if len(Args.csts[0]) == 1:
      Args.csts = [''.join(Args.csts)]

  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  # print 'Args.pdbs', Args.pdbs
  
  ### extrapolate from input cst to new extended backbones

  for i, Pdb in enumerate(Args.pdbs):
    print 'Pdb:', Pdb
    # load Pdb
    Pose = rosetta.pose_from_pdb(Pdb)
    Pose.pdb_info(rosetta.core.pose.PDBInfo( Pose ))
    # get name base for output pdbs 
    InputPdbStem = Pdb.replace('.pdb', '')
    
    # get cst file if csts files
    if Args.csts:
      Cst = Args.csts[i]
      print 'Cst:', Cst
      Constrainer = constraint_extrapolator(Cst)

    # Get repeat unit poses from function above
    if Args.man_repeat == False:
      ConsolidatedRepeatStarts, TandemRepeats = pose_repeat_unit_finder(Pose)
    else:
      RepeatChains = Args.man_repeat.split('__')
      RepeatChains = [ [ int(Number) for Number in Chain.split('_') ] for Chain in RepeatChains]
      # print 'RepeatChains', RepeatChains
      # sys.exit()
      ConsolidatedRepeatStarts, TandemRepeats = pose_repeat_unit_finder(Pose, RepeatChains)
    
    ConsolidatedRepeatStarts.extend([45,46,47])
    print 'ConsolidatedRepeatStarts', ConsolidatedRepeatStarts
    print 'TandemRepeats', TandemRepeats

    InputPoseRepeatNumber = len(TandemRepeats[ConsolidatedRepeatStarts[0]])
    MaxTurns = min([InputPoseRepeatNumber, Args.max_turns_per_repeat])
    BaseRepeatUnitLength = TandemRepeats[ConsolidatedRepeatStarts[0]][0]
    UniformLength = Args.repeat * BaseRepeatUnitLength
    # MaxReach = InputPoseRepeatNumber * BaseRepeatUnitLength
    # print 'BaseRepeatUnitLength', BaseRepeatUnitLength
    # print 'UniformLength:', UniformLength

    ExtrapolatedPoses = []
    for RepeatUnitCombo in itertools.combinations(ConsolidatedRepeatStarts, Args.max_turns_per_repeat):
      print 'RepeatUnitCombo', RepeatUnitCombo
      RepeatUnit1Start, RepeatUnit2Start = RepeatUnitCombo
      assert RepeatUnit1Start <= RepeatUnit2Start, ' RepeatUnit1 must begin before RepeatUnit2 '
      if (RepeatUnit1Start + 4) <= RepeatUnit2Start <= (RepeatUnit1Start + BaseRepeatUnitLength - 4):  
        for TurnsPerRepeat in range(1, MaxTurns+1):
          UnitLength = (BaseRepeatUnitLength * TurnsPerRepeat)
          UnitShift = UnitLength - 1
          AvailableShifts = InputPoseRepeatNumber - TurnsPerRepeat
          # print 'TurnsPerRepeat', TurnsPerRepeat
          # print 'UnitLength', UnitLength
          # print 'AvailableShifts', AvailableShifts
          # print 'MaxReach', MaxReach
          # print 'MaxReach1 %d %d'%(RepeatUnit1Start, RepeatUnit1Start+MaxReach)
          # print 'MaxReach2 %d %d'%(RepeatUnit2Start, RepeatUnit2Start+MaxReach)
          
          for Shift in range(AvailableShifts+1):
            ShiftedUnit1Start = RepeatUnit1Start + Shift
            ShiftedUnit1End = ShiftedUnit1Start + UnitShift
            ShiftedUnit2Start = RepeatUnit2Start + Shift
            ShiftedUnit2End = ShiftedUnit2Start + UnitShift
 
            Repeat1Unit = grafting.return_region(Pose, ShiftedUnit1Start, ShiftedUnit1End)
            Repeat2Unit = grafting.return_region(Pose, ShiftedUnit2Start, ShiftedUnit2End)

            # print 'ShiftedUnit1Start, ShiftedUnit1End: ', ShiftedUnit1Start, ShiftedUnit1End
            # print 'ShiftedUnit2Start, ShiftedUnit2End: ', ShiftedUnit2Start, ShiftedUnit2End

            # use function to extrapolate from a partial repeat 
            Extrapolation = extrapolate_repeat_pose(Repeat1Unit, Repeat2Unit, Args.repeat - 1)
            # trim down to uniform length 
            Extrapolation = grafting.return_region(Extrapolation, 1, UniformLength)

            ExtrapolatedPoses.append(Extrapolation)
            rosetta.dump_pdb( Extrapolation, '%s%s_extra_%d.pdb'%(Args.out, InputPdbStem, len(ExtrapolatedPoses)) )

            # if cst provided with pdb, extrapolate constraints to corrspond to new backbone
            if Args.csts:
              ExtrapolatedCsts = Constrainer.extrapolate_from_repeat_unit(ShiftedUnit1Start, ShiftedUnit2End, UnitLength, Extrapolation)
              with open('%s%s_extra_%d.cst'%(Args.out, InputPdbStem, len(ExtrapolatedPoses)), 'w') as CstFile:
                print>>CstFile, '\n'.join(ExtrapolatedCsts) 

            # return 
   
        #   print '\n'
        # print '\n'*3

    print 'Number of extrapolated poses:', len(ExtrapolatedPoses)

    # for Repeat1Start in ConsolidatedRepeatStarts:
    #   # for Repeat1Start in ConsolidatedRepeatStarts:

    #   EndShift = TandemRepeats[RepeatStart][0] - 1    
    #   RepeatEnd = RepeatStart + EndShift

    #   RepeatPoses.append(RepeatUnit)
    
    # return RepeatPoses
    
    # for RepeatUnit in RepeatUnitPoses:
    #   for Res in range(1, RepeatUnit.n_residue()+1):
    #     if RepeatUnit.residue(Res).name() == 'CYD':
    #       rosetta.core.conformation.change_cys_state(Res, 'CYS', RepeatUnit.conformation())

    # # RepeatPoses.append(RepeatUnitPoses)
    # return RepeatUnitPoses
    
    # extrapolate_repeat_pose(Repeat1Pose, Repeat2Pose, Duplications)

    # # Non-rosetta part, together these four lines finds primary sequence repeats to duplicate     
    # wDag = pose_wdag(Pose) # default settings are for beta solenoids !
    # wDag.find_downstream_neighbors()
    # RepeatChains = wDag.find_repeat_chains()
    # ConsolidatedRepeatStarts, TandemRepeats = solenoid_tools.consolidate_repeats(RepeatChains)

    # RepeatPoses = []
    # for RepeatStart in ConsolidatedRepeatStarts:
    # # for RepeatSegmentPose in RepeatSegments:
    #   RepeatSegmentPose = pose_repeat_unit_generator(Pose, RepeatStart, TandemRepeats[RepeatStart])
      
    #   for Res in range(1, RepeatSegmentPose.n_residue()+1):
    #     if RepeatSegmentPose.residue(Res).name() == 'CYD':
    #       rosetta.core.conformation.change_cys_state(Res, 'CYS', RepeatSegmentPose.conformation())

    #   RepeatPoses.append(RepeatSegmentPose)
    # return RepeatPoses

if __name__ == "__main__":
  sys.exit(main())
