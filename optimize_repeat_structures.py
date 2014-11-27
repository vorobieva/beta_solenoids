#! /usr/bin/env python
InfoString = ''' 
To optimize repeat backbone with NCS torsion constraint and other constraints'''

# '''
# # from repo 
# import solenoid_tools

# libraries
from multiprocessing import Pool
import numpy as np
import subprocess
import argparse
import glob
import time
import sys
import os
import re

if '-h' not in sys.argv:
  import rosetta
  # rosetta.init()
  rosetta.init(extra_options = " -mute protocols") # -mute basic -mute core
  from rosetta.protocols.simple_moves import symmetry

# from rosetta.protocols import grafting 

'''
sys.argv.extend(['-pdb_stem', '3ult_rep', '-thread', '18'])
# '''


def optimize_repeat_pdb( (Pdb, Cst, RepeatLength) ):
  if Cst:
    print ' Relaxing %s with %s '%(Pdb, Cst)
  else:
    print ' Relaxing %s without atom pair constraints '%(Pdb)

  # idealize peptide bonds with command line subprocess
  subprocess.check_output(['idealize_jd2.default.linuxgccrelease', '-s', Pdb])
  IdealizedPdb = Pdb.replace('.pdb', '_0001.pdb') 
  time.sleep(0.5)

  Pose = rosetta.pose_from_pdb(IdealizedPdb)
  rosetta.dump_pdb(Pose, Pdb.replace('.pdb', '_ideal_presymm.pdb') )

  subprocess.check_output(['rm', IdealizedPdb])

  PoseLength = Pose.n_residue()
  # print 'RepeatLength', RepeatLength
  # print 'PoseLength', PoseLength

  assert PoseLength % RepeatLength == 0, 'pdb input into optimize_repeat_pdb must have integer multiple of repeat_length number of residues'
  NumberRepeats = PoseLength / RepeatLength

  # print 'NumberRepeats', NumberRepeats

  Sequence = Pose.sequence()
  # print Sequence
  
  RepeatRanges = []
  Start = 1
  for Repeat in range(NumberRepeats):
    End = Start + RepeatLength - 1
    RepeatRanges.append((Start, End))
    Start += RepeatLength

  assert len(RepeatRanges) == NumberRepeats
  
  print 'RepeatRanges', RepeatRanges

  MidRepeat = ( NumberRepeats / 2 ) - 1  
  ReferenceRange = RepeatRanges[MidRepeat]

  print 'MidRepeat', MidRepeat
  print 'ReferenceRange', ReferenceRange

  SetupNCS = symmetry.SetupNCSMover()

  for TargetRange in RepeatRanges:
    if TargetRange != ReferenceRange:
      print 'OtherRange', TargetRange

      # skip first residue (not enougth atoms for torsion)
      if TargetRange[0] == 1:
        SetupNCS.add_group( "%dA-%dA"%(ReferenceRange[0]+1, ReferenceRange[1]), "%dA-%dA"%(TargetRange[0]+1, TargetRange[1]) )        
      # skip last residue (not enougth atoms for torsion)
      elif TargetRange[1] == PoseLength:
        SetupNCS.add_group( "%dA-%dA"%(ReferenceRange[0], ReferenceRange[1]-1), "%dA-%dA"%(TargetRange[0], TargetRange[1]-1) )
      else:
        SetupNCS.add_group( "%dA-%dA"%(ReferenceRange[0], ReferenceRange[1]), "%dA-%dA"%(TargetRange[0], TargetRange[1]) )

  # default talaris 2013 score function
  ScoreFunction = rosetta.getScoreFunction()
  
  # turning on constraint weights
  ScoreFunction.set_weight(rosetta.atom_pair_constraint, 1.0)
  ScoreFunction.set_weight(rosetta.angle_constraint, 1.0)
  ScoreFunction.set_weight(rosetta.dihedral_constraint, 1.0)

  if Cst:
    # make constraint mover
    Constrainer = rosetta.ConstraintSetMover()
    # get constraints from file
    Constrainer.constraint_file(Cst)  
    # Apply constraints to pose
    Constrainer.apply(Pose)

  SymmPose = Pose.clone()
  SetupNCS.apply(SymmPose)

  # print 'Pose', Pose.constraint_set()
  # print 'SymmPose', SymmPose.constraint_set()

  # FastRelax asymmetric
  # rosetta.relax_pose(Pose, ScoreFunction, 'tag')
  # FastRelax asymmetric
  rosetta.relax_pose(SymmPose, ScoreFunction, 'tag')

  # rosetta.dump_pdb(Pose, Pdb.replace('.pdb', '_asymm.pdb') )
  rosetta.dump_pdb(SymmPose, Pdb.replace('.pdb', '_relax.pdb') )
  # # Pose


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=" args for optimize_repeat_structures ")
  ArgParser.add_argument('-pdb_stem', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  ArgParser.add_argument('-thread', type=int, help=" number of threads to run simultaneously ", default=15 ) # with default, there can be only one !    
  Args = ArgParser.parse_args()
  
  Pdbs = glob.glob('*%s.pdb'%Args.pdb_stem)
  Csts = glob.glob('*%s.cst'%Args.pdb_stem)
  Pdbs.sort()
  Csts.sort()

  if len(Csts):
    InputConstraints = True
    assert len(Pdbs) == len(Csts), 'must be equal numbers of pdb and cst files'
  else:
    InputConstraints = False

  for ThreadChunkNumber in range( (len(Pdbs)/Args.thread) + 1):
  # for ThreadChunkNumber in range( 1 ):
    Start = ThreadChunkNumber*Args.thread
    End = Start+Args.thread
    # print Start, End 
    PdbSubset = Pdbs[Start: End]
    
    if InputConstraints:
      CstSubset = Csts[Start: End]
    
    OptimizationInputTuples = []
    print 'PdbSubset:', PdbSubset 
    for i, Pdb in enumerate(PdbSubset):
      RepeatLength = int(re.sub(r'.*rep(\d+)_.*', r'\1', Pdb))
      # print RepeatLength
      if InputConstraints:
        OptimizationInputTuples.append((Pdb, CstSubset[i], RepeatLength))
      else:
        OptimizationInputTuples.append((Pdb, False, RepeatLength))

    print 'OptimizationInputTuples:', OptimizationInputTuples
    # print Pdb_Cst_Tuples

    pool = Pool(processes=len(OptimizationInputTuples))
    pool.map(optimize_repeat_pdb, OptimizationInputTuples)

    # except:
    #   for InputTuple in OptimizationInputTuples:
    #     print 'Error in one of pooled multiprocessing threads; iterating sequentially for debugging '
    #     optimize_repeat_pdb(InputTuple)

if __name__ == "__main__":
  sys.exit(main())

# def mcmc_mover(Pose, ScoreFunction, kT=1.0, Gen=100):
#   ''' #threadable backbone breathing and side chain repacking with constraints
# ''' 
#   print 'Starting work on', Pose

#   movemap = rosetta.MoveMap()
#   movemap.set_bb(True)
#   small_mover = rosetta.SmallMover(movemap, kT, 1)
#   shear_mover = rosetta.ShearMover(movemap, kT, 1)

#   MontyCarlos = rosetta.MonteCarlo(Pose, ScoreFunction, kT)

#   small_mover(Pose)

#   for g in range(Gen):
#     small_mover(Pose)
#     MontyCarlos.boltzmann(Pose)
#     shear_mover(Pose)
#     MontyCarlos.boltzmann(Pose)


#   print MontyCarlos.show_scores()
#   print MontyCarlos.show_counters()
#   print MontyCarlos.show_state()

# '''
