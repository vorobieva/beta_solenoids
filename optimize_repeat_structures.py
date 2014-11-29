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
  rosetta.init(extra_options = " -ex1 -ex2 -no_optH false -use_input_sc ") # -mute basic -mute core
  from rosetta.protocols.simple_moves import symmetry
  from rosetta.utility import ostream

  from expand_constraints import set_all_weights_zero
# from rosetta.protocols import grafting 

'''
sys.argv.extend(['-pdb_stem', '3ult_rep', '-thread', '18'])
# '''

def minimize_pose_backbone( (Pose, ScoreFunction) ):
  MinPose = Pose.clone()
  Movemap = rosetta.MoveMap()
  Movemap.set_bb(True)
  Movemap.set_chi(True)
  MinimizationMover = rosetta.MinMover()
  MinimizationMover.movemap(Movemap)
  MinimizationMover.score_function(ScoreFunction)
  MinimizationMover.apply(MinPose)
  return MinPose

def optimize_repeat_pdb( (Pdb, CstSets, RepeatLength) ):
  ''' parallelizable '''

  # idealize peptide bonds with command line subprocess
  subprocess.check_output(['idealize_jd2.default.linuxgccrelease', '-s', Pdb])
  IdealizedPdbOldName = Pdb.replace('.pdb', '_0001.pdb') 
  IdealizedPdbNewName = Pdb.replace('.pdb', '_ideal.pdb')
  subprocess.check_output(['mv', IdealizedPdbOldName, IdealizedPdbNewName])
  time.sleep(0.5)

  Pose = rosetta.pose_from_pdb(IdealizedPdbNewName)
  PoseLength = Pose.n_residue()

  assert PoseLength % RepeatLength == 0, 'pdb input into optimize_repeat_pdb must have integer multiple of repeat_length number of residues'
  NumberRepeats = PoseLength / RepeatLength

  print 'NumberRepeats', NumberRepeats
  print 'RepeatLength', RepeatLength
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

  SetupNCS.apply(Pose)

  # default talaris 2013 score function
  ScoreFunction = rosetta.getScoreFunction()

  CstFunction = set_all_weights_zero( rosetta.getScoreFunction() )
  CstFunction.set_weight(rosetta.atom_pair_constraint, 1.0)
  CstFunction.set_weight(rosetta.angle_constraint, 1000.0)
  CstFunction.set_weight(rosetta.dihedral_constraint, 100.0)
  
  for Cst in CstSets:
    print 'Cst:', Cst
    CstPose = Pose.clone()
    CstStem = (Cst+'!').replace('.cst!', '').replace('!', '')
    
    # make constraint mover
    Constrainer = rosetta.ConstraintSetMover()
    # get constraints from file
    Constrainer.constraint_file(Cst)
    Constrainer.apply(CstPose)
    CstScore0 = CstFunction(CstPose)
    TalScore0 = ScoreFunction(CstPose)

    CstMinimizedPose = minimize_pose_backbone( (CstPose, CstFunction) )
    CstScore1 = CstFunction(CstMinimizedPose)
    TalScore1 = ScoreFunction(CstMinimizedPose)
    rosetta.dump_pdb( CstMinimizedPose, CstStem+'_CstMin.pdb' )

    CstTalMinimizedPose = minimize_pose_backbone( (CstMinimizedPose, CstFunction) )
    CstScore2 = CstFunction(CstTalMinimizedPose)
    TalScore2 = ScoreFunction(CstTalMinimizedPose)
    rosetta.dump_pdb( CstTalMinimizedPose, CstStem+'_CstTalMin.pdb' )

    with open(CstStem+'_score.txt', 'w') as ScoreFile:
      print>>ScoreFile, 'Stage\t\tTalaris\t\tCst'
      print>>ScoreFile, '%d\t\t%.1f\t\t%.1f'%(0,TalScore0,CstScore0)
      print>>ScoreFile, '%d\t\t%.1f\t\t%.1f'%(1,TalScore1,CstScore1)
      print>>ScoreFile, '%d\t\t%.1f\t\t%.1f'%(2,TalScore2,CstScore2)

    # print ' Cst start pose: '
    # ScoreFunction.show(CstPose)

    # rosetta.relax_pose(CstPose, ScoreFunction, 'tag')
    # print ' Cst end pose: '
    # ScoreFunction.show(CstPose)
 
    # rosetta.dump_pdb(CstPose, Pdb.replace('.pdb', '_CstRelax.pdb') )
  

  # print ' No cst start pose: '
  # ScoreFunction.show(NoCstPose)

  # rosetta.relax_pose(NoCstPose, ScoreFunction, 'tag')

  # print ' No cst end pose: '
  # ScoreFunction.show(NoCstPose)

  # # rosetta.dump_pdb(Pose, Pdb.replace('.pdb', '_asymm.pdb') )
  # rosetta.dump_pdb(NoCstPose, Pdb.replace('.pdb', '_NoCst.pdb') )

  # ### turning on constraint weights
  # ScoreFunction.set_weight(rosetta.atom_pair_constraint, 1.0)
  # ScoreFunction.set_weight(rosetta.angle_constraint, 1000.0)
  # ScoreFunction.set_weight(rosetta.dihedral_constraint, 100.0)


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=" args for optimize_repeat_structures ")
  ArgParser.add_argument('-pdb_stem', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  ArgParser.add_argument('-thread', type=int, help=" number of threads to run simultaneously ", default=15 ) # with default, there can be only one !    
  Args = ArgParser.parse_args()
  
  Pdbs = glob.glob('*%s.pdb'%Args.pdb_stem)

  CstHash = { }
  for Pdb in Pdbs:
    CstGlober = (Pdb+'!').replace('.pdb!', '*cst')
    # print 'Pdb, CstGlober'
    # print Pdb, CstGlober
    Csts = glob.glob(CstGlober)
    Csts = [ Cst for Cst in Csts if not Cst.endswith('Temp.cst') ]
    CstHash[Pdb] = Csts

  # print CstHash
  for ThreadChunkNumber in range( (len(Pdbs)/Args.thread) + 1):
  # for ThreadChunkNumber in range( 1 ):
    Start = ThreadChunkNumber*Args.thread
    End = Start+Args.thread
    # print Start, End 
    PdbSubset = Pdbs[Start: End]
    
    OptimizationInputTuples = []
    print 'PdbSubset:', PdbSubset 
    for i, Pdb in enumerate(PdbSubset):
      RepeatLength = int(re.sub(r'.*rep(\d+)_.*', r'\1', Pdb))
      OptimizationInputTuples.append( (Pdb, CstHash[Pdb], RepeatLength) )

    pool = Pool(processes=len(OptimizationInputTuples))
    pool.map(optimize_repeat_pdb, OptimizationInputTuples)

    # # # # except:
    # for InputTuple in OptimizationInputTuples:
    #   # print 'Error in one of pooled multiprocessing threads; iterating sequentially for debugging '
    #   optimize_repeat_pdb(InputTuple)

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
