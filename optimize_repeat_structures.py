#! /usr/bin/env python
InfoString = ''' 
To optimize repeat backbone with NCS torsion constraint and other constraints'''


# '''
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

# from rosetta.protocols import grafting 

'''
### '''


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

def responsibly_optimize_repeat_pdbs( ListOfPdbTuplesForThisThread ):

  for OptTuple in ListOfPdbTuplesForThisThread:
    print 'Starting new optimization trajectory'
    optimize_repeat_pdb(OptTuple)


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

  # print 'NumberRepeats', NumberRepeats
  # print 'RepeatLength', RepeatLength
  Sequence = Pose.sequence()
  # print Sequence
  
  RepeatRanges = []
  Start = 1
  for Repeat in range(NumberRepeats):
    End = Start + RepeatLength - 1
    RepeatRanges.append((Start, End))
    Start += RepeatLength

  assert len(RepeatRanges) == NumberRepeats
  # print 'RepeatRanges', RepeatRanges

  MidRepeat = ( NumberRepeats / 2 ) - 1  
  ReferenceRange = RepeatRanges[MidRepeat]
  # print 'MidRepeat', MidRepeat
  # print 'ReferenceRange', ReferenceRange

  SetupNCS = symmetry.SetupNCSMover()

  for TargetRange in RepeatRanges:
    if TargetRange != ReferenceRange:
      # print 'OtherRange', TargetRange
      # skip first three residue (not enougth atoms for torsion), and amino acid types allowed to vary
      if TargetRange[0] == 1:
        SetupNCS.add_group( "%dA-%dA"%(ReferenceRange[0]+3, ReferenceRange[1]), "%dA-%dA"%(TargetRange[0]+3, TargetRange[1]) )        
      # skip last residue (not enougth atoms for torsion)
      elif TargetRange[1] == PoseLength:
        SetupNCS.add_group( "%dA-%dA"%(ReferenceRange[0], ReferenceRange[1]-3), "%dA-%dA"%(TargetRange[0], TargetRange[1]-3) )
      else:
        SetupNCS.add_group( "%dA-%dA"%(ReferenceRange[0], ReferenceRange[1]), "%dA-%dA"%(TargetRange[0], TargetRange[1]) )

  SetupNCS.apply(Pose)

  # default talaris 2013 score function plus dihedral wieght for symmetry ncs mimization
  SymmTalaris = rosetta.getScoreFunction()
  SymmTalaris.set_weight(rosetta.dihedral_constraint, 1.0)

  TalarisPlusCst = rosetta.getScoreFunction()
  TalarisPlusCst.set_weight(rosetta.atom_pair_constraint, 10.0)
  TalarisPlusCst.set_weight(rosetta.angle_constraint, 5.0)
  TalarisPlusCst.set_weight(rosetta.dihedral_constraint, 3.0)

  TalarisPlusCstLowerFaRep = rosetta.getScoreFunction()
  TalarisPlusCstLowerFaRep.set_weight(rosetta.atom_pair_constraint, 10.0)
  TalarisPlusCstLowerFaRep.set_weight(rosetta.angle_constraint, 5.0)
  TalarisPlusCstLowerFaRep.set_weight(rosetta.dihedral_constraint, 3.0)
  TalarisPlusCstLowerFaRep.set_weight(rosetta.fa_rep, 0.25)
  print 'Pdb:', Pdb

  OptimizedPoses = []
  PoseIDs = []

  for Cst in CstSets:
    print 'Cst:', Cst
    CstPose = Pose.clone()
    CstStemName = re.sub(r'^(.*)\.cst$', r'\1', Cst)

    # make constraint mover
    Constrainer = rosetta.ConstraintSetMover()
    # get constraints from file
    Constrainer.constraint_file(Cst)
    Constrainer.apply(CstPose)

    FxnTags = [ 'TalCst', 'LowFaRep'  ]

    for i, ScoreFunction in enumerate( [ TalarisPlusCst, TalarisPlusCstLowerFaRep ] ):
      # for AbsoluteWeight in [1, 5, 10, 100]:

      RelaxPose = CstPose.clone()
      rosetta.relax_pose(RelaxPose, ScoreFunction, 'tag')
      rosetta.dump_pdb( RelaxPose, CstStemName+'_%s.pdb'%FxnTags[i] )
      # remove all constraints
      RelaxPose.remove_constraints()
      # reapply ncs constraints
      SetupNCS.apply(RelaxPose)

      rosetta.relax_pose(RelaxPose, SymmTalaris, 'tag')
      # Trekker.score(RelaxPose)
      rosetta.dump_pdb( RelaxPose, CstStemName+'_%s_Relax.pdb'%FxnTags[i] )

  JustRelaxPose = Pose.clone()
  SetupNCS.apply( JustRelaxPose )

  rosetta.relax_pose( JustRelaxPose, SymmTalaris, 'tag' )
  rosetta.dump_pdb( JustRelaxPose, CstStemName+'_JustRelax.pdb' )

#### sys.argv = [sys.argv[0], '-pdb_stem', 'rep40_2G0Y_Relax', '-thread', '5']

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

  JobOptTuplesSortedByThread = [[] for i in range(Args.thread)]

  print 'JobOptTuplesSortedByThread', JobOptTuplesSortedByThread
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

      # OptimizationInputTuples.append( (Pdb, CstHash[Pdb], RepeatLength) )
      JobOptTuplesSortedByThread[i].append( (Pdb, CstHash[Pdb], RepeatLength) )

    # for InputTuple in OptimizationInputTuples:
    #   # print 'Error in one of pooled multiprocessing threads; iterating sequentially for debugging '
    #   optimize_repeat_pdb(InputTuple)

    # pool = Pool(processes=len(OptimizationInputTuples))
    # pool.map(optimize_repeat_pdb, OptimizationInputTuples)


  # for JobSet in JobOptTuplesSortedByThread:
  #   responsibly_optimize_repeat_pdbs(JobSet)

  pool = Pool(processes=len(JobOptTuplesSortedByThread))
  pool.map(responsibly_optimize_repeat_pdbs, JobOptTuplesSortedByThread)
  


if __name__ == "__main__":
  sys.exit(main())

