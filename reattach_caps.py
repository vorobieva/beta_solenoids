#! /usr/bin/env python
InfoString = ''' 

Run on expanded pose with src range tag in pdb name like:
.*src13_26__15_28.*.pdb

'''

'''
### external
### libraries

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
  rosetta.init(extra_options = " -ex1 -ex2 -no_optH false -use_input_sc ") # "-mute basic -mute core -mute protocols"
  from rosetta.protocols import grafting 
  
  # from repo 
  import solenoid_tools
  from generate_backbones import fuse
  from generate_backbones import get_residue_array
  from expand_constraints import constraint_extrapolator 
  from expand_constraints import pose_has

# '''

sys.argv = [ sys.argv[0], '-ref_pdb', '1EZG_Relax.pdb', '-ref_cst', '1EZG_Relax_All.cst', '-repeat_tag', 'rep24_1EZG_Relax']


def cap_and_relax_pdb( (RepeatPdb, ReferencePdb, ReferenceCst) ):

  RepeatPose = rosetta.pose_from_pdb(RepeatPdb)
  TrimmedRepeatPose = grafting.return_region( RepeatPose, 3, RepeatPose.n_residue()-3 )
  TrimmedRepeatPose.pdb_info( rosetta.core.pose.PDBInfo( TrimmedRepeatPose ) )

  ReferencePose = rosetta.pose_from_pdb( ReferencePdb )
  ReferencePose.pdb_info( rosetta.core.pose.PDBInfo( ReferencePose ) )

  # rosetta.dump_pdb(TrimmedRepeatPose, 'Trimmed.pdb')

  RepeatLength = int(re.sub(r'.*rep(\d+).*pdb', r'\1', RepeatPdb))
  SourceRanges = re.sub(r'.*src(\d+_\d+__\d+_\d+).*pdb', r'\1', RepeatPdb)
  SourceRanges = SourceRanges.split('__')
  SourceRanges = [ [ int(Value) for Value in Range.split('_') ] for Range in SourceRanges ]
  SourceStart = SourceRanges[0][0]
  SourceEnd = SourceRanges[0][1]


  '''
   Add N terminal cap 
   '''
  NcapPose = grafting.return_region( ReferencePose, 1, SourceStart+5 )
  # rosetta.dump_pdb(NcapPose, 'Ncap.pdb')
  NcapLength = NcapPose.n_residue()
  
  NcapOverhangPositions = [ Position for Position in range(NcapLength-3, NcapLength+1) ]
  # print NcapOverhangPositions
  NcapOverhangArray = get_residue_array( NcapPose, NcapOverhangPositions )
  
  RepStartOverhangPositions = [1,2,3,4]
  RepStartOverhangArray = get_residue_array( TrimmedRepeatPose, RepStartOverhangPositions )
  # print RepStartOverhangArray

  RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( NcapOverhangArray, RepStartOverhangArray )
  rosetta.Pose.apply_transform_Rx_plus_v(TrimmedRepeatPose, rMtx, tVec)
  # rosetta.dump_pdb( TrimmedRepeatPose, 'TrimmedShifted.pdb' )
  NcapPlusRepeatPose, RMSD, NcapCorrespondingResidues = fuse(NcapPose, TrimmedRepeatPose)
  print 'Ncap attachment RMSD %f'%RMSD
  # rosetta.dump_pdb( NcapPlusRepeatPose, 'NcapPlusRepeat.pdb' )
  NcapPlusRepeatPose.pdb_info( rosetta.core.pose.PDBInfo( NcapPlusRepeatPose ) )    
  

  '''
   Add C terminal cap 
  '''
  Cshift = SourceEnd-6
  CcapPose = grafting.return_region( ReferencePose, Cshift, ReferencePose.n_residue() )
  # rosetta.dump_pdb(CcapPose, 'Ccap.pdb')
  CcapOverhangPositions = [1,2,3,4]
  CcapOverhangArray = get_residue_array( CcapPose, CcapOverhangPositions )

  RepEndOverhangPositions = [ Position for Position in range( NcapPlusRepeatPose.n_residue()-3, NcapPlusRepeatPose.n_residue()+1 ) ]
  # print 'RepEndOverhangPositions', RepEndOverhangPositions
  RepEndOverhangArray = get_residue_array( NcapPlusRepeatPose, RepEndOverhangPositions )
  
  RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( RepEndOverhangArray, CcapOverhangArray )
  rosetta.Pose.apply_transform_Rx_plus_v(CcapPose, rMtx, tVec)
  # rosetta.dump_pdb( CcapPose, 'CcapPose.pdb' )
  CappedRepeatPose, RMSD, CcapCorrespondingResidues = fuse(NcapPlusRepeatPose, CcapPose)
  print 'Ccap attachment RMSD %f'%RMSD

  CappedNamePdb = re.sub(r'(.*).pdb$', r'\1_Cap.pdb', RepeatPdb)
  assert CappedNamePdb != RepeatPdb, 'regular expression substitution failed!'
  rosetta.dump_pdb( CappedRepeatPose, CappedNamePdb )



  '''
   Generate csts for cap/repeat edges 
  '''
  CstExtrapolator = constraint_extrapolator(ReferenceCst)
  ConstraintSet = []
  
  ' N cap constraints are easy; no shifts are needed '

  # For catching when individual constraints have been considered already  
  Redundict = {} 
  for Position in range(1, SourceStart+6):
    # print 'Position', Position
    # Skip positions w/out constraints
    try:
      PositionCstDict = CstExtrapolator.Cst[Position]
    except KeyError:
      continue

    for AtomName in PositionCstDict:
      for Constraint in PositionCstDict[AtomName]:
        # unpack tuple values 
        AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
        
        # Redundancy check with redundict 
        try:
          Check = Redundict[CstLineNumber]
          # if cst considered already, skip it! 
          continue
        except KeyError:
          Redundict[CstLineNumber] = 1
        
        if pose_has(CappedRepeatPose, AtomResidueCoords):
          ConstraintSet.append(Constraint)

  ' C cap constraints are harder; need to shift due to pose expansion '

  # CstExtrapolator.output_cst(ConstraintSet, 'NcapConstraints.cst')\
  Redundict = {} 

  # print 'CcapCorrespondingResidues', CcapCorrespondingResidues
  RepeatCcapPositionStart = CcapCorrespondingResidues[0][0]
  # print 'RepeatCcapPositionStart', RepeatCcapPositionStart

  ShiftToRepeatPose = RepeatCcapPositionStart - Cshift
  # print 'ShiftToRepeatPose', ShiftToRepeatPose

  for Position in range( Cshift, ReferencePose.n_residue()+1 ):
    # Skip positions w/out constraints
    try:
      PositionCstDict = CstExtrapolator.Cst[Position]
    except KeyError:
      continue

    for AtomName in PositionCstDict:
      for Constraint in PositionCstDict[AtomName]:
        # unpack tuple values 
        AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
        
        # Redundancy check with redundict 
        try:
          Check = Redundict[CstLineNumber]
          # if cst considered already, skip it! 
          continue
        except KeyError:
          Redundict[CstLineNumber] = 1

        ExpandedPoseAtomResidueCoords = []
        # iterate through atom residue pairs
        for AtomResiduePair in AtomResidueCoords:
          # print 'AtomResiduePair', AtomResiduePair
          ExpandedPosePosition = (AtomResiduePair[1]) + ShiftToRepeatPose
          # print 'ExpandedPosePosition', ExpandedPosePosition
          ExpandedPoseAtomResidueCoords.append( ( AtomResiduePair[0], ExpandedPosePosition ) )

        ShiftedConstraint = ExpandedPoseAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType       

        if pose_has(CappedRepeatPose, ExpandedPoseAtomResidueCoords):
          ConstraintSet.append(ShiftedConstraint)  

  CapCstName = re.sub(r'(.*).pdb$', r'\1.cst', CappedNamePdb)
  CstExtrapolator.output_cst(ConstraintSet, CapCstName)

  '''
  idealize peptide bonds with command line subprocess
  '''
  subprocess.check_output(['idealize_jd2.default.linuxgccrelease', '-s', CappedNamePdb])
  IdealizedPdbOldName = re.sub(r'(.*).pdb$', r'\1_0001.pdb', CappedNamePdb)
  IdealizedPdbNewName = re.sub(r'(.*).pdb$', r'\1_Ideal.pdb', CappedNamePdb)
  
  subprocess.check_output(['mv', IdealizedPdbOldName, IdealizedPdbNewName])
  time.sleep(0.2)

  IdealizedCappedPose = rosetta.pose_from_pdb( IdealizedPdbNewName )

  # make constraint mover
  Constrainer = rosetta.ConstraintSetMover()
  # get constraints from file
  Constrainer.constraint_file(CapCstName)
  Constrainer.apply(IdealizedCappedPose)


  ''' SET UP WEIGHTS AS decided '''


  Talaris = rosetta.getScoreFunction()
  rosetta.relax_pose(IdealizedCappedPose, Talaris, 'tag')


  RelaxedPdbName = re.sub(r'(.*)_Ideal.pdb$', r'\1_Relax.pdb', IdealizedPdbNewName)
  rosetta.dump_pdb(IdealizedCappedPose, RelaxedPdbName)



def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' reattach_cap.py ( -help ) %s'%InfoString)
  # Required args
  ArgParser.add_argument('-ref_pdb', type=str, help=' Reference pdb ', required=True)
  ArgParser.add_argument('-ref_cst', type=str, help=' Reference pdb ', required=True)
  ArgParser.add_argument('-repeat_tag', type=str, help=' Input pdb tag ', required=True)
  ArgParser.add_argument('-thread', type=int, help=" number of threads to run simultaneously ", default=15 ) # with default, there can be only one !

  # Optional args
  ArgParser.add_argument('-out', type=str, help=' Output directory ', default='./')
  Args = ArgParser.parse_args()
  
  ReferencePdb = Args.ref_pdb
  ReferenceCst = Args.ref_cst

  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  GlobString = '*%s*.pdb'%Args.repeat_tag
  # print GlobString  
  Pdbs = glob.glob( GlobString ) 
  
  Pdbs = [ Pdb for Pdb in Pdbs if not 'Cap' in Pdb ]

  print ' Globbed %d pdb(s) with: %s '%( len(Pdbs), GlobString)

  for ThreadChunkNumber in range( (len(Pdbs)/Args.thread) + 1):

    Start = ThreadChunkNumber*Args.thread
    End = Start+Args.thread
    # print Start, End 
    PdbSubset = Pdbs[Start: End]
    
    # pass strings, not poses, for ez pickling
    ParallelizableInputTuples = []
    for Pdb in PdbSubset:
      ParallelizableInputTuples.append((Pdb, ReferencePdb, ReferenceCst))
    
    for InputTuple in ParallelizableInputTuples:
      # print 'Error in one of pooled multiprocessing threads; iterating sequentially for debugging '
      cap_and_relax_pdb(InputTuple)

    # pool = Pool(processes=len(ParallelizableInputTuples))
    # pool.map(cap_and_relax_pdb, ParallelizableInputTuples)



# if __name__ == "__main__":
#   sys.exit(main())

