#! /usr/bin/env python
InfoString = ''' 

Run on expanded pose with src range tag in pdb name like:
.*src13_26__15_28.*.pdb

'''

# '''
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
  rosetta.init(extra_options = " -ex1 -ex2 -no_optH false -use_input_sc -mute basic -mute core -mute protocols ") 
  from rosetta.protocols import grafting 
  
  # from repo 
  import solenoid_tools
  import generate_backbones
  import generate_cst
  import expand_cst

# '''


def cap_pdb_make_cst( RepeatPdbFileName, RepeatCstFileName, ReferencePdb, ReferenceCst, Ntrim=0, Ctrim=0, Step=0 ):
  if Step:
    if Ntrim:
      assert Ntrim % Step == 0
    if Ctrim:
      assert Ctrim % Step == 0 

  # Grep out repeat length and src ranges 
  RepeatLength = int(re.sub(r'.*rep(\d+).*pdb', r'\1', RepeatPdbFileName))
  SourceRanges = re.sub(r'.*src(\d+_\d+__\d+_\d+).*pdb', r'\1', RepeatPdbFileName)
  assert SourceRanges != RepeatPdbFileName, 'src string not found in pdb name '
  SourceRanges = SourceRanges.split('__')
  SourceRanges = [ [ int(Value) for Value in Range.split('_') ] for Range in SourceRanges ]
  SourceStart = SourceRanges[0][0]
  SourceEnd = SourceRanges[0][1]

  # Load repeat pose
  RepeatPose = rosetta.pose_from_pdb( RepeatPdbFileName )
  # Trim off floppy end residues
  TrimmedRepeatPose = grafting.return_region( RepeatPose, 3, RepeatPose.n_residue()-3 )
  TrimmedRepeatPose.pdb_info( rosetta.core.pose.PDBInfo( TrimmedRepeatPose ) )
  # rosetta.dump_pdb(TrimmedRepeatPose, 'Trimmed.pdb')
  # Load reference (native) pose
  ReferencePose = rosetta.pose_from_pdb( ReferencePdb )
  ReferencePose.pdb_info( rosetta.core.pose.PDBInfo( ReferencePose ) )

  PdbCstPairs = []

  ''' Loop through N terminal caps '''
  # print '(SourceStart-Ntrim, SourceStart, -1*Step)', (SourceStart-Ntrim, SourceStart, -1*Step)

  for NcapTrimBackSteps in range(0, (Ntrim/Step) + 1 ):
    # print 'Ntrils -ltrhm:', NcapTrimBackSteps * Step
    NcapLastRes = SourceStart - (NcapTrimBackSteps * Step)
    # print 'NcapLastRes:', NcapLastRes

    ### Get pose for n-terminal cap with overhang for superimpositions
    try:
      NcapPose = grafting.return_region( ReferencePose, 1, NcapLastRes+5 )
    except RuntimeError:
      print 'Requested end of n-terminal cap, %d, beyond range of reference protein. '%NcapLastRes
      continue    
    except OverflowError:
      print 'Requested end of n-terminal cap, %d, beyond range of reference protein. '%NcapLastRes
      continue    

    try:
      assert NcapPose.n_residue() > 4
    except AssertionError:
      print 'Too few residues to attach n-terminal cap ending at %d; skipping '%NcapLastRes
      continue

    # rosetta.dump_pdb(NcapPose, 'Ncap.pdb')
    NcapLength = NcapPose.n_residue()
    
    NcapOverhangPositions = [ Position for Position in range(NcapLength-3, NcapLength+1) ]
    # print NcapOverhangPositions
    NcapOverhangArray = generate_backbones.get_residue_array( NcapPose, NcapOverhangPositions )
    

    RepStartOverhangPositions = [1, 2, 3, 4]
    RepStartOverhangArray = generate_backbones.get_residue_array( TrimmedRepeatPose, RepStartOverhangPositions )
    # print RepStartOverhangArray

    RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( NcapOverhangArray, RepStartOverhangArray )
    rosetta.Pose.apply_transform_Rx_plus_v(TrimmedRepeatPose, rMtx, tVec)
    # rosetta.dump_pdb( TrimmedRepeatPose, 'TrimmedShifted.pdb' )
    
    try:
      NcapPlusRepeatPose, RMSD, NcapCorrespondingResidues = generate_backbones.fuse(NcapPose, TrimmedRepeatPose)
    except AssertionError:
      print 'Not enough structural similarity to attach n-terminal cap ending at %d; skipping '%NcapLastRes
      continue

    # print 'Ncap attachment RMSD %f'%RMSD
    # rosetta.dump_pdb( NcapPlusRepeatPose, 'NcapPlusRepeat.pdb' )
    NcapPlusRepeatPose.pdb_info( rosetta.core.pose.PDBInfo( NcapPlusRepeatPose ) )    

    RepeatCstExtrapolator = expand_cst.constraint_extrapolator(RepeatCstFileName)
    # print 'NcapLastRes', NcapLastRes
    # print NcapPlusRepeatPose

    ''' Shift repeat unit constraints to accomadiate numbering with n-cap length'''
    Redundict = {}
    RepeatCsts = []
    for RepeatPosition in range(1, RepeatPose.n_residue()+1 ):
      # print 'RepeatPosition', RepeatPosition
      try:
        RepeatPositionCstDict = RepeatCstExtrapolator.Cst[RepeatPosition]
      except KeyError:
        continue
      for AtomName in RepeatPositionCstDict:
        for Cst in RepeatPositionCstDict[AtomName]:
          ### unpack tuple values 
          AtomResidueCoords, CstParameters, CstLineNumber, CstType = Cst
          ### Redundancy check with redundict 
          try:
            Check = Redundict[CstLineNumber]
            ### if cst considered already, skip it! 
            continue
          except KeyError:
            Redundict[CstLineNumber] = 1

          ShiftedPoseAtomResidueCoords = []
          ### iterate through atom residue pairs
          for AtomResiduePair in AtomResidueCoords:
            # print 'AtomResiduePair', AtomResiduePair
            RepeatPosePosition = (AtomResiduePair[1]) + NcapLastRes - 1
            # print 'RepeatPosePosition', RepeatPosePosition
            ShiftedPoseAtomResidueCoords.append( ( AtomResiduePair[0], RepeatPosePosition ) )

          ShiftedCst = ShiftedPoseAtomResidueCoords, CstParameters, CstLineNumber, CstType       

          
          if expand_cst.pose_has(NcapPlusRepeatPose, ShiftedPoseAtomResidueCoords):
            RepeatCsts.append(ShiftedCst)
          try:
            assert expand_cst.pose_has(NcapPlusRepeatPose, ShiftedPoseAtomResidueCoords), ' Cst shifted from repeat pose not found in capped pose'
          except AssertionError:
            pass
            # print 'AtomResidueCoords', AtomResidueCoords
            # print 'ShiftedPoseAtomResidueCoords', ShiftedPoseAtomResidueCoords

  
    ''' Loop through C terminal caps '''
    for CcapTrimForwardSteps in range(0, (Ctrim/Step) + 1 ):
      # print 'CcapTrimForwardSteps', CcapTrimForwardSteps
      CcapFirstRes = SourceEnd + ( CcapTrimForwardSteps * Step )
      # print 'CcapFirstRes:', CcapFirstRes
      Cshift = CcapFirstRes-6
      print 'Cshift', Cshift
      print 'ReferencePose.n_residue()', ReferencePose.n_residue()
      
      try:
        CcapPose = grafting.return_region( ReferencePose, Cshift, ReferencePose.n_residue() )
      except RuntimeError:
        print 'Requested start of c-terminal, %d, beyond range of reference protein. '%CcapFirstRes
        continue        
      except OverflowError:
        print 'Requested start of c-terminal, %d, beyond range of reference protein. '%CcapFirstRes
        continue   

      # rosetta.dump_pdb(CcapPose, 'Ccap.pdb')

      try:
        assert CcapPose.n_residue() > 4
      except AssertionError:
        print 'Too few residues to attach c-terminal cap starting at %d; skipping '%CcapFirstRes
        continue

      CcapOverhangPositions = [1, 2, 3, 4]
      CcapOverhangArray = generate_backbones.get_residue_array( CcapPose, CcapOverhangPositions )

      RepEndOverhangPositions = [ Position for Position in range( NcapPlusRepeatPose.n_residue()-3, NcapPlusRepeatPose.n_residue()+1 ) ]
      # print 'RepEndOverhangPositions', RepEndOverhangPositions
      RepEndOverhangArray = generate_backbones.get_residue_array( NcapPlusRepeatPose, RepEndOverhangPositions )
      
      RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( RepEndOverhangArray, CcapOverhangArray )
      rosetta.Pose.apply_transform_Rx_plus_v(CcapPose, rMtx, tVec)
      # rosetta.dump_pdb( CcapPose, 'CcapPose.pdb' )

      try:
        CappedRepeatPose, RMSD, CcapCorrespondingResidues = generate_backbones.fuse(NcapPlusRepeatPose, CcapPose)
      except AssertionError:
        print 'Not enough structural similarity to attach c-terminal cap starting at %d; skipping '%CcapFirstRes
        continue

      CappedNamePdb = re.sub(r'(.*).pdb$', r'\1_%dCap%d.pdb'%(NcapLastRes, CcapFirstRes), RepeatPdbFileName)
      assert CappedNamePdb != RepeatPdbFileName, 'regular expression substitution failed!'
      
      rosetta.dump_pdb( CappedRepeatPose, CappedNamePdb )

      ''' Generate csts for cap/repeat edges '''
      CapCstExtrapolator = expand_cst.constraint_extrapolator(ReferenceCst)
      CapCsts = []
      
      ' N cap constraints are easy; no shifts are needed '

      # For catching when individual constraints have been considered already  
      Redundict = {} 
      for Position in range(1, NcapLastRes):
        # print 'Position', Position
        # Skip positions w/out constraints
        try:
          PositionCstDict = CapCstExtrapolator.Cst[Position]
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
            
            CapCsts.append(Constraint)

      ' C cap constraints are harder; need to shift due to pose expansion '
      CcapCstShift = CappedRepeatPose.n_residue() - ReferencePose.n_residue()

      # CapCstExtrapolator.output_cst(CapCsts, 'NcapConstraints.cst')\
      Redundict = {} 

      # print 'CcapCorrespondingResidues', CcapCorrespondingResidues
      RepeatCcapPositionStart = CcapCorrespondingResidues[0][0]
      # print 'RepeatCcapPositionStart', RepeatCcapPositionStart

      ShiftToRepeatPose = RepeatCcapPositionStart - Cshift
      # print 'ShiftToRepeatPose', ShiftToRepeatPose

      for Position in range( CcapFirstRes, ReferencePose.n_residue()+1 ):
        # Skip positions w/out constraints
        try:
          PositionCstDict = CapCstExtrapolator.Cst[Position]
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
              ExpandedPosePosition = (AtomResiduePair[1]) + CcapCstShift  
              # print 'ExpandedPosePosition', ExpandedPosePosition
              ExpandedPoseAtomResidueCoords.append( ( AtomResiduePair[0], ExpandedPosePosition ) )

            ShiftedConstraint = ExpandedPoseAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType       

            CapCsts.append(ShiftedConstraint)  


      CappedCstName = re.sub(r'(.*).pdb$', r'\1.cst', CappedNamePdb)
      
      with open(CappedCstName, 'w') as OverwriteExistingFile:
        pass
      
      FinalCstSet = []
      
      for Cst in CapCsts:
        if expand_cst.pose_has(CappedRepeatPose, Cst[0]):
          FinalCstSet.append(Cst)
      for Cst in RepeatCsts:
        if expand_cst.pose_has(CappedRepeatPose, Cst[0]):
          FinalCstSet.append(Cst)

      CapCstExtrapolator.output_cst(FinalCstSet, CappedCstName)        
      PdbCstPairs.append((CappedNamePdb, CappedCstName))
  
  return PdbCstPairs


def idealize_and_relax_pdb_set( PdbCstPairs ):
  
  for PdbName, CstName in PdbCstPairs:
    print '(PdbName, CstName) ', (PdbName, CstName) 
    '''  idealize peptide bonds with command line subprocess  '''
    subprocess.check_output([ 'idealize_jd2.default.linuxgccrelease', '-s', PdbName ])
    IdealizedPdbOldName = re.sub( r'(.*).pdb$', r'\1_0001.pdb', PdbName )
    IdealizedPdbNewName = re.sub( r'(.*).pdb$', r'\1_Ideal.pdb', PdbName )
    
    subprocess.check_output(['mv', IdealizedPdbOldName, IdealizedPdbNewName])
    time.sleep(0.2)

    IdealizedCappedPose = rosetta.pose_from_pdb( IdealizedPdbNewName )

    # make constraint mover
    Constrainer = rosetta.ConstraintSetMover()
    # get constraints from file
    Constrainer.constraint_file(CstName)
    Constrainer.apply(IdealizedCappedPose)


    ''' SET UP WEIGHTS '''
    Talaris = rosetta.getScoreFunction()
    TalarisPlusCst = rosetta.getScoreFunction()

    TalarisPlusCst.set_weight(rosetta.atom_pair_constraint, 10.0)
    TalarisPlusCst.set_weight(rosetta.angle_constraint, 5.0)
    TalarisPlusCst.set_weight(rosetta.dihedral_constraint, 2.5)
    print 'relaxing %s with %s'%(IdealizedPdbNewName, CstName) 

    # relax w/ cst
    rosetta.relax_pose(IdealizedCappedPose, TalarisPlusCst, 'tag')
    # relax w/o cst
    rosetta.relax_pose(IdealizedCappedPose, Talaris, 'tag')

    RelaxedPdbName = re.sub(r'(.*)_Ideal.pdb$', r'\1_Relax.pdb', IdealizedPdbNewName)
    rosetta.dump_pdb(IdealizedCappedPose, RelaxedPdbName)

    # yield IdealizedCappedPose


def relax_and_score():
  MinimizedPose = idealize_and_relax_pdb_set(PdbName, CstName)
  LogName = re.sub(r'^(.*).pdb$', r'\1.log', PdbName)
  # score_pose(MinimizedPose, LogName)

# def score_pose(Pose):
#   AtomPairCst = expand_cst.set_all_weights_zero( rosetta.getScoreFunction() )
#   AtomPairCst.set_weight(rosetta.atom_pair_constraint, 1.0)
  # Scores = {}
  # Scores["RosettaScore"] = Talaris(IdealizedCappedPose) 
  # Scores["AtomPair"] = AtomPairCst(IdealizedCappedPose) 


# sys.argv = [ sys.argv[0], '-ref_pdb', '2G0Y_Relax.pdb', '-ref_cst', '2G0Y_Relax_All.cst', '-repeat_tag', '_src22_41__35_54_rep40_2rwAsp_hfGlu']
# sys.argv[4] = '2G0Y_Relax_Disulf.cst'

def main(argv = None):
  if argv != None:
    print 'sys.argv:', sys.argv
    sys.argv = [sys.argv[0]] + [ arg for arg in argv ]

  ArgParser = argparse.ArgumentParser(description=' reattach_cap.py ( -help ) %s '%InfoString)
  # Required args
  ArgParser.add_argument('-ref_pdb', type=str, help=' Reference pdb ', required=True)
  ArgParser.add_argument('-ref_cst', type=str, help=' Reference pdb ', required=True)
  ArgParser.add_argument('-repeat_tag', type=str, help=" Input pdb tag, used like:  glob.glob('*%s*.pdb'%Args.repeat_tag ", required=True )
  # Optional args
  ArgParser.add_argument('-thread', type=int, help=" number of threads to run simultaneously ", default=1 ) # with default, there can be only one !
  ArgParser.add_argument('-n_trim', type=int, help=" max number of residues to trim off of end n terminal cap ", default=0 ) 
  ArgParser.add_argument('-c_trim', type=int, help=" max number of residues to trim off begining of c terminal cap ", default=0 ) 
  ArgParser.add_argument('-step', type=int, help=" number of residues to step back for each trim ", default=1 )
  # ArgParser.add_argument('-out', type=str, help=' Output directory ', default='./')
  Args = ArgParser.parse_args()
  
  ReferencePdb = Args.ref_pdb
  ReferenceCst = Args.ref_cst
 
  # if Args.out [-1] != '/':
  #   Args.out = Args.out + '/'

  GlobString = '*%s*.pdb'%Args.repeat_tag
  # print GlobString  
  Pdbs = glob.glob( GlobString ) 
  Pdbs = [ Pdb for Pdb in Pdbs if not 'Cap' in Pdb ]

  print ' Globbed %d pdb(s) with: %s '%( len(Pdbs), GlobString )

  InputTuples = []
  for Pdb in Pdbs:
    # print 'Pdb', Pdb

    Cst = re.sub(r'^(.*).pdb$', r'\1_All.cst', Pdb)
    # print 'Cst', Cst
    
    try:
      with open(Cst, 'r') as CheckingIfExists:
        pass 
    except IOError:
      generate_cst.main(['-pdbs', Pdb])

    CappedPdbCstPairs = cap_pdb_make_cst( Pdb, Cst, ReferencePdb, ReferenceCst, Args.n_trim, Args.c_trim, Args.step )
    # print 'CappedPdbCstPairs:', CappedPdbCstPairs
    for CapPdb, CapCst in CappedPdbCstPairs:
      InputTuples.append( (CapPdb, CapCst) )

  # print 'InputTuples', InputTuples

  ParallelInputTuples = [ [] for Thread in range(Args.thread)]
  # print '\n'*3
  # print 'len(InputTuples)', len(InputTuples)
  LengthOfParallelList = len(InputTuples)/Args.thread
  if len(InputTuples)%Args.thread != 0:
    LengthOfParallelList += 1
  # print 'LengthOfParallelList', LengthOfParallelList

  for ThreadChunkNumber in range( Args.thread ):
    Start = ThreadChunkNumber*LengthOfParallelList
    End = Start+LengthOfParallelList
    
    ParallelInputTuples[ThreadChunkNumber] = InputTuples[Start: End]
 
  # print 'ParallelInputTuples', ParallelInputTuples
  # print 'InputTupleSubset', InputTupleSubset
  ## Iterative for debug    
  # for InputTuple in ParallelInputTuples:
  #   # print 'Error in one of pooled multiprocessing threads; iterating sequentially for debugging '
  #   idealize_and_relax_pdb_set(InputTuple)

  ### Parallel for running
  pool = Pool(processes=len(ParallelInputTuples))
  pool.map(idealize_and_relax_pdb_set, ParallelInputTuples)

if __name__ == "__main__":
  sys.exit(main())
