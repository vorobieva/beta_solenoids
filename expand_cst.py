#! /usr/bin/env python
InfoString = ''' 
expand_cst.py -ref_pdb 3fy3.pdb -ref_cst 3fy3_All.cst -repeat_pdb_tag _3fy3

!!! Be careful and considerate with this script as runs things in parrallel !!!
optimize_repeat_structures.py -pdb_stem _3fy3 -thread 10
OR
!!! Now your also generating silent processes !!!
nohup optimize_repeat_structures.py -pdb_stem _3fy3 -thread 10 > log.txt &
'''

# '''
# libraries
# from scipy import spatial
# import itertools
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
  rosetta.init()
  # from beta_solenoids repo 
  import solenoid_tools
  ## debug
  # rosetta.init(extra_options = "-out:level 1000")

  ## runtime
  # rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
  from rosetta.protocols import grafting 
  from rosetta.utility import ostream

# ''' 
# sys.argv = [sys.argv[0], '-repeat_pdb_tag', '_1EZG', '-ref_pdb', '1EZG.pdb', '-ref_cst', '1EZG_ON_All.cst']

class constraint_extrapolator:
  """ constraint extrapolator """
  def __init__(self, CstFilename):
    ''' Makes hash: self.Cst[Residue] = { AtomName: (AllAtomResiduePairs, ConstraintParameters, LineNumber) } '''

    self.CstFilename = CstFilename
    with open(CstFilename, 'r') as CstFile:
      self.CstLines = CstFile.readlines()
    
    # print 'self.CstLines:', self.CstLines
    self.Cst = {}

    for iLine, Line in enumerate(self.CstLines):
      Line = Line.replace('\n','')
      LineList = Line.split()
      if Line.startswith('AtomPair'):
        CstType, Name1, Res1, Name2, Res2 = tuple( LineList[:5] )
        AtomResidueCoords = [ (Name1, int(Res1) ), (Name2, int(Res2) ) ]
        ConstraintParameters = LineList[5:]

      elif Line.startswith('Angle'):
        CstType, Name1, Res1, Name2, Res2, Name3, Res3 = tuple( LineList[:7] )
        AtomResidueCoords = [ (Name1, int(Res1) ), (Name2, int(Res2) ), (Name3, int(Res3) ) ]
        ConstraintParameters = LineList[7:]

      elif Line.startswith('Dihedral'):        
        CstType, Name1, Res1, Name2, Res2, Name3, Res3, Name4, Res4  = tuple( LineList[:9] )
        AtomResidueCoords = [ (Name1, int(Res1) ), (Name2, int(Res2) ), (Name3, int(Res3) ), (Name4, int(Res4) ) ]
        ConstraintParameters = LineList[9:]
      
      elif Line.startswith('#'):  
        pass
      
      else:
        assert len(Line) < 1, ' Cst line... \n%s\n ...not recognized '%Line

      for Atom, Residue in AtomResidueCoords:
        # these check for residue entry in constraint dict
        try:
          check = self.Cst[Residue]
        except KeyError:
          self.Cst[Residue] = {}

        # PartnerAtomResidues = []
        # for OtherAtom, OtherResidue in AtomResidueCoords:
        #   if OtherAtom == Atom and OtherResidue == Residue:
        #     continue
        #   else:
        #     PartnerAtomResidues.append( (OtherAtom, OtherResidue) )

        # add atom name entries for residue subdict as needed
        try:
          self.Cst[Residue][Atom].append( (AtomResidueCoords, ConstraintParameters, iLine+1, LineList[0]) )
        except KeyError:
          self.Cst[Residue][Atom] = [ (AtomResidueCoords, ConstraintParameters, iLine+1, LineList[0]) ] 

      # return AtomResidueCoords
  def output_cst(self, ConstraintTuples, FileName):
    with open(FileName, 'w') as CstFile:  
      for CstTuple in ConstraintTuples:
        CstString = self.reassemble_cst( CstTuple[3], CstTuple[0], CstTuple[1] )
        print>>CstFile, CstString 
      # pass

  def reassemble_cst(self, ConstraintType, AtomNameResNumberPairs, ConstraintParameters):
    # print 'AtomNameResNumberPairs', AtomNameResNumberPairs

    # assert len(AtomNameResNumberPairs) % 2 == 0, ' Atom name / residue number list must be equal length ! '
    assert len(AtomNameResNumberPairs) <= 4, ' more than 4 Atom name / residue pairs '

    AtomResidueStringFormatingList = []
    for AtomResTuple in AtomNameResNumberPairs:
      for Item in AtomResTuple:
      	AtomResidueStringFormatingList.append( str(Item) )
    
    AtomResidueString = ' '.join( AtomResidueStringFormatingList )
    return '%s %s %s'%(ConstraintType, AtomResidueString , ' '.join(ConstraintParameters) )

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

  def shift_and_sort_constraints(self, ReferenceStart, ReferenceEnd, RepeatUnitLength):  
    # One of two edge sets will be chosen based on the ability of existing amino acids to host the interactions
    # and the scores of the constraints on the extrapolated pose
    
    # For catching when individual constraints have been considered already  
    Redundict = {}

    # Occur between first turn of repeat and upstream residues
    Edge1Cst = []
    # Occur between last turn of repeat and downstream residues 
    Edge2Cst = []
    # Involve residues in repeat, and both upstream and downstream residues 
    BothEdges = []
    # Occur in middle of repeat, will be included either way
    MiddleCst = []

    ShiftToRepeat = ReferenceStart - 1

    for ExtrapolationIndex, ReferencePosition in enumerate( range(ReferenceStart, ReferenceEnd+1) ):
      ExtrapolationPosition = ExtrapolationIndex + 1
      # Skip positions w/out constraints
      try:
        # print 'trying: ', ReferencePosition
        ReferencePositionCstDict = self.Cst[ReferencePosition]
      except KeyError:
        continue

      for AtomName in ReferencePositionCstDict:
        # print 'AtomName: ', AtomName
        # Loop through constraint tuples in hash
        for Constraint in ReferencePositionCstDict[AtomName]:

          # unpack tuple values 
          AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
          
          # Redundancy check with redundict 
          try:
            Check = Redundict[CstLineNumber]
            # if cst considered already, skip it! 
            continue
          except KeyError:
            # to prevent examining individual constraint again in subsequent iterations
            Redundict[CstLineNumber] = 1

          # print 'ReferenceStart, ReferenceEnd', ReferenceStart, ReferenceEnd
          # print 'ReferencePosition:', ReferencePosition
          # print 'ExtrapolationPosition:', ExtrapolationPosition
          # print 'Constraint:', Constraint
          # print 'AtomResidueCoords:', AtomResidueCoords

          # ResidueNumber are shifted to correspond with first repeat unit in the extrapolated pose
          RepeatAtomResidueCoords = []
          E1 = 0
          E2 = 0
          Mid = 0

          # iterate thorugh atom reisdue pairs
          for AtomResiduePair in AtomResidueCoords:
            # print 'AtomResiduePair', AtomResiduePair
            RepeatPosition = (AtomResiduePair[1]) - ShiftToRepeat
            # print 'RepeatPosition', RepeatPosition
            RepeatAtomResidueCoords.append( ( AtomResiduePair[0], RepeatPosition ) )
            
            if RepeatPosition < 1:
              E1 += 1
            elif RepeatPosition > RepeatUnitLength:
              E2 += 1
            else:
              Mid += 1

          ShiftedConstraint = RepeatAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType
          # print 'E1', E1
          # print 'E2', E2
          # print 'Mid', Mid
          # print          
          
          if Mid:
            if E1 and E2:
              BothEdges.append(ShiftedConstraint)
            elif E1:
              Edge1Cst.append(ShiftedConstraint)
            elif E2:
              Edge2Cst.append(ShiftedConstraint)
            else:
              MiddleCst.append(ShiftedConstraint)
            # print 'Edge1 and Edge2'
            # sys.exit()
    
    return ( Edge1Cst, Edge2Cst, BothEdges, MiddleCst )

  def parse_ostream_cst(self, String):
    ModString = String.replace('Angle', '!Angle').replace('Dihedral', '!Dihedral').replace('AtomPair', '!AtomPair') 
    List = ModString.split('!')
    ModList = [ String for String in List if len(String) > 2 ]
    return ModList
  # def inbound(self, UpEdge, DownEdge, AtomResidueCoords):
    # AtomResidueCoords = [ UpEdge ] 


  def extrapolate_from_repeat_unit(self, ReferenceStart, ReferenceEnd, RepeatUnitLength, NewPose, FinalCstName, PdbTag):
    ''' renumbers based on repeat unit pose '''

    # Loop through positions in range of archetype
    # To avoid double counting first only add constraints from archetype residues to 
    # more C-terminal residues 
    NewLength = NewPose.n_residue()
    self.Range = (1, NewLength)
    self.NewPoseStartShift = ReferenceStart - 1 # for 1 indexing

    UnitShiftMultiples = (NewLength / RepeatUnitLength)
    UnitShiftList = [ RepeatUnitLength * Multiple for Multiple in range( UnitShiftMultiples ) ] 
    
    Edge1Cst, Edge2Cst, BothEdgeCst, MiddleCst = self.shift_and_sort_constraints(ReferenceStart, ReferenceEnd, RepeatUnitLength)
    
    # self.output_cst(Edge1Cst, 'Edge1.cst')
    # self.output_cst(Edge2Cst, 'Edge2.cst')
    # self.output_cst(BothEdgeCst, 'BothEdgeCst.cst')
    # self.output_cst(MiddleCst, 'Middle.cst')
    # print 'Edge1Cst:', Edge1Cst, '\n'
    # print 'Edge2Cst:', Edge2Cst, '\n'
    # print 'BothEdgeCst:', BothEdgeCst, '\n'
    # print 'MiddleCst:', MiddleCst, '\n'
    print 'UnitShiftList:', UnitShiftList
    print 'RepeatUnitLength:', RepeatUnitLength
    
    MiddleRepeatCstList = []
    MiddleSkippedCst = 0
    for Constraint in MiddleCst:
      AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
      # Loops through all repeat positions corresponding to reference position
      for Shift in UnitShiftList:
        # print 'Shift:', Shift
        # print 'AtomResidueCoords:', AtomResidueCoords
        ShiftedAtomResidueCoords = [ (AtomName, ResidueNumber+Shift) for AtomName, ResidueNumber in AtomResidueCoords ]
        if pose_has(NewPose, ShiftedAtomResidueCoords):
          MiddleRepeatCstList.append( ( ShiftedAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType ) )
        else:
          MiddleSkippedCst += 1
          # print 'Skipping constraint involving:', ShiftedAtomResidueCoords

    Edge1RepeatCstList = []
    Edge1SkippedCst = 0
    for Constraint in Edge1Cst:
      AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
      for Shift in UnitShiftList[1:]:
        ShiftedAtomResidueCoords = [ (AtomName, ResidueNumber+Shift) for AtomName, ResidueNumber in AtomResidueCoords ]
        if pose_has(NewPose, ShiftedAtomResidueCoords):        
          Edge1RepeatCstList.append( ( ShiftedAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType ) )
        else:
          Edge1SkippedCst += 1
          # print 'Skipping constraint involving:', ShiftedAtomResidueCoords
    
    Edge2RepeatCstList = []     
    Edge2SkippedCst = 0
    for Constraint in Edge2Cst:
      AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
      for Shift in UnitShiftList[:-1]:
        ShiftedAtomResidueCoords = [ (AtomName, ResidueNumber+Shift) for AtomName, ResidueNumber in AtomResidueCoords ]
        if pose_has(NewPose, ShiftedAtomResidueCoords):       
          Edge2RepeatCstList.append( ( ShiftedAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType ) )
        else:
          Edge2SkippedCst += 1
          # print 'Skipping constraint involving:', ShiftedAtomResidueCoords

    BothEdgeRepeatCstList = []
    BothEdgeSkippedCst = 0
    for Constraint in BothEdgeCst:
      AtomResidueCoords, ConstraintParameters, CstLineNumber, CstType = Constraint
      
      for Shift in UnitShiftList:
        ShiftedAtomResidueCoords = [ (AtomName, ResidueNumber+Shift) for AtomName, ResidueNumber in AtomResidueCoords ]
        if pose_has(NewPose, ShiftedAtomResidueCoords):
          BothEdgeRepeatCstList.append( ( ShiftedAtomResidueCoords, ConstraintParameters, CstLineNumber, CstType ) )
        else:
          BothEdgeSkippedCst += 1
          # print 'Skipping constraint involving:', ShiftedAtomResidueCoords

    # RepPose.constraint_set().show_definition(ostream(sys.stdout), RepPose )

    self.output_cst(MiddleRepeatCstList, '%s_MidRepTemp.cst'%PdbTag)
    self.output_cst(Edge1RepeatCstList, '%s_Edge1RepTemp.cst'%PdbTag)
    self.output_cst(Edge2RepeatCstList, '%s_Edge2RepTemp.cst'%PdbTag)
    self.output_cst(BothEdgeRepeatCstList, '%s_BothEdgeRepTemp.cst'%PdbTag)
    
    AllRepeatCst = Edge1RepeatCstList[:]
    AllRepeatCst.extend(Edge1RepeatCstList)
    AllRepeatCst.extend(Edge2RepeatCstList)
    AllRepeatCst.extend(BothEdgeRepeatCstList)
    self.output_cst(AllRepeatCst, 'AllRepeatCst.cst')

    ''' trying out constraints to pick between edge 1 and edge 2 (and filter?) '''
    # print
    # print 'MiddleSkippedCst', MiddleSkippedCst
    # print 'Edge1SkippedCst', Edge1SkippedCst
    # print 'Edge2SkippedCst', Edge2SkippedCst
    # print 'BothEdgeSkippedCst', BothEdgeSkippedCst
    # print
    # print 'MiddleRepeatCst' 
    NumberMiddleRepeatCst = len(MiddleRepeatCstList)
    # print 'Edge1RepeatCst'
    NumberEdge1RepeatCst = len(Edge1RepeatCstList)
    # print 'Edge2RepeatCst'
    NumberEdge2RepeatCst = len(Edge2RepeatCstList)
    # print 'BothEdgeRepeatCst'
    NumberBothEdgeRepeatCst = len(BothEdgeRepeatCstList)

    NumberAllRepeatCst = len(AllRepeatCst)
    
    # # All default talaris 2013 non zero weights set to zero
    CstScoreFunction = set_all_weights_zero( rosetta.getScoreFunction() )
    # # turning on constraint weights
    CstScoreFunction.set_weight( rosetta.atom_pair_constraint, 1.0 )
    CstScoreFunction.set_weight( rosetta.angle_constraint, 1.0 )
    CstScoreFunction.set_weight( rosetta.dihedral_constraint, 1.0 )

    print 'MiddlePose should have %d constraints !!! '%NumberMiddleRepeatCst 
    MiddlePose = NewPose.clone()
    if NumberEdge1RepeatCst:
      ConstraintSetter = rosetta.ConstraintSetMover()
      ConstraintSetter.constraint_file('%s_MidRepTemp.cst'%PdbTag) 
      ConstraintSetter.apply(MiddlePose) 
      # return ConstraintSetter
      # return MiddlePose
      CstScoreFunction.show(MiddlePose)
      # MiddlePose.constraint_set().show_definition(ostream(sys.stdout), MiddlePose )
      print

    print 'Edge1Pose should have %d constraints !!! '%NumberEdge1RepeatCst   
    Edge1Pose = NewPose.clone()
    if NumberEdge1RepeatCst:
      ConstraintSetter = rosetta.ConstraintSetMover()
      ConstraintSetter.constraint_file('%s_Edge1RepTemp.cst'%PdbTag) 
      ConstraintSetter.apply(Edge1Pose) 
      CstScoreFunction.show(Edge1Pose)
      Edge1Score = CstScoreFunction(Edge1Pose)
      Edge1ScoreNorm = Edge1Score / NumberEdge1RepeatCst
      # Edge1Pose.constraint_set().show_definition(ostream(sys.stdout), Edge1Pose )
      print

    print 'Edge2Pose should have %d constraints !!! '%NumberEdge2RepeatCst   
    Edge2Pose = NewPose.clone()
    if NumberEdge2RepeatCst:
      ConstraintSetter = rosetta.ConstraintSetMover()
      ConstraintSetter.constraint_file('%s_Edge2RepTemp.cst'%PdbTag) 
      ConstraintSetter.apply(Edge2Pose) 
      CstScoreFunction.show(Edge2Pose)
      Edge2Score = CstScoreFunction(Edge2Pose)
      Edge2ScoreNorm = Edge2Score / NumberEdge2RepeatCst
      # Edge2Pose.constraint_set().show_definition(ostream(sys.stdout), Edge2Pose )
      print 
    
    # print 'BothEdgePose should have %d constraints !!! '%NumberBothEdgeRepeatCst   
    # BothEdgePose = NewPose.clone()
    # if NumberBothEdgeRepeatCst:
    #   ConstraintSetter = rosetta.ConstraintSetMover()
    #   ConstraintSetter.constraint_file('%s_AllRepeatCstTemp.cst'%PdbTag) 
    #   ConstraintSetter.apply(BothEdgePose) 
    #   CstScoreFunction.show(BothEdgePose)
    #   BothEdgeScore = CstScoreFunction(BothEdgePose)
    #   BothEdgeScoreNorm = BothEdgeScore / NumberBothEdgeRepeatCst
    #   # BothEdgePose.constraint_set().show_definition(ostream(sys.stdout), BothEdgePose )
    #   print 

    # print 'AllCstPose should have %d constraints !!! '%NumberAllRepeatCst   
    # AllCstPose = NewPose.clone()
    # ConstraintSetter = rosetta.ConstraintSetMover()
    # ConstraintSetter.constraint_file('%s_AllRepeatCstTemp.cst'%PdbTag) 
    # ConstraintSetter.apply(AllCstPose) 
    # CstScoreFunction.show(AllCstPose)
    # # AllCstPose.constraint_set().show_definition(ostream(sys.stdout), AllCstPose )
    # print 

    CuratedRepeatCst = MiddleRepeatCstList[:]
    ## whether these should be included or not is up in the air!!
    CuratedRepeatCst.extend(BothEdgeRepeatCstList)

    if NumberEdge1RepeatCst and NumberEdge2RepeatCst:
      if Edge1ScoreNorm <= Edge2ScoreNorm:
        CuratedRepeatCst.extend(Edge1RepeatCstList)
      else:
        CuratedRepeatCst.extend(Edge2RepeatCstList)

    elif NumberEdge1RepeatCst:
      CuratedRepeatCst.extend(Edge1RepeatCstList)
    elif NumberEdge2RepeatCst:
      CuratedRepeatCst.extend(Edge2RepeatCstList)

    # CuratedRepeatCst
    # print 'Edge1ScoreNorm, Edge2ScoreNorm', Edge1ScoreNorm, Edge2ScoreNorm
    # self.output_cst(CuratedRepeatCst, FinalCstName)

    AllWithEdge1RepeatCst = MiddleRepeatCstList[:]
    ## whether these should be included or not is up in the air!!
    # AllWithEdge1RepeatCst.extend(BothEdgeRepeatCstList)
    AllWithEdge1RepeatCst.extend(Edge1RepeatCstList)

    AllWithEdge2RepeatCst = MiddleRepeatCstList[:]
    ## whether these should be included or not is up in the air!!
    # AllWithEdge2RepeatCst.extend(BothEdgeRepeatCstList)
    AllWithEdge2RepeatCst.extend(Edge2RepeatCstList)

    ModFinalCstName = (FinalCstName+'!').replace('.cst!', '')
    self.output_cst(AllWithEdge1RepeatCst, ModFinalCstName+'_e1.cst')
    self.output_cst(AllWithEdge2RepeatCst, ModFinalCstName+'_e2.cst')

    RemainingTempFiles = glob.glob( '%s_*Temp.cst'%PdbTag )
    for File in RemainingTempFiles:
      subprocess.check_output(['rm', File])

def pose_has(Pose, AtomResidueCoords):
  for AtomName, ResidueNumber in AtomResidueCoords:
    if ResidueNumber < 1:
      # print AtomResidueCoords
      # print 'check now!!!'
      # time.sleep(2)
      return False
    try:
      Residue = Pose.residue(ResidueNumber)
    except RuntimeError:
      return False
    # Continue only if
    # residue has atom
    if Residue.has(AtomName):
      # and this atom is not virtual
      if Residue.is_virtual(Residue.atom_index(AtomName)):
        return False
      else:
        continue

    else:
      return False
  return True


def set_all_weights_zero(ScoreFunction):
  ScoreFunction.set_weight(rosetta.fa_atr, 0.0) 
  ScoreFunction.set_weight(rosetta.fa_rep, 0.0) 
  ScoreFunction.set_weight(rosetta.fa_sol, 0.0) 
  ScoreFunction.set_weight(rosetta.fa_intra_rep, 0.0) 
  ScoreFunction.set_weight(rosetta.fa_elec, 0.0) 
  ScoreFunction.set_weight(rosetta.pro_close, 0.0) 
  ScoreFunction.set_weight(rosetta.hbond_sr_bb, 0.0) 
  ScoreFunction.set_weight(rosetta.hbond_lr_bb, 0.0) 
  ScoreFunction.set_weight(rosetta.hbond_bb_sc, 0.0) 
  ScoreFunction.set_weight(rosetta.hbond_sc, 0.0) 
  ScoreFunction.set_weight(rosetta.dslf_fa13, 0.0) 
  ScoreFunction.set_weight(rosetta.rama, 0.0) 
  ScoreFunction.set_weight(rosetta.omega, 0.0) 
  ScoreFunction.set_weight(rosetta.fa_dun, 0.0) 
  ScoreFunction.set_weight(rosetta.p_aa_pp, 0.0) 
  ScoreFunction.set_weight(rosetta.ref, 0.0) 
  return ScoreFunction

# # embrassingly parallel function for  
# def extrapolate_constraints(Pose, CstFile):
#   ConstrainerInstance = constraint_extrapolator(CstFile)


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  # Arg block
  ArgParser = argparse.ArgumentParser(description=' expand_constraints.py ( -help ) %s'%InfoString)
  # Required args
  ArgParser.add_argument('-ref_pdb', type=str, help=' reference pdb ', required=True)
  ArgParser.add_argument('-ref_cst', type=str, help=' corresponding to reference pdb ', required=True)
  ArgParser.add_argument('-repeat_pdb_tag', type=str, help=' input pdb tag ', required=True)
  # Optional args
  ArgParser.add_argument('-out', type=str, help=' Output directory ', default='./')
  Args = ArgParser.parse_args()
  if Args.out [-1] != '/':
    Args.out = Args.out + '/'


  # default talaris 2013 score function
  ScoreFunction = rosetta.getScoreFunction()
  # turning on constraint weights
  ScoreFunction.set_weight(rosetta.atom_pair_constraint, 1.0)
  ScoreFunction.set_weight(rosetta.angle_constraint, 1.0)
  ScoreFunction.set_weight(rosetta.dihedral_constraint, 1.0)

  RefPdb = Args.ref_pdb
  # print RefPdb
  ReferencePose = rosetta.pose_from_pdb( RefPdb )
  print 'ReferencePose', ReferencePose

  # modify rosetta cst w/o rosetta
  Constrainer = constraint_extrapolator(Args.ref_cst)

  # RefCst = Args.ref_cst
  # # make constraint mover
  # Constrainer = rosetta.ConstraintSetMover()
  # # get constraints from file
  # Constrainer.constraint_file(RefCst)  
  # # Apply constraints to pose
  # Constrainer.apply(ReferencePose)

  # return Constrainer

  Pdbs = glob.glob( '*%s*.pdb'%Args.repeat_pdb_tag ) 
  assert len(Pdbs), r"No pdbs found with glob: \n %s \n '* % s *.pdb' % Args.repeat_pdb_tag "%Args.repeat_pdb_tag
  
  for Pdb in Pdbs:
    ## For debug put pdb of interest here:
    # if Pdb == 'src15_38__22_45_rep24_1EZG.pdb':

    print 'Pdb:', Pdb 
    Pose = rosetta.pose_from_pdb(Pdb)

    try: 
      SourceRangeString = re.sub(r'.*src(\d+_\d+__\d+_\d+).*pdb', r'\1', Pdb)
      SourceRanges = [ [ int(Number) for Number in Range.split('_') ]   for Range in SourceRangeString.split('__') ]
    except ValueError:
      print 'No src range tag, skipping: %s '%Pdb
      continue

    print 'SourceRanges:', SourceRanges
    RepeatLength = int( re.sub(r'.*rep(\d+).*pdb', r'\1', Pdb) )
    print 'RepeatLength', RepeatLength
    print
    
    # print [Pdb]
    PdbTag = (Pdb+'!').replace('.pdb!', '').replace('!', '')
    CstName = PdbTag+'.cst'
    ExtrapolatedConstraints = Constrainer.extrapolate_from_repeat_unit(SourceRanges[0][0], SourceRanges[0][1], RepeatLength, Pose, CstName, PdbTag)
    

if __name__ == "__main__":
  sys.exit(main())
