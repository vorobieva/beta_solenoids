#! /usr/bin/env python

# '''
import numpy as np
import subprocess
import itertools
import argparse
import glob
import time
import copy
import sys
import os
import re

if '-h' not in sys.argv:
  import rosetta
  # rosetta.init()
  rosetta.init(extra_options = " -ex1 -ex2 -no_optH false -use_input_sc -mute basic -mute core -mute protocols")
  from rosetta.protocols import grafting 
  # from repo 
  import solenoid_tools
  from expand_cst import set_all_weights_zero

# '''

class plotly_plotter:
  def __init__(self, User, Key, RefPdb, ScoreFxns=[], FxnNames=[], PerResidue=True):
    ''' Track scores of design trajectories for plotly plots '''

    self.User = User
    self.ApiKey = Key

    import plotly.graph_objs as Graph
    import plotly.plotly as py

    self.Graph = Graph
    self.py = py 

    self.ScoreFxns = ScoreFxns
    self.FxnNames = FxnNames
    assert len(self.ScoreFxns) == len(self.FxnNames)

    self.PerRes = PerResidue

    self.RefPdb = RefPdb
    
    self.RefPose = rosetta.pose_from_pdb( RefPdb )

    self.Score2dComboTraces = {}
    # self.ColorIterator = 0
    # self.Colors = [] 

    self.MaxScores = [ 0 for Fxn in self.ScoreFxns ]
    self.MinScores = [ 999 for Fxn in self.ScoreFxns ]

    # Scores ordered in all score lists in same order for ploting
    self.TaggedPoseScores = {}
    self.PoseTags = []
    # Later keyed with index of self.ScoreFxns
    self.ScoreFunctionScoredPdbs = {}

    self.CstDict = {}

  def clear_traces(self):
    self.Score2dComboTraces = {}

  def add_fxn(self, OtherScoreFxn, FxnNames):
    ''' add addional score functions '''
    if type(OtherScoreFxn) == list:
      self.ScoreFxns.extend(OtherScoreFxn)
      self.FxnNames.extend(FxnNames)
    else:
      self.ScoreFxns.append(OtherScoreFxn)
      self.FxnNames.append(FxnNames)

  def glob_cst_file(self, PdbName):
    StemName = re.sub( r'(.*)\.pdb$', r'\1', PdbName )
    assert StemName != PdbName, ' re.sub for pdb stem failed '

    # print 'entering while loop 1', StemName
    while len( glob.glob('%s.cst'%StemName) ) == 0:
      StemName = StemName[:-1]
      try:
        assert len(StemName), 'No cst file found for %s'%PdbName
      except AssertionError:
       return None

    assert len( glob.glob('%s.cst'%StemName) ) == 1, 'ambigous cst'
    # print 'leaving while loop 1'

    CstName = glob.glob('%s.cst'%StemName)[0]
    return CstName

  def score_poses(self, Poses, Cst=0, Tag='', Color=''):
    ''' Give a list of rosetta poses '''

    self.PdbNames = [ Pose.pdb_info().name() for Pose in Poses ] 

    if type(Cst) == str or Cst == 1:
      if type(Cst) == str:
        CstNames = [ Cst for Pdb in self.PdbNames ]
      elif Cst == 1:
        CstNames = [ self.glob_cst_file(Pdb) for Pdb in self.PdbNames]

      for i, Cst in enumerate(CstNames):
        self.CstDict[self.PdbNames[i]] = Cst

      for i, Pose in enumerate(Poses):
        # make constraint mover
        Constrainer = rosetta.ConstraintSetMover()
        # get constraints from file
        # print CstNames[i]
        if CstNames[i] == None:
          return False
        Constrainer.constraint_file(CstNames[i])
        Constrainer.apply(Pose)

    
    # Loop through all score functions and score poses    
    for i, Fxn in enumerate(self.ScoreFxns):
      if self.PerRes:
        PoseScores = [ Fxn(Pose) / Pose.n_residue() for Pose in Poses ]
      else:
        PoseScores = [ Fxn(Pose) for Pose in Poses ]

      try:
        DictCheck = self.TaggedPoseScores[Tag]
      except KeyError:
        self.TaggedPoseScores[Tag] = {}
        self.PoseTags.append(Tag)

      try:
        self.TaggedPoseScores[Tag][i].extend(PoseScores)
      except KeyError:
        self.TaggedPoseScores[Tag][i] = PoseScores 

      ScoredPdbList = [ (PoseScores[j], Pose.pdb_info().name() ) for j, Pose in enumerate(Poses)]
      try:
        self.ScoreFunctionScoredPdbs[i].extend( ScoredPdbList )
      except KeyError:
        self.ScoreFunctionScoredPdbs[i] = ScoredPdbList 


  def plot_2d_score_combinations(self):

    FxnCombos = itertools.combinations( range(0, len(self.ScoreFxns) ), 2 )
    AllTraces = []
    TraceInfo = []

    print 'self.Score2dComboTraces', 2, self.Score2dComboTraces

    for Combo in FxnCombos:
  
      #            ///
      ComboKey = '%d_%d' % (Combo[0], Combo[1])
      #            | |
      #             M
      #           kevin is 2-D
      #       but he can show you many
      #          two by two by two
  
      print 'ComboKey:', ComboKey
  
      for Tag in self.PoseTags:

        ComboTrace = self.Graph.Scatter(
          x=self.TaggedPoseScores[Tag][ Combo[0] ],
          y=self.TaggedPoseScores[Tag][ Combo[1] ],
          text=self.PdbNames,
          mode='markers',
          name=Tag                    )
  
        try:
          self.Score2dComboTraces[ComboKey].append( ComboTrace )
        except:
          self.Score2dComboTraces[ComboKey] = [ ComboTrace ]

        if max( self.TaggedPoseScores[Tag][ Combo[0] ] ) > self.MaxScores[ Combo[0] ]:
          self.MaxScores[ Combo[0] ] = max( self.TaggedPoseScores[Tag][ Combo[0] ] )
        if max( self.TaggedPoseScores[Tag][ Combo[1] ] ) > self.MaxScores[ Combo[1] ]:
          self.MaxScores[ Combo[1] ] = max( self.TaggedPoseScores[Tag][ Combo[1] ] )
        
        if min( self.TaggedPoseScores[Tag][ Combo[0] ] ) < self.MinScores[ Combo[0] ]:
          self.MinScores[ Combo[0] ] = min( self.TaggedPoseScores[Tag][ Combo[0] ] )
        if min( self.TaggedPoseScores[Tag][ Combo[1] ] ) < self.MinScores[ Combo[1] ]:
          self.MinScores[ Combo[1] ] = min( self.TaggedPoseScores[Tag][ Combo[1] ] )

    print 'self.Score2dComboTraces', 2.5, self.Score2dComboTraces


  def draw_comparisons(self):
    
    self.RefScores = [ Fxn(self.RefPose) for Fxn in self.ScoreFxns ]
    if self.PerRes:
      self.RefScores = [ ( Score / self.RefPose.n_residue() ) for Score in self.RefScores ]

    print 'self.Score2dComboTraces', 3, self.Score2dComboTraces

    for ComboKey in self.Score2dComboTraces:
      Combo = [ int(i) for i in ComboKey.split('_') ]
      Xaxis = [ self.MinScores[Combo[0]], self.MaxScores[Combo[0]] ]
      NativePoseScores = [ self.RefScores[Combo[1]] for j in range(len(Xaxis)) ]

      CompareTrace = self.Graph.Scatter(
      x = Xaxis,
      y = NativePoseScores,
      name = self.RefPdb,
      mode ='lines' )
      self.Score2dComboTraces[ComboKey].append(CompareTrace)

      Yaxis = [ self.MinScores[Combo[1]], self.MaxScores[Combo[1]] ]
      NativePoseScores = [ self.RefScores[Combo[0]] for j in range(len(Yaxis)) ]

      CompareTrace = self.Graph.Scatter(
      x = NativePoseScores,
      y = Yaxis,
      name = self.RefPdb,
      mode ='lines' )
      self.Score2dComboTraces[ComboKey].append(CompareTrace)

    print 'self.Score2dComboTraces', 4, self.Score2dComboTraces


  def render_scatter_plot(self, PlotName=''):
    ''' Plot premade plotly traces. Used within plot_scores, but 
    can also be used on larger collection of traces '''
    # Import plotly for ploting 
    import plotly.plotly as py
    
    # Sign in with class login info. Change in def __init__ above
    py.sign_in( self.User, self.ApiKey )
    print 'out of loop'
    print 'self.Score2dComboTraces', 4, self.Score2dComboTraces
    for ComboKey in self.Score2dComboTraces:
      print 'plotting combination: ', ComboKey
      print 'traces: ', self.Score2dComboTraces[ComboKey]
      data = self.Graph.Data(self.Score2dComboTraces[ComboKey])
      ComboList = [ int(i) for i in ComboKey.split('_') ]
      x = ComboList[0]
      y = ComboList[1]
      if len(PlotName):
        ComboName = '%s vs %s; %s'%(self.FxnNames[y], self.FxnNames[x], PlotName)
      else:
        ComboName = '%s vs %s; %s'%(self.FxnNames[y], self.FxnNames[x], PlotName)

      layout = self.Graph.Layout(

          xaxis=self.Graph.XAxis(
              title=self.FxnNames[x], 
              showgrid=True,
              zeroline=False
          ),
          yaxis=self.Graph.YAxis(
              title=self.FxnNames[y], 
              showgrid=True,
              zeroline=False
          )
      )
      fig = self.Graph.Figure(data=data, layout=layout)
      print 'assembled figure'
      plot_url = self.py.plot(fig, filename=ComboName)
      print 'sent figure to plotly'

# sys.argv = [ sys.argv[0], '-pdb_glob', 'A011__src144_163__151_170_rep40_3petA_2rASP_plusGLU_ladderASN_MetLeuCore_plusArgBkUpII_Relax00*.pdb', '-param', 'CA9.params', 'CO3.params', '-native', 'A011__src144_163__151_170_rep40_3petA_2rASP_plusGLU_ladderASN_MetLeuCore_plusArgBkUpII_SolutionState.pdb','-cst', 'A011__src144_163__151_170_rep20_3petA_2AspRow.cst', '-name', 'interface_test']
# score_and_select_2d.py -pdb_glob 'A011__src144_163__151_170_rep*Relax*pdb' -native A011__src144_163__151_170_rep40_3petA_2rASP_plusGLU_ladderASN_MetLeuCore_plusArgBkUpII_SolutionState.pdb -cst A011__src144_163__151_170_rep20_3petA_2AspRow.cst

def main(ExtraResidues=0, ipython=0):
  ### Required args
  ArgParser = argparse.ArgumentParser(description=" for plotting pdb scores and selecting subsets based on absolute or per residue scores ")
  ArgParser.add_argument('-pdb_glob', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  ArgParser.add_argument('-native', type=str, help=" pdb to compare designs against ", required=True )    
  ### Default args
  ArgParser.add_argument('-cst', type=str, help=" to provide cst manually, will apply to all globed pdbs!!! ", default=False )
  ArgParser.add_argument('-param', type=str, nargs='+', help=" params ", default=[] )
  ArgParser.add_argument('-norm', type=int, help=" 0|(1) normalize scores by residue ", default=1 )

  ### following args are for plotly:
  ### change if you use this script!!!
  ArgParser.add_argument('-plotly_id', type=str, help=" ", default="pylesharley") # required=True )    
  ArgParser.add_argument('-plotly_key', type=str, help="  ", default="cc5z4a8kst") # required=True )    
  ArgParser.add_argument('-plot', type=int, help=" 0|(1) plot scores with plotly ", default=1 )
  ArgParser.add_argument('-name', type=str, help=" plot tag ", default='' )
  ArgParser.add_argument('-and_or', type=str, help=" And/Or logic for score cutoffs. Default = 'and'  ", default='and' )
  ArgParser.add_argument('-multi', type=int, help=" 0|(1) plot different methods together on same plot ", default=1 )
  
  Args = ArgParser.parse_args()
  Pdbs = glob.glob( Args.pdb_glob )
  print 'globed %d pdbs'%len(Pdbs)

  if ExtraResidues == 0 and len(Args.param) > 0:
    try: 
      ExtraParams = rosetta.Vector1( Args.param )
      ExtraResidues = rosetta.generate_nonstandard_residue_set( ExtraParams )
    except:
      ExtraParams = rosetta.Vector1( Args.param )
      ExtraResidues = rosetta.generate_nonstandard_residue_set( ExtraParams )
    ### for ipython mode
    if ipython: 
      return ExtraResidues

  Args.and_or = Args.and_or.lower()
  assert Args.and_or == 'and' or Args.and_or == 'or', " -and_or must equal 'and' or 'or' "

  RepeatLengths = []
  ProcessTags = {}
  TagList = []
  TagByPdbName = {}

  # better to find out of native pdb is wrong before waiting for pdb scoring
  Check = open(Args.native, 'r')

  # print ' first loop '
  OverlapStarts = []
  for Pdb in Pdbs:
    Tag = re.sub(r'^.*rep\d+(.*)\.pdb$', r'\1', Pdb)
    for OtherPdb in Pdbs:
      OtherTag = re.sub(r'^.*rep\d+(.*)\.pdb$', r'\1', Pdb)
      i = 0
      if Pdb != OtherPdb:
        while Pdb[:i] == OtherPdb[:i]:
          i+=1
        Overlap = OtherPdb[:i-1]
        OverlapStarts.append( ( len(Overlap), Overlap ) )

  OverlapStarts.sort()
  ShortestOverlap = OverlapStarts[0][1]

  # print 'OverlapStarts', OverlapStarts
  # print 'ShortestOverlap', ShortestOverlap
  
  for Pdb in Pdbs:
    try:
      RepeatLength = int(re.sub(r'^.*rep(\d+).*pdb$', r'\1', Pdb))
    except ValueError:
      RepeatLength = 0
    # SourceStart = int(re.sub(r'^.*src(\d+).*pdb$', r'\1', Pdb))
    assert RepeatLength != Pdb, " regular expression extraction of 'rep' (repeat length) value failed on %s "%Pdb 
    # assert SourceStart != Pdb and RepeatLength != Pdb, ' regular expression extraction of rep or src value failed on %s '%Pdb 
    RepeatLengths.append(RepeatLength)    

    #### re.sub out tag from design process
    Tag = re.sub(r'^.*rep\d+(.*)\.pdb$', r'\1', Pdb)
    Tag = re.sub(r'^%s(.*)\.pdb$'%(ShortestOverlap), r'\1', Tag)
    
    TagByPdbName[Pdb] = Tag
    try:
      TagNumber = ProcessTags[Tag] 
    except:
      TagNumber = len(ProcessTags) + 1
      ProcessTags[Tag] = TagNumber
    TagList.append(TagNumber)

  # Scoring is redundant, once for sorting outside plotter, then again in plotter
  # making not redundant not a priority. 
  # Scoring in the plotter object is so multiple score functions can be plotted easily

  # Sort by repeat length, then score
  if Args.multi:
    # Sort by repeat length, then method tag, then score
    MultiPoseSortingTuples = []
  else:
    PoseSortingTuples = []

  Talaris = rosetta.getScoreFunction()
  for i, Pdb in enumerate(Pdbs):
    RepeatLength = RepeatLengths[i]
    ProcessNumber = TagList[i]
    Pose = rosetta.pose_from_pdb(Pdb)
    if Args.norm:
      Score = Talaris(Pose) / Pose.n_residue()
    else:
      Score = Talaris(Pose) 
    
    # print 'Pdb', Pdb
    if Args.multi:
      MultiPoseSortingTuples.append( (RepeatLength, ProcessNumber, Score, Pose) )
    else:
      PoseSortingTuples.append( (RepeatLength, Score, Pose) )


  if Args.multi:
    # Sort by repeat length, then method tag, then score
    MultiPoseSortingTuples.sort()
  else:
    # sorts by repeat length (shortest to longest) then score (best to worst)
    PoseSortingTuples.sort()

  if Args.multi:
    # print 'MultiPoseSortingTuples', MultiPoseSortingTuples
    SortedTuples = MultiPoseSortingTuples
  else:
    # print 'PoseSortingTuples', PoseSortingTuples
    SortedTuples = PoseSortingTuples

  LastLength = 0
  LastTag = 0
  AllGroups = []
  CurrentGroup = []

  for PoseTuple in SortedTuples:
    Length = PoseTuple[0]
    if Args.multi:
      Tag = PoseTuple[1]
    
    if LastLength and Length != LastLength:
      AllGroups.append(CurrentGroup)
      CurrentGroup = []
    
    if Args.multi:
      if LastTag and Tag != LastTag:
        AllGroups.append(CurrentGroup)
        CurrentGroup = [] 
    
    CurrentGroup.append(PoseTuple)
    LastLength = Length
    if Args.multi: 
      LastTag = Tag

  # for last repeat length
  AllGroups.append(CurrentGroup)

  ''' Build score functions here: '''

  Talaris = rosetta.getScoreFunction()

  # This line returns a talaris function with all default weights set to 0
  CstScore = set_all_weights_zero( rosetta.getScoreFunction() )
  CstScore.set_weight(rosetta.atom_pair_constraint, 10.0)
  CstScore.set_weight(rosetta.angle_constraint, 5.0)
  CstScore.set_weight(rosetta.dihedral_constraint, 3.0)

  HbondScore = set_all_weights_zero( rosetta.getScoreFunction() )
  HbondScore.set_weight(rosetta.hbond_sr_bb, 1.170)
  HbondScore.set_weight(rosetta.hbond_lr_bb, 1.170)
  HbondScore.set_weight(rosetta.hbond_bb_sc, 1.170)
  HbondScore.set_weight(rosetta.hbond_sc, 1.100)

  Disulfide = set_all_weights_zero( rosetta.getScoreFunction() )
  Disulfide.set_weight(rosetta.dslf_fa13, 1.0)

  if Args.plot:
    if Args.norm:
      PerRes = True
    else:
      PerRes = False
    ''' Add and remove score functions here '''
    Plotter = plotly_plotter( Args.plotly_id, Args.plotly_key, Args.native,
                              ScoreFxns=[ CstScore, Talaris, HbondScore, Disulfide ],
                              FxnNames=[ 'ConstraintScore', 'Talaris2013', 'H-bond', 'Disulfide' ],
                              PerResidue=PerRes )

  XaxisSortingTuples = []

  for PoseGroup in AllGroups:
  # for PoseGroup in [SortedTuples]:
    if len(PoseGroup):
      # print 
      # print 'Group:', PoseGroup
      Poses = [ PoseTuple[-1] for PoseTuple in PoseGroup ]
      # print PoseGroup
      RepeatLength = PoseGroup[0][0]
      # print '\n'.join( [ Pose.pdb_info().name() for Pose in Poses ] ) 
      # print 'Zero index pose tuple:'
      # print PoseGroup[0]
     
      if Args.plot:
        GroupPdbName = PoseGroup[0][-1].pdb_info().name()
        if Args.multi:
          Tag = TagByPdbName[GroupPdbName] 
          
          if Args.cst:
            Plotter.score_poses( Poses, Args.cst, Tag )
          else:
            Plotter.score_poses( Poses, 1, Tag )
  
  # return Plotter
  Plotter.plot_2d_score_combinations()
  print 'Plotter.Score2dComboTraces', 3, Plotter.Score2dComboTraces

  Plotter.draw_comparisons()

  print 'plotting...'
  if len(Args.name):
    Name = Args.name
  else:
    Name = '%s based %d res '%( Args.native, RepeatLength )
  Plotter.render_scatter_plot( PlotName=Name )
  
  while 1:

    ScoreFunctionScoreCutoffs = []
    for i, Name in enumerate( Plotter.FxnNames ):
      while 1:
        try:
          Cutoff = float( raw_input('\tEnter cutoff value (maximum) for %s function: '%Name) ) 
          break
        except ValueError:
          pass  
      ScoreFunctionScoreCutoffs.append(Cutoff)

    print 'Cutoff values set at:'
    for i, Name in enumerate( Plotter.FxnNames ):
      # print Name, ScoreFunctionScoreCutoffs[i]
      Plotter.ScoreFunctionScoredPdbs[i].sort()

    PassingPdbs = []
    for i, Name in enumerate( Plotter.FxnNames ):
      PassThisFxn = []
      Cutoff = ScoreFunctionScoreCutoffs[i]
      # print Plotter.ScoreFunctionScoredPdbs[i]
      for Score, Pdb in Plotter.ScoreFunctionScoredPdbs[i]:
        if Score <= Cutoff:
          PassThisFxn.append(Pdb)
        else:
          break
      PassingPdbs.append( PassThisFxn )

    PdbsPassingAll = PassingPdbs[0]
    if Args.and_or == 'and':
      for OtherSet in PassingPdbs[1:]:
        PdbsPassingAll = list( set(PdbsPassingAll) & set(OtherSet) )
    else:
      for OtherSet in PassingPdbs[1:]:
        PdbsPassingAll = list( set(PdbsPassingAll + OtherSet) )
    
    Outdir = raw_input( '\tEnter folder to copy pdbs that pass these thresholds (%s logic) to: '%Args.and_or ) 

    if not os.path.isdir(Outdir):
      subprocess.check_output(['mkdir', Outdir])
    if Outdir [-1] != '/':
      Outdir = Outdir + '/'

    for Pdb in PdbsPassingAll:
      subprocess.check_output([ 'cp', Pdb, Outdir ])
      if Plotter.CstDict[Pdb] != None:
        subprocess.check_output([ 'cp', Plotter.CstDict[Pdb], Outdir ])

    Continue = str( raw_input( '\tEnter Y to add another set of selection threshold, or anything else to quit: ') ).upper()
    if Continue == 'Y':
      pass
    else:
      break

if __name__ == "__main__":
  sys.exit(main())

