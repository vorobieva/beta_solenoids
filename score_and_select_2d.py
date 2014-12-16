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
  from expand_constraints import set_all_weights_zero

# '''

class plotly_plotter:
  def __init__(self, User, Key, ScoreFxns=[], FxnNames=[], PerResidue=True):
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
    self.ScoreTraces = []

    self.ColorIterator = 0
    # self.Colors = [] ##### add crayons here

    self.MaxX = 0

    # Later keyed with index of self.ScoreFxns
    self.ScoreFunctionScoredPdbs = {}

  def clear_traces(self):
    self.ScoreTraces = []

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

    while len( glob.glob('%s.cst'%StemName) ) == 0:
      StemName = StemName[:-1]
    assert len( glob.glob('%s.cst'%StemName) ) == 1, 'ambigous cst'

    CstName = glob.glob('%s.cst'%StemName)[0]
    return CstName

  def score_poses(self, Poses, Cst=0, Name='', Color=''):
    ''' Give a list of rosetta poses '''

    PdbNames = [ Pose.pdb_info().name() for Pose in Poses ] 
    
    if Cst == 1:
      CstNames = [ self.glob_cst_file(Pdb) for Pdb in PdbNames]
      print 
      print 'PdbNames', PdbNames
      print 'CstNames', CstNames
      for i, Pose in enumerate(Poses):
        # make constraint mover
        Constrainer = rosetta.ConstraintSetMover()
        # get constraints from file
        Constrainer.constraint_file(CstNames[i])
        Constrainer.apply(Pose)

    # self.FxnNames
    ScoreFunctionScores = []
    
    for i, Fxn in enumerate(self.ScoreFxns):
      if self.PerRes:
        PoseScores = [ Fxn(Pose) / Pose.n_residue() for Pose in Poses ]
      else:
        PoseScores = [ Fxn(Pose) for Pose in Poses ]

      if len(Name):
        PlotName = Name+'  '+self.FxnNames[i]
      else:
        PlotName = self.FxnNames[i]

      ScoreFunctionScores.append(PoseScores)
      ScoredPdbList = [ (PoseScores[j], Pose.pdb_info().name() ) for j, Pose in enumerate(Poses)]

      try:
        self.ScoreFunctionScoredPdbs[i].extend( ScoredPdbList )
      except:
        self.ScoreFunctionScoredPdbs[i] = ScoredPdbList 

    FxnCombos = itertools.combinations( range(0, len(ScoreFunctionScores) ), 2 )
    AllTraces = []
    TraceInfo = []

    for Combo in FxnCombos:
      x = Combo[0]
      y = Combo[1] 

      ComboTrace = self.Graph.Scatter(
          x=ScoreFunctionScores[x],
          y=ScoreFunctionScores[y],
          text=PdbNames,
          mode='markers',
          name=Name
      )

      self.ScoreTraces.append(ComboTrace)

      MaxX = max( ScoreFunctionScores[x] )

      if MaxX > self.MaxX:
        self.MaxX = MaxX

  def add_comparsion_threshold(self, Pose, Function, Xaxis):
    ''' for adding native pose energy line to compare designs against'''

    # for i, Fxn in enumerate(self.ScoreFxns):
    if self.PerRes:
      NativeScore = Function(Pose) / Pose.n_residue()
    else:
      NativeScore = Function(Pose)

    NativePoseScores = [ NativeScore for j in range(len(Xaxis)) ]

    CompareTrace = self.Graph.Scatter(
    x = Xaxis,
    y = NativePoseScores,
    name = Pose.pdb_info().name()+' native',
    mode ='lines' )
    self.ScoreTraces.append(CompareTrace)


  def plot_traces(self, PlotName='Unamed', Xaxis='', Yaxis=''):
    ''' Plot premade plotly traces. Used within plot_scores, but 
    can also be used on larger collection of traces '''
    # Import plotly for ploting 
    import plotly.plotly as py
    
    # Sign in with class login info. Change in def __init__ above
    py.sign_in( self.User, self.ApiKey )

    data = self.Graph.Data(self.ScoreTraces)
    layout = self.Graph.Layout(

        xaxis=self.Graph.XAxis(
            title=Xaxis, 
            showgrid=True,
            zeroline=False
        ),
        yaxis=self.Graph.YAxis(
            title=Yaxis, 
            showgrid=True,
            zeroline=False
        )
    )
    fig = self.Graph.Figure(data=data, layout=layout)
    plot_url = self.py.plot(fig, filename=PlotName)


# sys.argv = [ sys.argv[0], '-pdb_glob', 'src*_Relax*Relax.pdb', '-native_pdb', '1M8N_Relax.pdb', '-out', 'LowNrgRepeats', '-score', '-1.5' ]


def main(argv=None):
  if argv is None:
    argv = sys.argv
  ArgParser = argparse.ArgumentParser(description=" for plotting pdb scores and selecting subsets based on absolute or per residue scores ")
  ArgParser.add_argument('-pdb_glob', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  ArgParser.add_argument('-native_pdb', type=str, help=" pdb to compare designs against ", required=True )    
  ArgParser.add_argument('-out', type=str, help=" folder to move files to ", required=True )    
  ArgParser.add_argument('-score', type=float, help=" select all structures with less than this REU / residue ", default=None )
  ArgParser.add_argument('-norm', type=int, help=" 0|(1) normalize scores by residue ", default=1 )
  # following args are for plotly:
  ArgParser.add_argument('-plot', type=int, help=" 0|(1) plot scores with plotly ", default=1 )
  ArgParser.add_argument('-plotly_id', type=str, help=" pdb stem, start of globs for pdbs and csts ", default="pylesharley") # required=True )    
  ArgParser.add_argument('-plotly_key', type=str, help=" pdb stem, start of globs for pdbs and csts ", default="cc5z4a8kst") # required=True )    
  ArgParser.add_argument('-name', type=str, help=" plot tag ", default='' )
  ArgParser.add_argument('-multi', type=int, help=" 0|(1) plot different methods together on same plot ", default=1 )
  Args = ArgParser.parse_args()
  Pdbs = glob.glob( Args.pdb_glob )
  print 'globed %d pdbs'%len(Pdbs)

  if not os.path.isdir(Args.out):
    subprocess.check_output(['mkdir', Args.out])
  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  if Args.name != '':
    Args.out = Args.out + ' '

  NativePose = rosetta.pose_from_pdb( Args.native_pdb )

  RepeatLengths = []
  ProcessTags = {}
  TagList = []
  TagByPdbName = {}

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
    RepeatLength = int(re.sub(r'^.*rep(\d+).*pdb$', r'\1', Pdb))
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

  # set up scorefunctions for plotting
  Talaris = rosetta.getScoreFunction()
  CstScore = set_all_weights_zero( rosetta.getScoreFunction() )
  CstScore.set_weight(rosetta.atom_pair_constraint, 10.0)
  CstScore.set_weight(rosetta.angle_constraint, 5.0)
  CstScore.set_weight(rosetta.dihedral_constraint, 3.0)

  # Add more score functions as wanted
  if Args.plot:
    Plotter = plotly_plotter( Args.plotly_id, Args.plotly_key, ScoreFxns=[ CstScore, Talaris ], FxnNames=[ 'ConstraintScore', 'Talaris2013' ], PerResidue=True )

  XaxisSortingTuples = []

  for PoseGroup in AllGroups:
  # for PoseGroup in [SortedTuples]:
    if len(PoseGroup):
      print 
      print 'Group:', PoseGroup
      Poses = [ PoseTuple[-1] for PoseTuple in PoseGroup ]
      print PoseGroup
      RepeatLength = PoseGroup[0][0]
      print '\n'.join( [ Pose.pdb_info().name() for Pose in Poses ] ) 
      # print 'Zero index pose tuple:'
      # print PoseGroup[0]
     
      if Args.plot:
        GroupPdbName = PoseGroup[0][-1].pdb_info().name()
        if Args.multi:
          Tag = TagByPdbName[GroupPdbName] 
          Plotter.score_poses( Poses, 1, Tag )
    
  Plotter.add_comparsion_threshold( NativePose, Talaris, [0, Plotter.MaxX] )
  Plotter.plot_traces( PlotName='%s%s based %d res '%( Args.name, Args.native_pdb, RepeatLength ) )

  # print 'Plotter.ScoreFunctionScoredPdbs'
  # print Plotter.ScoreFunctionScoredPdbs

  ScoreFunctionScoreCutoffs = []

  for i, Name in enumerate( Plotter.FxnNames ):
    while 1:
      try:
        Cutoff = float( raw_input('Enter cutoff value (maximum) for %s function: '%Name) ) 
        break
      except ValueError:
        pass  
    ScoreFunctionScoreCutoffs.append(Cutoff)

  print 'Cutoff values set at:'
  for i, Name in enumerate( Plotter.FxnNames ):
    print Name, ScoreFunctionScoreCutoffs[i]
    Plotter.ScoreFunctionScoredPdbs[i].sort()

  PassingPdbs = []
  for i, Name in enumerate( Plotter.FxnNames ):
    PassThisFxn = []
    Cutoff = ScoreFunctionScoreCutoffs[i]
    print Plotter.ScoreFunctionScoredPdbs[i]
    for Score, Pdb in Plotter.ScoreFunctionScoredPdbs[i]:
      if Score <= Cutoff:
        PassThisFxn.append(Pdb)
      else:
        break
    PassingPdbs.append( PassThisFxn )

  PdbsPassingAll = PassingPdbs[0]
  for OtherSet in PassingPdbs[1:]:
    PdbsPassingAll = list( set(PdbsPassingAll) & set(OtherSet) )
  
  for Pdb in PdbsPassingAll:
    subprocess.check_output([ 'cp', Pdb, Args.out ])

  #   if Args.score != None:
  #     with open('%sScores.log'%Args.out, 'a') as Log:  
  #       for PoseTuple in PoseGroup:
          
  #         print 'PoseTuple', PoseTuple
  #         #### RepLen, Score, Pose
  #         # if Score > Args.score:
  #         #   break
  #         # PdbName = Pose.pdb_info().name()
  #         # subprocess.check_output([ 'cp', PdbName, Args.out ])
  #         # print>>Log, '%s\t%.3f'%(PdbName, Score)

  # if Args.multi:
  #   XaxisSortingTuples.sort()
  #   print 'here', XaxisSortingTuples[-1]
  #   Xaxis = XaxisSortingTuples[-1][1]
  #   Plotter.add_comparsion_threshold( NativePose, Xaxis )
  #   Plotter.plot_traces( PlotName='All %s%s based globed with %s'%(Args.name, Args.native_pdb, Args.pdb_glob) )


###### sys.exit() CstGlober = (Pdb+'!').replace('.pdb!', '*cst')
###### CapCstName = re.sub(r'(.*).pdb$', r'\1.cst', CappedNamePdb)


if __name__ == "__main__":
  sys.exit(main())

