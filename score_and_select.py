#! /usr/bin/env python

# '''
import numpy as np
import subprocess
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
# '''

class plotly_plotter:
  def __init__(self, ScoreFxns=[], FxnNames=[], EnergyPerResidue=True):
    ''' Track scores of design trajectories for plotly plots '''

    self.User = "pylesharley"
    self.ApiKey = "cc5z4a8kst"

    import plotly.graph_objs as Graph
    self.Graph = Graph

    self.ScoreFxns = ScoreFxns
    self.FxnNames = FxnNames

    assert len(self.ScoreFxns) == len(self.FxnNames)

    self.NrgPerRes = EnergyPerResidue
    self.ScoreTraces = []

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

  def score_poses(self, Poses):
    ''' Give a list of rosetta poses '''
    
    Xaxis = [ Pose.pdb_info().name() for Pose in Poses ]
    # Xaxis = [ i for i in range(len(Poses)) ]
    
    for i, Fxn in enumerate(self.ScoreFxns):
      if self.NrgPerRes:
        PoseScores = [ Fxn(Pose) / Pose.n_residue() for Pose in Poses ]
      else:
        PoseScores = [ Fxn(Pose) for Pose in Poses ]

      # Y coordinates is current iteration's score fxn scores
      FxnTrace = self.Graph.Scatter(
      x = Xaxis,
      y = PoseScores,
      name = self.FxnNames[i],
      mode = 'markers' )

      self.ScoreTraces.append(FxnTrace)

    return Xaxis

  def add_comparsion_threshold(self, Pose, Xaxis):
    ''' for adding native pose energy line to compare designs against'''

    for i, Fxn in enumerate(self.ScoreFxns):
      if self.NrgPerRes:
        NativeScore = Fxn(Pose) / Pose.n_residue()
      else:
        NativeScore = Fxn(Pose)

      NativePoseScores = [ NativeScore for j in range(len(Xaxis)) ]

      CompareTrace = self.Graph.Scatter(
      x = Xaxis,
      y = NativePoseScores,
      name = Pose.pdb_info().name()+'_'+self.FxnNames[i],
      mode ='lines' )
      self.ScoreTraces.append(CompareTrace)


  def plot_traces(self, PlotName='Unamed'):
    ''' Plot premade plotly traces. Used within plot_scores, but 
    can also be used on larger collection of traces '''
    # Import plotly for ploting 
    import plotly.plotly as py
    
    # Sign in with class login info. Change in def __init__ above
    py.sign_in( self.User, self.ApiKey )

    PlotlyData = self.Graph.Data(self.ScoreTraces)
    plot_url = py.plot(PlotlyData, filename=PlotName)


# sys.argv = [ sys.argv[0], '-pdb_glob', 'src*Relax*Relax.pdb', '-native_pdb', '4DT5_Relax.pdb', '-out', 'Test']


def main(argv=None):
  if argv is None:
    argv = sys.argv
  ArgParser = argparse.ArgumentParser(description=" for plotting pdb scores and selecting subsets based on absolute or per residue scores ")
  ArgParser.add_argument('-pdb_glob', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  ArgParser.add_argument('-native_pdb', type=str, help=" pdb to compare designs against ", required=True )    
  ArgParser.add_argument('-out', type=str, help=" folder to move files to ", required=True )    
  ArgParser.add_argument('-score', type=float, help=" select all structures with less than this REU / residue ", default=None )
  ArgParser.add_argument('-plot', type=int, help=" 0|(1) plot scores with plotly ", default=1 )
  ArgParser.add_argument('-norm', type=int, help=" 0|(1) normalize scores by residue ", default=1 )
  ArgParser.add_argument('-name', type=str, help=" plot tag ", default='' )

  Args = ArgParser.parse_args()
  print Args
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
  for Pdb in Pdbs:
    RepeatLength = int(re.sub(r'^.*rep(\d+).*pdb$', r'\1', Pdb))
    # SourceStart = int(re.sub(r'^.*src(\d+).*pdb$', r'\1', Pdb))
    assert RepeatLength != Pdb, " regular expression extraction of 'rep' (repeat length) value failed on %s "%Pdb 
    # assert SourceStart != Pdb and RepeatLength != Pdb, ' regular expression extraction of rep or src value failed on %s '%Pdb 
    RepeatLengths.append(RepeatLength)
    # RepeatLengths.append(SourceStart)


  PoseSortingTuples = []
  # Scoring is redundant, once for sorting outside plotter, then again in plotter
  # making not redundant not a priority. 
  # Scoring in the plotter object is so multiple score functions can be plotted easily
  Talaris = rosetta.getScoreFunction()
  for i, Pdb in enumerate(Pdbs):
    RepeatLength = RepeatLengths[i]
    Pose = rosetta.pose_from_pdb(Pdb)
    if Args.norm:
      Score = Talaris(Pose) / Pose.n_residue()
    else:
      Score = Talaris(Pose) 
    PoseSortingTuples.append( (RepeatLength, Score, Pose) )
  # sorts by repeat length (shortest to longest) then score (best to worst)
  PoseSortingTuples.sort()

  # print 'PoseSortingTuples', PoseSortingTuples

  AllRepeatLengthGroups = []
  RepeatRepeatLengthGroup = []
  LastLength = 0
  for PoseTuple in PoseSortingTuples:
    Length = PoseTuple[0]
    if LastLength and Length != LastLength:
      AllRepeatLengthGroups.append(RepeatRepeatLengthGroup)
      RepeatRepeatLengthGroup = []
    RepeatRepeatLengthGroup.append(PoseTuple)
    LastLength = Length
  # for last repeat length
  AllRepeatLengthGroups.append(RepeatRepeatLengthGroup)

  # print 'AllRepeatLengthGroups', AllRepeatLengthGroups

  # Add more score functions as wanted
  if Args.plot:
    Plotter = plotly_plotter(ScoreFxns=[ Talaris ], FxnNames=[ 'Talaris' ], EnergyPerResidue=True )

  for RepeatLengthGroup in AllRepeatLengthGroups:
    print 'RepeatLengthGroup', RepeatLengthGroup
    Poses = [ PoseTuple[2] for PoseTuple in RepeatLengthGroup ]
    RepeatLength = RepeatLengthGroup[0][0]
    if Args.plot:
      Plotter.clear_traces()
      Xaxis = Plotter.score_poses( Poses )
      Plotter.add_comparsion_threshold( NativePose, Xaxis )
      Plotter.plot_traces( PlotName='%s%s based %d res repeats globed with %s'%(Args.name, Args.native_pdb, RepeatLength, Args.pdb_glob) )

    if Args.score != None:
      with open('%sScores.log'%Args.out, 'a') as Log:  
        for RepLen, Score, Pose in RepeatLengthGroup:
          if Score > Args.score:
            break
          PdbName = Pose.pdb_info().name()
          subprocess.check_output([ 'cp', PdbName, Args.out ])
          print>>Log, '%s\t%.3f'%(PdbName, Score)


if __name__ == "__main__":
  sys.exit(main())

