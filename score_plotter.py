#! /usr/bin/env python
InfoString = ''' 
Contains light wrapper object for plotly to track design trajectories '''

class tracker:
  def __init__(self, InputScoreFxnList=[], FxnNames=[]):
    ''' Track scores of design trajectories for plotly plots '''

    self.User = "pylesharley"
    self.ApiKey = "cc5z4a8kst"

    import plotly.graph_objs as Graph
    # self.Graph = Graph

    self.Iterations = []
    self.ScoreFxns = InputScoreFxnList
    self.FxnNames = FxnNames
    self.Scores = [ [] for Fxn in self.ScoreFxns ]
    self.Notes = [ [] for Fxn in self.ScoreFxns ]

    assert len(self.ScoreFxns) == len(self.FxnNames)


  def add_fxn(self, OtherScoreFxn, FxnNames):
    ''' add addional score functions
     WILL DELETE EXISTING SCORES ! '''
    if type(OtherScoreFxn) == list:
      self.ScoreFxns.extend(OtherScoreFxn)
      self.FxnNames.extend(FxnNames)
    else:
      self.ScoreFxns.append(OtherScoreFxn)
      self.FxnNames.append(FxnNames)
    self.Scores = [ [] for Fxn in self.ScoreFxns ]
    self.Notes = [ [] for Fxn in self.ScoreFxns ]


  def score(self, Pose, UserDefIteration=False):
    ''' Give single pose '''
    if UserDefIteration == False:
      if len(self.Iterations): self.Iterations.append( self.Iterations[-1]+1 )
      else: self.Iterations.append(1)
    else:
      self.Iterations.append(UserDefIteration)

    for i, Fxn in enumerate(self.ScoreFxns):
      Score = Fxn(Score)
      self.Scores[i].append(Score)


  def plot_traces(PlotName, ScoreTraces):
    ''' Plot premade plotly traces. Used within plot_scores, but 
    can also be used on larger collection of traces '''
    # Import plotly for ploting 
    import plotly.plotly as py
    
    # Sign in with class login info. Change in def __init__ above
    py.sign_in( self.User, self.ApiKey )

    PlotlyData = Graph.Data(ScoreTraces)
    plot_url = py.plot(PlotlyData, filename=PlotName)


  def plot_scores(self, PlotName, TraceNamePrefix=False, GrabDontPlot=False):
    ''' Uses plotly api to make scatter line plots '''

    # X coordinates are always the design iterations
    Xaxis = np.array( self.Iterations )

    if TraceNamePrefix:
      assert type(TraceNamePrefix) == str, ' Name prefix must be string! '
      PartialName = TraceNamePrefix+'_%s'
    else:
      PartialName = '%s'
    
    ScoreTraces = []
    for i in range(len( self.ScoreFxns )):
      Name = self.FxnNames[i]
      # Y coordinates are score fxn scores
      FxnTrace = Graph.Scatter(
        x = Xaxis,
        y = self.Scores[i],
        name=PartialName%self.FxnNames[i],
        mode='lines+markers'
        )
      ScoreTraces.append(FxnTrace)

    if GrabDontPlot:
      return ScoreTraces

    self.plot_traces( PlotName, ScoreTraces )

