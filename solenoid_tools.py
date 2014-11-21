#! /usr/bin/env python

import rosetta
rosetta.init()

from scipy import spatial
import numpy as np

class alpha_carbon:
  """ Calpha node for wDAG searches of proteins for repeats """
  def __init__(self, Number, CoordinateArray):
    if type(CoordinateArray) == list:
      CoordinateArray = np.array(CoordinateArray)
    self.CoordinateArray = CoordinateArray
    self.Number = int(Number)
    self.DownstreamNeighbors = {} # keyed with residue number, value is displacement vector to that residue's own CA

class pose_wdag:
  """ Init method takes protein and spawns Calpha instances to populate """
  def __init__( self, Pose, MinRepeats=3, DistanceTarget=4.7, DistanceFlex=0.5, AngleFlex=15.0 ):
    # Pose = Pose
    self.MinRepeats = MinRepeats
    self.DistanceTarget = DistanceTarget
    self.DistanceFlex = DistanceFlex
    self.AngleFlex = np.radians(AngleFlex) # converts to radians

    # make instance of alpha_carbon for each row
    self.CalphaInstances = []
    CalphaCoords = []
    for P in range( 1, Pose.n_residue() + 1):
      # print P
      # print Pose.residue(P)
      # print Pose.residue(P).xyz('CA')
      CoordList = [Value for Value in Pose.residue(P).xyz('CA') ]
      # print CoordList
      CoordList = list( Pose.residue(P).xyz('CA') )
      # print CoordList
      # print '\n*2'
      self.CalphaInstances.append( alpha_carbon(P, CoordList) )
      CalphaCoords.append( CoordList )

    self.CalphaArray = np.array(CalphaCoords)


  def find_repeat_chains(self):
    CheckRedundancyChecker = {}
    Repeats = []
    Count = 0
    for Calpha in self.CalphaInstances:
      for NeighborNumber in Calpha.DownstreamNeighbors:
        '''
        Repeat chains start with each of Calphas downstream neighbors 
        and are recursively extended so long as VecA ~= VecB
        
           VecA      VecB
        CA- - - > CA- - - > CA
        '''
        Count += 1
        RepeatChain = [ Calpha.Number ]
        RecursiveCalpha = Calpha
        DownstreamNumber = NeighborNumber
        
        try:
          CheckRedundancyChecker['%d_%d'%(Calpha.Number, DownstreamNumber)]

        except KeyError:
          ChainExtended = 1

        while ChainExtended:
          RepeatChain.append(DownstreamNumber)
          DownstreamCalpha = self.CalphaDict[DownstreamNumber]          
          VectorA = RecursiveCalpha.DownstreamNeighbors[DownstreamNumber]
          
          ChainExtended = 0 # not extended unless valid extension found
          for NextDownstreamNumber in DownstreamCalpha.DownstreamNeighbors:
            
            # To prevent redundant calculations, upstream downstream pairs already investigated aren't considered         
            try:
              CheckRedundancyChecker['%d_%d'%(DownstreamNumber, NextDownstreamNumber)]
            # only residue pairs not in redundancy checking dictionary are considering
            except KeyError:
              VectorB = DownstreamCalpha.DownstreamNeighbors[NextDownstreamNumber]
            
              if spatial.distance.cosine(VectorA, VectorB) < self.AngleFlex:
                CheckRedundancyChecker['%d_%d'%(DownstreamNumber, NextDownstreamNumber)] = True
                DownstreamNumber = NextDownstreamNumber
                RecursiveCalpha = DownstreamCalpha
                ChainExtended = 1
                break

        if len(RepeatChain) >= self.MinRepeats:
          Repeats.append(RepeatChain)
          # print 'select repeat%d, resi %s'%(Count,'+'.join([str(Number).split('.')[0] for Number in RepeatChain]))

    return Repeats

  def find_downstream_neighbors(self):
    ''' To set up DAG Calphas become nodes'''
    self.CalphaDict = {}
    for i, Instance in enumerate(self.CalphaInstances):
      # To make the graph directed, repeat chains will always be arranged from 
      for Subsequent in self.CalphaInstances[i:]:
        # Displacement from lower numbered residue to higher number
        Displacement = Subsequent.CoordinateArray - Instance.CoordinateArray
        # Magnitute of displacement vector
        Distance = (np.sum([Component**2 for Component in Displacement]))**0.5
        # check if displacement of given Instance (upstream) and Subsequent (downstream) has magnitute within target range
        if np.abs(Distance - self.DistanceTarget) <= self.DistanceFlex:
          # if in range record Subsequent as downstream neighbor 
          Instance.DownstreamNeighbors[Subsequent.Number] = Displacement
      # add residue alpha carbon instance to dictionary
      self.CalphaDict[Instance.Number] = Instance

def better_consolidate_repeats(ListOfRepeatPositions):
  
  TandemIndenticalSpacings = {0:[]}
  
  for RepeatPositions in ListOfRepeatPositions:
    RepeatPositions = RepeatPositions[:]
    RepeatSpacings = [ RepeatPositions[i+1] - RepeatPositions[i] for i in range(len(RepeatPositions)-1) ]

    # appending buffer value for added last repeat chain to TandemIndenticalSpacings
    RepeatSpacings.append(0)
    RepeatPositions.append(0)

    LastSpacing = 0
    RepeatChain = [RepeatSpacings[0]]
    Start = RepeatPositions[0]
    
    for i, Spacing in enumerate(RepeatSpacings):
      if Spacing == LastSpacing:
        RepeatChain.append(Spacing)
      else:
        TandemIndenticalSpacings[Start] = RepeatChain
        Start = RepeatPositions[i]
        RepeatChain = [RepeatSpacings[i]]

      LastSpacing = Spacing
  
  # print 'TandemIndenticalSpacings', TandemIndenticalSpacings
  # MaxNumberStart = 0
  # EqualLengthStarts = []
  SortedTandemRepeatKeys = [ Key for Key in TandemIndenticalSpacings ]
  SortedTandemRepeatKeys.sort()

  RepeatLengthFreqHash = {}

  for Position in TandemIndenticalSpacings:
    if Position:
      try: RepeatLengthFreqHash[ TandemIndenticalSpacings[Position][0] ] += 1
      except: RepeatLengthFreqHash[ TandemIndenticalSpacings[Position][0] ] = 1

  print 'RepeatLengthFreqHash', RepeatLengthFreqHash
  RepeatStretchesByLengthHash = {}

  # iterate through lengths in length frequency hash
  for TargetLength in RepeatLengthFreqHash:
    # print 'TargetLength', TargetLength
    RepeatStretches = []
    Stretch = []
    LastPosition = -1000
    # iterate through repeat starts in 
    for Position in SortedTandemRepeatKeys:  
      if Position:
        PositionRepeatLength = TandemIndenticalSpacings[Position][0]
        if PositionRepeatLength == TargetLength:
          if Position < LastPosition + TargetLength:
            # print 'if'
            Stretch = Stretch[:-1]
            Stretch.extend([LastPosition, Position])
          else:
            if len(Stretch) > 1:
              RepeatStretches.append(Stretch)
              Stretch = []
            # print 'else'
            Stretch = [Position]
          # print 'Stretch', Stretch
          LastPosition = Position
    
    if len(Stretch) > 1:
      RepeatStretches.append(Stretch)
      Stretch = []
    # except NameError:
    #   print ' Strech doesnt exist yet '
    #   pass

    # print 'TargetLength', TargetLength
    # print 'RepeatStretches', RepeatStretches
    RepeatStretchesByLengthHash[TargetLength] = RepeatStretches

  return RepeatStretchesByLengthHash, TandemIndenticalSpacings
  # for Start in TandemIndenticalSpacings:
  #   if len(TandemIndenticalSpacings[Start]) > len(TandemIndenticalSpacings[MaxNumberStart]):
  #     MaxNumberStart = Start
  #     EqualLengthStarts = []
  #   elif len(TandemIndenticalSpacings[Start]) == len(TandemIndenticalSpacings[MaxNumberStart]):
  #     EqualLengthStarts.append(Start)

  # for RepeatStart in EqualLengthStarts:
  #   try:
  #     assert TandemIndenticalSpacings[MaxNumberStart][0] == TandemIndenticalSpacings[RepeatStart][0], ' different repeat spacings have same max copy number ' 
  #   # This is explicted (semi) silenced to prevent a large job from stopping at some later date
  #   except AssertionError:
  #     print '\n ERROR: multiple different repeat spacings have max copy number.\n'

  # MaxNumberRepeatStarts = [MaxNumberStart] + EqualLengthStarts

  # return MaxNumberRepeatStarts, TandemIndenticalSpacings


def consolidate_repeats(ListOfRepeatPositions):
  
  TandemIndenticalSpacings = {0:[]}
  
  for RepeatPositions in ListOfRepeatPositions:
    RepeatPositions = RepeatPositions[:]
    RepeatSpacings = [ RepeatPositions[i+1] - RepeatPositions[i] for i in range(len(RepeatPositions)-1) ]

    # appending buffer value for added last repeat chain to TandemIndenticalSpacings
    RepeatSpacings.append(0)
    RepeatPositions.append(0)

    LastSpacing = 0
    RepeatChain = [RepeatSpacings[0]]
    Start = RepeatPositions[0]
    
    for i, Spacing in enumerate(RepeatSpacings):
      if Spacing == LastSpacing:
        RepeatChain.append(Spacing)
      else:
        TandemIndenticalSpacings[Start] = RepeatChain
        Start = RepeatPositions[i]
        RepeatChain = [RepeatSpacings[i]]

      LastSpacing = Spacing
  
  print 'TandemIndenticalSpacings', TandemIndenticalSpacings
  MaxNumberStart = 0
  EqualLengthStarts = []

  for Start in TandemIndenticalSpacings:
    if len(TandemIndenticalSpacings[Start]) > len(TandemIndenticalSpacings[MaxNumberStart]):
      MaxNumberStart = Start
      EqualLengthStarts = []
    elif len(TandemIndenticalSpacings[Start]) == len(TandemIndenticalSpacings[MaxNumberStart]):
      EqualLengthStarts.append(Start)

  for RepeatStart in EqualLengthStarts:
    try:
      assert TandemIndenticalSpacings[MaxNumberStart][0] == TandemIndenticalSpacings[RepeatStart][0], ' different repeat spacings have same max copy number ' 
    # This is explicted (semi) silenced to prevent a large job from stopping at some later date
    except AssertionError:
      print '\n ERROR: multiple different repeat spacings have max copy number.\n'

  MaxNumberRepeatStarts = [MaxNumberStart] + EqualLengthStarts

  return MaxNumberRepeatStarts, TandemIndenticalSpacings


def vector_magnitude(Vector):
   ''' ND pythagorean theorem '''
   return ( np.sum( [ Component**2 for Component in Vector ] ) )**0.5

def derosettafy(RosettaVector):
  print ''' STOP USING derosettafy USE LIST '''
  assert 0 ,''' STOP USING derosettafy USE LIST '''
  List = [float(Value) for Value in RosettaVector]
  return np.array(List)

def angle(Vector1, Vector2):
  return np.arccos( np.dot(Vector1, Vector2) / ( vector_magnitude(Vector1) * vector_magnitude(Vector2) ) )

def rmsd_2_np_arrays(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array
  Thanks Daniel! """
  ##Corrected to account for removal of the COM
  COM1 = np.sum(crds1,axis=0) / crds1.shape[0]
  COM2 = np.sum(crds2,axis=0) / crds2.shape[0]
  crds1-=COM1
  crds2-=COM2
  n_vec = np.shape(crds1)[0]
  correlation_matrix = np.dot(np.transpose(crds1), crds2)
  v, s, w_tr = np.linalg.svd(correlation_matrix)
  is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
  if is_reflection:
          s[-1] = - s[-1]
          v[:,-1] = -v[:,-1]
  E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  #Calculate the rotation matrix
  rMtx=np.dot(v, w_tr)
  #Calculate the translation Vector
  tVec=COM1-(np.dot(COM2, np.linalg.inv(rMtx)))

  return np.sqrt(rmsd_sq), rMtx, tVec

def rmsd_2_np_arrays_rosetta(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array"""
  #D assert(crds1.shape[1] == 3)
  #D assert(crds1.shape == crds2.shape)

  ##Corrected to account for removal of the COM
  COM1 = np.sum(crds1,axis=0) / crds1.shape[0]
  COM2 = np.sum(crds2,axis=0) / crds2.shape[0]
  crds1-=COM1
  crds2-=COM2
  n_vec = np.shape(crds1)[0]
  correlation_matrix = np.dot(np.transpose(crds1), crds2)
  v, s, w_tr = np.linalg.svd(correlation_matrix)
  is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
  if is_reflection:
          s[-1] = - s[-1]
          v[:,-1] = -v[:,-1]
  E0 = sum(sum(crds1 * crds1)) + sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  #Calculate the rotation matrix
  rMtx=np.dot(v, w_tr)
  #Calculate the translation Vector
  tVec=COM1-(np.dot(COM2, np.linalg.inv(rMtx)))
  ###Are you kidding me??? Is this the correct way to build arrays inside 
  #of Rosetta core? No indexing?? Who was the creator of this???
  rMtx_xyzM=rosetta.numeric.xyzMatrix_double()
  rMtx_xyzM.xx=rMtx[0,0]
  rMtx_xyzM.xy=rMtx[0,1]
  rMtx_xyzM.xz=rMtx[0,2]
  rMtx_xyzM.yx=rMtx[1,0]
  rMtx_xyzM.yy=rMtx[1,1]
  rMtx_xyzM.yz=rMtx[1,2]
  rMtx_xyzM.zx=rMtx[2,0]
  rMtx_xyzM.zy=rMtx[2,1]
  rMtx_xyzM.zz=rMtx[2,2]
  tVec_xyzV=rosetta.numeric.xyzVector_double()
  tVec_xyzV.x=tVec[0]
  tVec_xyzV.y=tVec[1]
  tVec_xyzV.z=tVec[2]
  return np.sqrt(rmsd_sq), rMtx_xyzM, tVec_xyzV

def parse_motif_pose_coords(Pose):
   # Coordinates keyed with position 
   CoordHash = {}
   for Position in xrange( 1, Pose.n_residue()+1 ):
        # de-rosettafying things for using scipy.spatial later
      CoordHash[Position] = np.array( list(Pose.residue(Position).xyz('CA')) )

   # makes array from hash. seems stupid, will break if gaps in sequence, shouldn't be a problem if you use a pose
   CoordArray = np.array([ CoordHash[Position] for Position in xrange( 1, Pose.n_residue()+1 ) ])

   return CoordArray, CoordHash

def match_superimposed_pose_residues(Pose1, Pose2, MaxDist=1.0):
  ''' uses scipy.spatial.KDTree to find motif's w/in 1.0 of docked pose positions 
  returns list or array of lists
     If x is a single point, returns a list of the indices of the neighbors of x. 
     If x is an array of points, returns an object array of shape tuple containing lists of neighbors. " '''
  Pose1CoordArray, Pose1CoordHash = parse_motif_pose_coords(Pose1)
  Pose2Array, Pose2Hash = parse_motif_pose_coords(Pose2)

  Pose2Tree = spatial.KDTree(Pose2Array)
  # KDtree documentation is from http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html#scipy.spatial.KDTree.query_ball_point

  # KDTree.query_ball_point(x, r, p = 2.0, eps = 0)
  # Find all points within distance r of point(s) x.
  #
  # Parameters:   
  # x : array_like, shape tuple + (self.m,)
  # The point or points to search for neighbors of.
  # r : positive float
  #
  # Returns the radius of points to return.
  Pose1_Coords_by_Pose2_Coords = Pose2Tree.query_ball_point(Pose1CoordArray, MaxDist)
  # from 0 indexing to 1 indexing
  Pose1_Coords_by_Pose2_Coords = { PositionIndex+1 :[ Index+1 for Index in PositionResidues ] for PositionIndex, PositionResidues in enumerate(Pose1_Coords_by_Pose2_Coords) }

  return Pose1_Coords_by_Pose2_Coords


  def pdb_coord_array(self, Pdb):

    # gets full path to pdb in lab database
    PdbfullPath = ''.join( ['/lab/databases/pdb_clean/', Pdb[1:3].lower(), '/', Pdb, '.pdb'] )

    CalphaLines = subprocess.check_output(['grep', 'ATOM.........CA......%s'%Chain, PdbfullPath]).strip('\n').split('\n')
    CalphaValues = [ [ int(Line[22:26]), float(Line[30:38]), float(Line[38:46]), float(Line[46:54]) ] for Line in CalphaLines]

    self.ResidueIndentityDict = { int(Line[22:26]) : ThreeToOne[Line[17:20]] for Line in CalphaLines}
    return np.array( CalphaValues )


def main(argv=None):
  print 'no main function'

if __name__ == "__main__":
  sys.exit(main())