#! /usr/bin/env python

from scipy import spatial
import numpy as np

def vector_magnitude(Vector):
   ''' ND pythagorean theorem '''
   return ( np.sum( [ Component**2 for Component in Vector ] ) )**0.5

def derosettafy(RosettaVector):
  List = [float(Value) for Value in RosettaVector]
  return np.array(List)

def angle(Vector1, Vector2):
  return np.arccos( np.dot(Vector1, Vector2) / ( vector_magnitude(Vector1) * vector_magnitude(Vector2) ) )

def rmsd_2_np_arrays(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] numpy array
  Thanks Daniel! """
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

  return np.sqrt(rmsd_sq), rMtx, tVec

def parse_motif_pose_coords(Pose):
   # Coordinates keyed with position 
   CoordHash = {}
   for Position in xrange( 1, Pose.n_residue()+1 ):
        # de-rosettafying things for using scipy.spatial later
      CoordHash[Position] = np.array( [ Number for Number in Pose.residue(Position).xyz('CA') ] )

   # makes array from hash. seems stupid, will break if gaps in sequence, shouldn't be a problem if you use a pose
   CoordArray = np.array([ CoordHash[Position] for Position in xrange( 1, Pose.n_residue()+1 ) ])

   return CoordArray, CoordHash

def match_hotspots_and_motifs_to_dock(DockedPose, MotifPose):
   ''' uses scipy.spatial.KDTree to find motif's w/in 1.0 of docked pose positions 
   returns list or array of lists
       If x is a single point, returns a list of the indices of the neighbors of x. 
       If x is an array of points, returns an object array of shape tuple containing lists of neighbors. " '''

   DockCoordArray, DockCoordHash = parse_motif_pose_coords(DockedPose)
   MotifCoordArray, MotifCoordHash = parse_motif_pose_coords(MotifPose)
 
   # print MotifCoordArray

   MotifCoordTree = spatial.KDTree(MotifCoordArray)
   # KDtree documentation is from http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.KDTree.query_ball_point.html#scipy.spatial.KDTree.query_ball_point

   # KDTree.query_ball_point(x, r, p=2.0, eps=0)[source]
   # Find all points within distance r of point(s) x.
   # Parameters:   
   # x : array_like, shape tuple + (self.m,)
   # The point or points to search for neighbors of.
   # r : positive float
   # The radius of points to return.
   MotifsByPosition = MotifCoordTree.query_ball_point(DockCoordArray, 1.0)

   # from 0 indexing to 1 indexing
   MotifsByPosition = { PositionIndex+1 :[ Index+1 for Index in PositionMotifs ] for PositionIndex, PositionMotifs in enumerate(MotifsByPosition) }

   # print MotifsByPosition

   return MotifsByPosition
