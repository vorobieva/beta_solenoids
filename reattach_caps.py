#! /usr/bin/env python
InfoString = ''' 

Run on expanded pose with src range tag in pdb name like:
.*src13_26__15_28.*.pdb

'''

'''
# from repo 
import solenoid_tools

from generate_backbones import fuse

# external
# libraries
#
# from scipy import spatial
# import itertools
import numpy as np
import subprocess
import argparse
import glob
import sys
import os
import re

import rosetta
# rosetta.init()
rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
from rosetta.protocols import grafting 

# '''     # # get cst file if csts files
    # if Args.csts:
    #   Cst = Args.csts[i]
    #   print 'Cst:', Cst
    #   Constrainer = constraint_extrapolator(Cst)
def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' reattach_cap.py ( -help ) %s'%InfoString)
  # Required args
  ArgParser.add_argument('-pdb_stem', type=str, nargs='+', help=' Input pdbs ', required=True)
  ArgParser.add_argument('-reference_pdb', type=str, nargs='+', help=' Reference pdb ', required=True)
  ArgParser.add_argument('-repeat_pdb_tag', type=str, nargs='+', help=' Input pdb tag ', required=True)
  # Optional args
  ArgParser.add_argument('-out', type=str, help=' Output directory ', default='./')
  Args = ArgParser.parse_args()

  Pdbs = glob.glob( '*%s*.pdb'%Args.repeat_pdb_tag ) 

  for Pdb in Pdbs:
    print Pdb 

    SourceRanges = re.sub(r'.*src(\d+_\d+__\d+_\d+).*pdb', r'\1', Pdb)

    SourceRanges[0]

  if Args.out [-1] != '/':
    Args.out = Args.out + '/'


  # EndMatch.reverse()
  # EndMatchArray = get_residue_array(Repeat2Pose, EndMatch)

  # # print 'EndMatch', EndMatch
  # RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays_rosetta( EndMatchArray, NterminalUnmatchedArray )
  # # print 'Corse RMSD:  ', RMSD
  # # print 'rMtx:  ', rMtx
  # # print 'tVec:  ', tVec

  # DuplicatePose = rosetta.Pose()
  # DuplicatePose.assign(FusionBasePose)
  # rosetta.Pose.apply_transform_Rx_plus_v(DuplicatePose, rMtx, tVec)
  # # rosetta.dump_pdb(DuplicatePose, 'Dup%d_DuplicatePose.pdb'%Duplication)

  # FusionPose, RMSD, CorrespondingResidues = fuse(FusionPose, DuplicatePose)
  # # print 'Refined RMSD', RMSD
  # # rosetta.dump_pdb(FusionPose, 'Dup%d_FusionPose.pdb'%Duplication)
  

# if __name__ == "__main__":
#   sys.exit(main())

