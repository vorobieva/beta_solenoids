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
  rosetta.init(extra_options = " -ex1 -ex2 -no_optH false -use_input_sc ") # "-mute basic -mute core -mute protocols"
  from rosetta.protocols import grafting 
  
  # from repo 
  import solenoid_tools

# '''

# sys.argv = [ sys.argv[0], '-pdb_stem', '_2qiv_Relax' ]

def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=" args for optimize_repeat_structures ")
  ArgParser.add_argument('-pdb_stem', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  Args = ArgParser.parse_args()
  
  Pdbs = glob.glob('*%s.pdb'%Args.pdb_stem)
  
  PdbSortTuples = []
  Skipped = []

  for Pdb in Pdbs:
    RepeatLength = int(re.sub(r'.*rep(\d+).*pdb', r'\1', Pdb))
    SourceStart = int(re.sub(r'.*src(\d+).*pdb', r'\1', Pdb))
    
    try:
      assert SourceStart != Pdb and RepeatLength != Pdb, 'regular expression substitution failed' 
    except AssertionError:
      Skipped.append(Pdb)
      continue

    PdbSortTuples.append( (RepeatLength, SourceStart, Pdb) )
  
  print 'Skipped:'
  print Skipped
  print
  
  PdbSortTuples.sort()

  LastPdb = PdbSortTuples[0][2]
  Pose = rosetta.pose_from_pdb(LastPdb)
  LastArray = np.array([ list( Pose.residue(P).xyz('CA') ) for P in range(1, Pose.n_residue()+1) ])

  subprocess.check_output(['mkdir', 'Redundant'])

  for PdbTup in PdbSortTuples[1:]:
    Pdb = PdbTup[2]
    Pose = rosetta.pose_from_pdb(Pdb)

    CA_Array = np.array([ list( Pose.residue(P).xyz('CA') ) for P in range(1, Pose.n_residue()+1) ])
  
    if len(CA_Array) == len(LastArray):

      RMSD, rMtx, tVec = solenoid_tools.rmsd_2_np_arrays(CA_Array, LastArray)
      print
      print 'LastPdb, Pdb'
      print LastPdb
      print Pdb
      print 'RMSD:', RMSD

      if RMSD < 0.001:
        PdbStem = re.sub(r'(.*).pdb$', r'\1', Pdb)
        GlobString = '%s*'%PdbStem

        PdbAssociatedFiles = glob.glob(GlobString)
        # print PdbAssociatedFiles

        for File in PdbAssociatedFiles:
          subprocess.check_output(['mv', File, 'Redundant/'])


    LastArray = copy.deepcopy(CA_Array)
    LastPdb = copy.deepcopy(Pdb)


if __name__ == "__main__":
  sys.exit(main())

