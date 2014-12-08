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

# sys.argv = [ sys.argv[0], '-pdb_stem', '_1M8N_Relax', '-out', 'Unrun', '-run_tag', 'TalCst0.200_Relax']

def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=" args for optimize_repeat_structures ")
  ArgParser.add_argument('-pdb_stem', type=str, help=" pdb stem, start of globs for pdbs and csts ", required=True )    
  ArgParser.add_argument('-out', type=str, help=" folder to move files to to ", required=True )    
  ArgParser.add_argument('-run_tag', type=str, help=" consider jobs unfinished if run with tag is not present ", default=None )
  Args = ArgParser.parse_args()
  
  Pdbs = glob.glob('*%s.pdb'%Args.pdb_stem)

  if not os.path.isdir(Args.out):
    subprocess.check_output(['mkdir', Args.out])
  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  for Pdb in Pdbs:
    
    PdbStem = re.sub(r'(.*).pdb$', r'\1', Pdb)
    MovedPdb = False
    GlobString = '%s*'%PdbStem
    PdbAssociatedFiles = glob.glob(GlobString)
    PdbCsts = [ File for File in PdbAssociatedFiles if File.endswith('cst') ]
    AllRelevantPdbs = [ File for File in PdbAssociatedFiles if File.endswith('pdb') ] 

    for Cst in PdbCsts:
      Stem = re.sub(r'(.*).cst$', r'\1', Cst)
      RunsWithCst = [ File for File in PdbAssociatedFiles if File.endswith('pdb') ]
      CstRelevantPdbs = [ File for File in AllRelevantPdbs if File.startswith(Stem) ]

      if Args.run_tag:
        found = 0
        for PdbFileName in CstRelevantPdbs:
          if Args.run_tag in PdbFileName:
            found = 1

      else:
        if len(CstRelevantPdbs):
          found = 1
        else:
          found = 0

      if found:
        print '%s ran with %s already'%( Pdb, Cst )
      else:
        print '\t %s not yet run with %s'%( Pdb, Cst )
        subprocess.check_output( ['mv', Cst, '%s'%Args.out] )
        if not MovedPdb:
          subprocess.check_output( ['mv', Pdb, '%s'%Args.out] )
          MovedPdb = True

if __name__ == "__main__":
  sys.exit(main())

