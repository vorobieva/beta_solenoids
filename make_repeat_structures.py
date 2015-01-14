#! /usr/bin/env python


InfoString = '''
1) generate_cst.py -pdbs 3fy3.pdb

1) generate_backbones.py -pdbs 3fy3.pdb -repeat 4 -max_turns_per_repeat 2

2) remove_redundant_structures.py -pdb_stem _3fy3

3) expand_cst.py -ref_pdb 3fy3.pdb -ref_cst 3fy3_All.cst -repeat_pdb_tag _3fy3
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
 
  ## debug
  # rosetta.init(extra_options = "-out:level 1000")

  ## runtime
  # rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
  from rosetta.protocols import grafting 
  from rosetta.utility import ostream

  # from repo 
  import solenoid_tools
  import generate_backbones
  import generate_cst
  import expand_cst
  import remove_redundant_structures
  
# ''' 

# sys.argv = [ sys.argv[0], '-ref_pdb', '2G0Y_Relax.pdb', '-ref_cst', '2G0Y_Relax_All.cst', '-repeat_tag', '_src22_41__35_54_rep40_2rwAsp_hfGlu']


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  # Arg block
  ArgParser = argparse.ArgumentParser(description=' make_repeat_structures.py ( -help ) %s'%InfoString)
  # Required args
  ArgParser.add_argument('-pdb', type=str, help=' reference pdb ', required=True)
  # Optional args
  ArgParser.add_argument('-repeat', type=str, help=' input pdb tag ', default='4')
  ArgParser.add_argument('-max_turns_per_repeat', type=str, help=' corresponding to reference pdb ', default='2')

  Args = ArgParser.parse_args()

  GenBbArgs = ['-pdbs', Args.pdb, '-repeat', Args.repeat, '-max_turns_per_repeat', Args.max_turns_per_repeat]
  print '\nRunning generate_backbones: ', GenBbArgs
  generate_backbones.main( GenBbArgs )

  GenCstArgs = ['-pdbs', Args.pdb]
  print '\nRunning generate_cst: ', GenCstArgs
  generate_cst.main( GenCstArgs )

  PdbBase = re.sub(r'^(.*).pdb$', r'\1', Args.pdb)

  RemoveRedunStrucArgs = ['-pdb_stem', '_%s'%PdbBase ]
  print '\nRunning remove_redundant_structures: ', RemoveRedunStrucArgs
  remove_redundant_structures.main( RemoveRedunStrucArgs )

  ExpandCstArgs = ['-ref_pdb', Args.pdb, '-ref_cst', '%s_All.cst'%PdbBase, '-repeat_pdb_tag', '_%s'%PdbBase]
  print '\nRunning expand_cst: ', ExpandCstArgs
  expand_cst.main( ExpandCstArgs )



if __name__ == "__main__":
  sys.exit(main())
