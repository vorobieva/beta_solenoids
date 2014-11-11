#! /usr/bin/env python

import argparse
import sys

def main(argv = None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' cst_to_pymol.py -cst contstraint_file.cst ')
  ArgParser.add_argument('-cst', type=str, help=' input cst file ', required=True)
  Args = ArgParser.parse_args()

  with open(Args.cst, 'r') as CstInput:
    for Index, Line in enumerate(CstInput.readlines()):
      if Line.startswith('AtomPair'):
        CstType, Name1, Res1, Name2, Res2 = tuple( Line.split()[:5] )
        print 'dist %s, name %s and resi %s, name %s and resi %s'%('Line_%d'%(Index+1), Name1, Res1, Name2, Res2)
      
      elif Line.startswith('Angle'):
        CstType, Name1, Res1, Name2, Res2, Name3, Res3 = tuple( Line.split()[:7] )
        print 'angle %s, name %s and resi %s, name %s and resi %s, name %s and resi %s'%('Line_%d'%(Index+1), Name1, Res1, Name2, Res2, Name3, Res3)


if __name__ == "__main__":
  sys.exit(main())