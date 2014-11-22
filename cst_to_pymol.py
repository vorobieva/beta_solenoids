#! /usr/bin/env python

# 8/14/14
# Baker Lab UW

# From Will's header from pymol_util
#####################################
import sys,os,inspect
newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
if not newpath in sys.path:
  sys.path.append(newpath)
# end of elements from Will's header 

# print 'newpath','/'.join(newpath.split('/')[:-1])

PATH_TO_REPO = '/'.join(newpath.split('/')[:-1])
sys.path.append(PATH_TO_REPO)

import argparse
import sys

def main(argv = None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' cst_to_pymol.py -cst contstraint_file.cst ')
  ArgParser.add_argument('-cst', type=str, help=' input cst file ', required=True)
  Args = ArgParser.parse_args()

  with open(Args.cst, 'r') as CstInput:
    CstLines = CstInput.readlines()

  for Index, Line in enumerate(CstLines):
    if Line.startswith('AtomPair'):
      CstType, Name1, Res1, Name2, Res2 = tuple( Line.split()[:5] )
      print 'dist %s, name %s and resi %s, name %s and resi %s'%('cst_%d_dist'%(Index+1), Name1, Res1, Name2, Res2)
    
    elif Line.startswith('Angle'):
      CstType, Name1, Res1, Name2, Res2, Name3, Res3 = tuple( Line.split()[:7] )
      print 'angle %s, name %s and resi %s, name %s and resi %s, name %s and resi %s'%('cst_%d_angle'%(Index+1), Name1, Res1, Name2, Res2, Name3, Res3)

    elif Line.startswith('Dihedral'):
      CstType, Name1, Res1, Name2, Res2, Name3, Res3, Name4, Res4 = tuple( Line.split()[:9] )
      print 'dihedral %s, name %s and resi %s, name %s and resi %s, name %s and resi %s, name %s and resi %s'%('cst_%d_tor'%(Index+1), Name1, Res1, Name2, Res2, Name3, Res3, Name4, Res4)

  print 'set dash_width, 1'
  print 'color red, cst*angle'
  print 'set dash_width, 2, cst*angle'
  print 'color cyan, cst*dist'
  print 'set dash_width, 5, cst*dist'

if __name__ == "__main__":
  sys.exit(main())