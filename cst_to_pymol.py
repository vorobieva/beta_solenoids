#! /usr/bin/env python

# 8/14/14
# Baker Lab UW

# # From Will's header from pymol_util
# #####################################
# import sys,os,inspect
# newpath = os.path.dirname(inspect.getfile(inspect.currentframe())) # script directory
# if not newpath in sys.path:
#   sys.path.append(newpath)
# # end of elements from Will's header 

# print 'newpath','/'.join(newpath.split('/')[:-1])
# PATH_TO_REPO = '/'.join(newpath.split('/')[:-1])
# sys.path.append(PATH_TO_REPO)

import argparse
import sys

def main(argv = None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' cst_to_pymol.py -cst contstraint_file.cst ')
  ArgParser.add_argument('-cst', type=str, nargs='+', help=' input cst file ', required=True)
  ArgParser.add_argument('-color', type=str, help=' color all cst one color ', default=None)
  ArgParser.add_argument('-multi', type=bool, help=' for viewing multiple cst in pymol ', default=False)
  ArgParser.add_argument('-name', type=bool, help=' include file name in pymol name ', default=False)
  Args = ArgParser.parse_args()

  for cst in Args.cst:
    with open(cst, 'r') as CstInput:
      CstLines = CstInput.readlines()

    NameTag = (cst+'!').replace('.cst!', '').replace('!', '')
    if Args.name or Args.multi:
      PymolName = '_'+NameTag[:]
    else:
      PymolName = ''

    for Index, Line in enumerate(CstLines):
      if Line.startswith('AtomPair'):
        CstType, Name1, Res1, Name2, Res2 = tuple( Line.split()[:5] )
        if not Args.multi:
          print 'dist %s, name %s and resi %s, name %s and resi %s'%('cst_%d_dist%s'%(Index + 1, PymolName), Name1, Res1, Name2, Res2)
        else:
          print 'dist %s, name %s and resi %s and %s, name %s and resi %s and %s'%('cst_%d_dist%s'%(Index + 1, PymolName), Name1, Res1, NameTag, Name2, Res2, NameTag)

      elif Line.startswith('Angle'):
        CstType, Name1, Res1, Name2, Res2, Name3, Res3 = tuple( Line.split()[:7] )
        if not Args.multi:
          print 'angle %s, name %s and resi %s, name %s and resi %s, name %s and resi %s'%('cst_%d_angle%s'%(Index + 1, PymolName), Name1, Res1, Name2, Res2, Name3, Res3)
        else:
          print 'angle %s, name %s and resi %s and %s, name %s and resi %s and %s, name %s and resi %s and %s'%('cst_%d_angle%s'%(Index + 1, PymolName), Name1, Res1, NameTag, Name2, Res2, NameTag, Name3, Res3, NameTag)

      elif Line.startswith('Dihedral'):
        CstType, Name1, Res1, Name2, Res2, Name3, Res3, Name4, Res4 = tuple( Line.split()[:9] )
        if not Args.multi:
          print 'dihedral %s, name %s and resi %s, name %s and resi %s, name %s and resi %s, name %s and resi %s'%('cst_%d_tor%s'%(Index + 1, PymolName), Name1, Res1, Name2, Res2, Name3, Res3, Name4, Res4)
        else:
          print 'dihedral %s, name %s and resi %s and %s, name %s and resi %s and %s, name %s and resi %s, and %s name %s and resi %s and %s'%('cst_%d_tor%s'%(Index + 1, PymolName), Name1, Res1, NameTag, Name2, Res2, NameTag, Name3, Res3, NameTag, Name4, Res4, NameTag)

    print 'set dash_width, 1'
    print 'set dash_width, 2, cst*angle*'
    print 'set dash_width, 5, cst*dist*'

    if Args.color:
      print 'color %s, cst_*%s'%(Args.color, PymolName)
    else:
      print 'color red, cst*angle*'
      print 'color cyan, cst*dist*'

if __name__ == "__main__":
  sys.exit(main())
