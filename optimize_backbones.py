#! /usr/bin/env python
InfoString = ''' 
TO GENERATE BACKBONES AROUND STARTING STRUCTURE 
    WHICH MAY BE MODIFIED 
BY STRECHING REPEAT REGION DESIGNATED BY NUMBER
'''

'''
# from repo 
import solenoid_tools

# libraries
from scipy import spatial
import multiprocessing
import numpy as np
import subprocess
import argparse
import glob
import sys
import os
import re

import rosetta
rosetta.init()
# rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
from rosetta.protocols import grafting 

from rosetta.core.pose import initialize_atomid_map
# from rosetta.core.id import AtomID, AtomID_Map_Real, AtomID_Map_bool

# pymol_link = rosetta.PyMolMover()

sys.argv.extend(['-pdbs', '1EZG.pdb', '-out', './' ])
# '''


def mcmc_mover(Pose, ScoreFunction, kT=1.0, Gen=100):
  ''' threadable backbone breathing and side chain repacking with constraints '''
  print 'Starting work on', Pose

  movemap = rosetta.MoveMap()
  movemap.set_bb(True)
  small_mover = rosetta.SmallMover(movemap, kT, 1)
  shear_mover = rosetta.ShearMover(movemap, kT, 1)

  MontyCarlos = rosetta.MonteCarlo(Pose, ScoreFunction, kT)

  small_mover(Pose)

  for g in range(Gen):
    small_mover(Pose)
    MontyCarlos.boltzmann(Pose)
    shear_mover(Pose)
    MontyCarlos.boltzmann(Pose)


  print MontyCarlos.show_scores()
  print MontyCarlos.show_counters()
  print MontyCarlos.show_state()



def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' optimize_backbone.py arguments ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=list, help=' input pdbs ', required=True)
  ArgParser.add_argument('-start', type=list, help=' repeat start residue ', default=False)
  # ArgParser.add_argument('-end', type=list, help=' input pdbs ', required=True)
  ArgParser.add_argument('-repeat', type=list, help=' number of repeats to make ', default=3)
  ArgParser.add_argument('-dat', type=str, help=' input symm dat file ', default=None)
  ArgParser.add_argument('-cst', type=str, help=' input symm cst file ', default=None)
  ArgParser.add_argument('-out', type=str, help=' output directory ', default='./')
  Args = ArgParser.parse_args()
  if len(Args.pdbs[0]) == 1:
    Args.pdbs = [''.join(Args.pdbs)]

'''
  # return Args
  # Super list to hold multithreaded design trajectories from all input pdbs
  PdbJobs = []

  # ScoreFunction = rosetta.create_score_function("standard")
  # ScoreFunction.set_weight(rosetta.atom_pair_constraint, 5.0)

      # Rosetta part
    # mcmc_mover(Pose)

    # Sub list to hold multithreaded design trajectories from all frames within current itertaions pdb
    # FrameJobs = []
    # if Arg.start:
    #   Arg.start
    # else:
      # for Frame in Frames:
      #   SingleProcess = multiprocessing.Process(target=worker, args=(i,))
      #   FrameJobs.append(SingleProcess)
      #   SingleProcess.start()

    # PdbJobs.append(FrameJobs)
'''

# if __name__ == "__main__":
#   sys.exit(main())