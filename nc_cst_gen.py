#! /usr/bin/env python

InfoString = ''' 
This script is to generate ATOMPAIR constraints for Rosetta,
 by default for all nitrogens and oxygens
                    within 3 angstroms 
                    on non-neighbor residues 
                    within an input pose. 

'''

'''
# This is used to give buried N - C contacts more weight. Thanks Alex Ford!
from interface_fragment_matching.utility.analysis import AtomicSasaCalculator

import tools

from multiprocessing import Process
from scipy import spatial

import numpy as np
import subprocess
import argparse
import sys
import os

import rosetta
rosetta.init(extra_options = "-mute basic -mute core -mute protocols")
'''

ThreeToOne = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','MET':'M','PRO':'P','PHE':'F','TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','TYR':'Y','CYS':'C','CYD':'C','LYS':'K','ARG':'R','HIS':'H','ASP':'D','GLU':'E','STO':'*','UNK':'U'}
ChainAlphabetIndices = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10, 'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 'S':19, 'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26 }

def pymol_commands(Pdb, Repeat, ReportedRepeatCount):
  return 'fetch %s\tselect rep%d, resi %s'%( Pdb, ReportedRepeatCount, '+'.join([str(Res) for Res in Repeat]) )

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
      print '\n\n\n LOOKOUT ERROR: multiple different repeat spacings have max copy number. \n\n\n'

  MaxNumberRepeatStarts = [MaxNumberStart] + EqualLengthStarts
  return MaxNumberRepeatStarts, TandemIndenticalSpacings


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' nc_cst_gen.py arguments ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=list, help=' input pdbs ', required=True)
  ArgParser.add_argument('-output', type=str, help=' output directory ', required=True)
  # Optional arguments:
  ArgParser.add_argument('-oxy', type=str, default='O\w?\d?', help=' grep for oxygen atoms ')
  ArgParser.add_argument('-nit', type=str, default='N\w?\d?', help=' grep for nitrogen atoms ')
  ArgParser.add_argument('-max_dist', type=float, default=3.0, help=' distance between the oxygens and nitrogens ')
  ArgParser.add_argument('-min_seq_sep', type=int, default=3, help=' minimum seperation in primary sequece ')
  ArgParser.add_argument('-num_repeats', type=int, default=5, help=' number of repeats to extrapolate contacts for ')

  ArgParser.add_argument('-min_sasa',  type=float, default=0.0,  help=' floor for weighting nitrogen oxygen contacts ')
  ArgParser.add_argument('-max_sasa',  type=float, default=500.0,  help=' ceiling for cst weighting nitrogen oxygen contacts ')
  ArgParser.add_argument('-sasa_probe_radius', type=float, default=2.2,  help=' probe radius for sasa calculations ')
  args = ArgParser.parse_args()
  
  OxygenGrep = Args.oxy
  NitrogenGrep = Args.nit

  # for making full atom kd tree
  ResAtmCoordLists = []
  # for translating from kd tree index to ( residue, atom ) coord
  ResAtmRecordLists = []

  # loop through all residue numbers
  for Res in range(1, Pose.n_residue() + 1):
    # remade for each residue
    AtmRecordList = []
    AtmCoordList = []
    # loop through residue's atom numbers
    for Atm in range(1, Pose.residue(Res).natoms() + 1):
      # add (residue, atom) coord to residue's list
      AtmRecordList.append((Res, Atm))
      # add atom xyz coord to residue's list
      AtmCoordList.append(Pose.residue(Res).atom(Atm).xyz())
    # add residue's lists to respective global lists
    ResAtmRecordLists.append(AtmRecordList)
    ResAtmCoordLists.append(AtmCoordList)

  ResidueAtomArray = np.array( ResAtmCoordLists )
  ResidueAtomKDTree = spatial.KDTree( ResidueAtomArray )

  ReportedRepeatCount = 0
  TotalPdbs = len(Args.pdbs)

  for iPdb, Pdb in enumerate(Args.pdbs):
    print ' Working with %s; %d of %d total pdbs '%(Pdb, iPdb+1, TotalPdbs)

    PdbWDAG = pdb_wdag(Pdb, 3, 4.7, 0.25, 15 )
    PdbWDAG.find_downstream_neighbors()
    RepeatChains = PdbWDAG.find_repeat_chains()
    
    ConsolidatedRepeatStarts, TandemRepeats = consolidate_repeats(RepeatChains)

    # Starting rosetta  
    Pose = rosetta.pose_from_pdb(Pdb)

    # thanks Alex!!!!
    ResidueAtomSasa = SasaCalculator.calculate_per_atom_sasa(Pose)
    
    # holds constraints before printing
    Constraints = []
    Oxygens = []
    # Nitrogens = {}
    # Other = {}

    # loop through residues, checks for oxygens in these outer
    # this is not 'import this' like :(
    for OxyRes in range( 1, Pose.n_residue() + 1 ):
      # loop through atoms
      for OxyAtm in range( 1, Pose.residue(OxyRes).natoms() + 1 ):

        print Pose.residue(OxyRes).atom_names(OxyAtm)
        
        OxyName = Pose.residue(OxyRes).atom_names(OxyAtm)
        # checks oxygen name
        if re.match(OxygenGrep, OxyName ): # <- this guy 
          # gets the coordinates from 
          oxygen_xyz_coords = tools.derosettafy( Pose.residue(OxyRes).atom(OxyAtm).xyz() )
          Oxygens.append((OxyRes, OxyAtm, OxyName, oxygen_xyz_coords))

    # assembles array with coordinates of all oxygens
    All_Oxygen_Array = np.array([ Oxy[3] for Oxy in Oxygens ])
    # Query KDtree for atoms within -max_dist from oxygens
    Neighbors_of_Oxygens = ResidueAtomKDTree.query_ball_point( All_Oxygen_Array, Args.max_dist )
    
    # sort out results
    for o, Oxy in enumerate(Oxygens):
      # unpack Oxy tuple
      OxyRes, OxyAtm, OxyName, oxygen_xyz_coords = Oxy
      Neighbor_Set = Neighbors_of_Oxygens[o]
      # prep for loop
      Neighbor_Coordinates = []     
      Neighbor_Names = []
      for i in Neighbor_Set:
        Res, Atm = ResAtmRecordLists[i]
        Neighbor_Coordinates.append( tools.derosettafy( Pose.residue(Res).atom(Atm).xyz() ) )
        Neighbor_Names.append( Pose.residue(Res).atom_name(Atm) )

      # Query KDtree for atoms closer than -max_dist from oxygens
      Neighbors_of_Oxygen_Array = np.array(Neighbor_Coordinates)
      Nearer_Distance = Args.max_dist * 0.9
      Neighbors_of_Neighbors = ResidueAtomKDTree.query_ball_point( Neighbors_of_Oxygen_Array, Nearer_Distance )

      for o, Oxy in enumerate(Oxygens):
        Neighbors Neighbors_of_Neighbors[o]
      
      Constraints.append('AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %f SUMFUNC 2 HARMONIC %.2f 1.0 CONSTANTFUNC -0.5' %( OxyName, OxyRes, NitroName, NitroRes, Weight, Distance ))

    # Nitrogen
    #   # Make new array from KDtree hits
    #   near_xyz_coords = np.array([ ResidueAtomArray[hit] for hit in oxygen_neighbors ])
    #   # this is to capture covalent and/or ionic bond lengths
    #   # Neighbors of neighbors, for angle constraints 
    #   secondary_hits = ResidueAtomKDTree.query_ball_point( near_xyz_coords, Nearer_Distance )

    #   # Loops through neighbors to look for nitrogens
    #   for i, hit in enumerate(oxygen_neighbors):
    #     NitroRes, NitroAtm = AtmRecordList[hit]
    #     NitroName = Pose.residue(NitroRes).atom_names(NitroAtm)

    #     if re.match(NitrogenGrep, NitroName ): 
    #       NitroCoords = Pose.residue(NitroRes).xyz(NitroAtm)

    #       # gets the coordinates from 
    #       NitrogenArray = tools.derosettafy( Pose.residue(NitroRes).atom(NitroAtm).xyz() )
          
    #       # almost certainly some rosetta way of doing this
    #       Displacement = OxygenArray - NitrogenArray
    #       Distance = tools.vector_magnitude(Displacement)

    #       # 'AtomPair %s %d %s %d SCALARWEIGHTEDFUNC 2.0 SUMFUNC 2 SIGMOID %.2f %.2f CONSTANTFUNC -0.5' %( OxyName, OxyRes, NitroName, NitroRes, Distance, Slope )

    #       # Loop through atoms that are closer to the nitrogen than the oxygen is
    #       for i, hit in enumerate(secondary_hits):
    #         print hit


if __name__ == "__main__":
  sys.exit(main())