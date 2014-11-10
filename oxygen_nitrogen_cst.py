#! /usr/bin/env python

InfoString = ''' 
This script is to generate ATOMPAIR constraints for Rosetta,
 by default for all nitrogens and oxygens
                    within 3 angstroms 
                    on non-neighbor residues 
                    within an input pose. 

'''

# comment just next line and copy block in 
# multiline string for ipython mode

# '''
# This is used to give buried N - O contacts more weight. Thanks Alex Ford!
try:
  from interface_fragment_matching.utility.analysis import AtomicSasaCalculator
except ImportError:
  ' Error: SASA weighting of contacts requires interface_fragment_matching from Alex Ford '

import tools

from multiprocessing import Process
from scipy import spatial
import numpy as np
import subprocess
import argparse
import sys
import os
import re

import rosetta
rosetta.init(extra_options = "-mute basic -mute core -mute protocols")

# ''', sys.argv.extend(['-pdbs', '1EZG.pdb', '-out', './' ])

ThreeToOne = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','MET':'M','PRO':'P','PHE':'F','TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','TYR':'Y','CYS':'C','CYD':'C','LYS':'K','ARG':'R','HIS':'H','ASP':'D','GLU':'E','STO':'*','UNK':'U'}
ChainAlphabetIndices = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10, 'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 'S':19, 'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26 }

def pymol_commands(Pdb, Repeat, ReportedRepeatCount):
  return 'fetch %s\tselect rep%d, resi %s'%( Pdb, ReportedRepeatCount, '+'.join([str(Res) for Res in Repeat]) )

class sasa_scale:
  ''' gives weight in arbitray range based on sasa score '''
  def __init__(self, FloorSasa, FloorWeight, CeilingSasa, CeilingWeight):
    self.FloorSasa = FloorSasa
    self.FloorWeight = FloorWeight
    self.CeilingSasa = CeilingSasa
    self.CeilingWeight = CeilingWeight
    
    self.SasaRangeMag = CeilingSasa - FloorSasa
    self.WeightRangeMag = CeilingWeight - FloorWeight

  def weigh(self, SasaToWeigh):
    if SasaToWeigh <= self.FloorSasa:
      return self.FloorWeight
    elif SasaToWeigh >= self.CeilingSasa:
      return self.CeilingWeight
    else:
      Proportion =  (SasaToWeigh - self.FloorSasa) / self.SasaRangeMag
      return (self.WeightRangeMag * Proportion) + self.FloorWeight


def consolidate_repeats(ListOfRepeatPositions):
  ''' takes lists of repeat positions (based on xyz coords) and returns subsets that split the pose longest equal residue number repeat chains '''  
  TandemIndenticalSpacings = {} #{0:[]}
  
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

def get_pose_nc_constraints(Pose, SasaScale):
    '''  '''
    # SasaCalculator is from Alex's interface_fragment_matching 
    # thanks Alex!
    #
    try:
      ResidueAtomSasa = SasaCalculator.calculate_per_atom_sasa(Pose)
    except:
      ResidueAtomSasa = SasaCalculator.calculate_per_atom_sasa(Pose)

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
        AtmCoordList.append( list(Pose.residue(Res).atom(Atm).xyz()) )
      
      # add residue's lists to respective global lists
      ResAtmRecordLists.extend(AtmRecordList)
      ResAtmCoordLists.extend(AtmCoordList)

    ResidueAtomArray = np.array( ResAtmCoordLists )
    ResidueAtomKDTree = spatial.KDTree( ResidueAtomArray )

    Oxygens = []

    # loop through residues storing info on oxygens
    for OxyRes in range( 1, Pose.n_residue() + 1 ):
      # loop through atoms
      for OxyAtm in range( 1, Pose.residue(OxyRes).natoms() + 1 ):
        OxyName = Pose.residue(OxyRes).atom_name(OxyAtm).replace(' ', '')

        #                   this code makes me feel like this guy 
        #                               | 
        # checks oxygen name            v
        if re.match(OxygenGrep, OxyName ):
          # print 'found oxygen %s'%OxyName
          # gets the coordinates from 
          OxyXyzCoords = np.array( list( Pose.residue(OxyRes).atom(OxyAtm).xyz ) )
          Oxygens.append( (OxyRes, OxyAtm, OxyName, OxyXyzCoords) )

    # assembles array with coordinates of all oxygens
    All_Oxygen_Array = np.array([ Oxy[3] for Oxy in Oxygens ])
    # Query KDtree for atoms within -max_dist from oxygens
    Neighbors_of_Oxygens = ResidueAtomKDTree.query_ball_point( All_Oxygen_Array, Args.max_dist )

    # holds constraints before printing
    AllConstraints = []    
    # ConstraintsByResidue = []

    # Loop through oxygens
    for o, Oxy in enumerate(Oxygens):
      # unpack Oxy tuple
      OxyRes, OxyAtm, OxyName, OxyXyzCoords = Oxy

      Neighbor_Set = Neighbors_of_Oxygens[o]
      # print 'Neighbor number: ', len(Neighbor_Set)
      # prep for loop
      Nitrogens = []

      for i in Neighbor_Set:
        NitroRes, NitroAtm = ResAtmRecordLists[i]
        NitroName = Pose.residue(NitroRes).atom_name(NitroAtm).replace(' ', '')
        # print NitroName
        # checks nitrogen name
        if re.match( NitrogenGrep, NitroName ):
          # checks primary sequence spacing
          # to remove set:  -min_seq_sep 0
          if np.abs(OxyRes - NitroRes) >= Args.min_seq_sep:
            # print 'found nitrogen neighbor %s'%NitroName
            NitroXyzCoords = list(Pose.residue(NitroRes).atom(NitroAtm).xyz())
            # print 'NitroRes, NitroName', NitroRes, NitroName
            # print 'NitroXyzCoords', NitroXyzCoords
            Nitrogens.append((NitroRes, NitroAtm, NitroName, NitroXyzCoords))
        elif 'N' in NitroName:
          print 'skipping atom named: %s , make sure it is not a nitrogen'%NitroName

      Distance = tools.vector_magnitude(NitroXyzCoords - OxyXyzCoords)

      if len(Nitrogens):
        Neighbor_Coordinates = [Nitro[3] for Nitro in Nitrogens]
        # print 'Neighbor_Coordinates', Neighbor_Coordinates[0:4], Neighbor_Coordinates[-4:]
        # Query KDtree for atoms closer than -max_dist from oxygens
        Neighbors_of_Oxygen_Array = np.array(Neighbor_Coordinates)
        Nearer_Distance = Distance * 0.9
        # print 'Neighbors_of_Oxygen_Array', Neighbors_of_Oxygen_Array
        Neighbors_of_Nitrogens = ResidueAtomKDTree.query_ball_point( Neighbors_of_Oxygen_Array, Nearer_Distance )
      
      Constraints = []
      #  Loop through all nitrogen neighbors of oxygen
      for n, Nitro in enumerate(Nitrogens):
        # unpack Nitro tuple
        NitroRes, NitroAtm, NitroName, NitroXyzCoords = Nitro

        # these trys / excepts seperate 
        # backbone-backbone from 
        # backbone-sidechain from
        # sidechain-sidechain interactions
        # 
        # in future maybe sort into seperate lists
        try:
          OxygenSasa = ResidueAtomSasa[OxyRes][OxyName]
          NitrogenSasa = ResidueAtomSasa[NitroRes][NitroName]
          AverageSasa = np.mean([OxygenSasa, NitrogenSasa])        
        except KeyError:
          # These lines handle backbone to sidechain interactions
          # set weight equal to the most buried 
          try:
            OxygenSasa = ResidueAtomSasa[OxyRes][OxyName]
            AverageSasa = SasaScale.FloorSasa
          except KeyError:
            try:
              NitrogenSasa = ResidueAtomSasa[NitroRes][NitroName]
              AverageSasa = SasaScale.FloorSasa 
            
            # set weight of side chain side chain equal to the most buried             
            except KeyError:
              AverageSasa = SasaScale.CeilingSasa 

        # Neighbors_of_Nitrogens[n]

        # use instance of sasa_scale to calculate weight based on avg sasa of N and O
        SasaBasedWeight = SasaScale.weigh(AverageSasa)
        # print 
        # print 'AverageSasa', AverageSasa
        # print 'SasaBasedWeight', SasaBasedWeight
        # adds atom pair constraint to list of constraints
        Constraints.append('AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %f SUMFUNC 2 HARMONIC %.2f 1.0 CONSTANTFUNC -0.5' %( OxyName, OxyRes, NitroName, NitroRes, SasaBasedWeight, Distance ))
      
      AllConstraints.extend(Constraints)
      # ConstraintsByResidue.append(Constraints)

    return AllConstraints


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' nc_cst_gen.py arguments ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=list, help=' input pdbs ', required=True)
  ArgParser.add_argument('-out', type=str, help=' output directory ', required=True)
  # Optional arguments:
  ArgParser.add_argument('-max_dist', type=float, default=3.2, help=' distance between the oxygens and nitrogens ')
  ArgParser.add_argument('-min_seq_sep', type=int, default=3, help=' minimum seperation in primary sequece ')
  ArgParser.add_argument('-oxy', type=str, default='O\w?\d?', help=' grep for oxygen atoms ')
  ArgParser.add_argument('-nitro', type=str, default='N\w?\d?', help=' grep for nitrogen atoms ')
  ArgParser.add_argument('-num_repeats', type=int, default=5, help=' number of repeats to extrapolate contacts for ')
  ArgParser.add_argument('-min_sasa',  type=float, default=0.0,  help=' floor for weighting nitrogen oxygen contacts ')
  ArgParser.add_argument('-min_sasa_weight',  type=float, default=1.0,  help=' weight of floor for nitrogen oxygen contacts ')
  ArgParser.add_argument('-max_sasa',  type=float, default=5.0,  help=' ceiling for cst weighting nitrogen oxygen contacts ')
  ArgParser.add_argument('-max_sasa_weight',  type=float, default=0.1,  help=' weight of ceiling for nitrogen oxygen contacts ')
  ArgParser.add_argument('-sasa_probe_radius', type=float, default=0.8,  help=' probe radius for sasa calculations ')
  ArgParser.add_argument('-renumber_pose', type=bool, default=True, help='True|False renumber pdb residues ' )
  Args = ArgParser.parse_args()
  
  if len(Args.pdbs[0]) == 1:
    Args.pdbs = [''.join(Args.pdbs)]

  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  OxygenGrep = Args.oxy
  NitrogenGrep = Args.nitro

  ReportedRepeatCount = 0
  TotalPdbs = len(Args.pdbs)

  # Instance of class to convert sasas to cst weight
  SasaScale = sasa_scale( Args.min_sasa, Args.min_sasa_weight, Args.max_sasa, Args.max_sasa_weight )
                       #(FloorSasa, FloorWeight, CeilingSasa, CeilingWeight)

  for iPdb, Pdb in enumerate(Args.pdbs):
    print ' Working with %s; %d of %d total pdbs '%(Pdb, iPdb+1, TotalPdbs)

    # Starting rosetta  
    Pose = rosetta.pose_from_pdb(Pdb)

    # Sets pdb info so residues in dumped pdbs are same as index 
    Pose.pdb_info(rosetta.core.pose.PDBInfo( Pose ))
    if Args.renumber_pose:
      rosetta.dump_pdb(Pose, Pdb)
    else:
      rosetta.dump_pdb(Pose, Pdb.replace('.pdb', '_renumbered.pdb'))

    AllConstraints = get_pose_nc_constraints(Pose, SasaScale)

    CstName = Pdb.replace('.pdb', '_nc.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(AllConstraints) 

    # # Non-rosetta part, this finds primary sequence repeats to duplicate 
    # PdbWDAG = pdb_wdag(Pdb, 3, 4.7, 0.25, 15 )
    # PdbWDAG.find_downstream_neighbors()
    # RepeatChains = PdbWDAG.find_repeat_chains()
    
    # ConsolidatedRepeatStarts, TandemRepeats = consolidate_repeats(RepeatChains)

    # rosetta.dump_pose(Pose, '.base.pdb')
    # FirstResidue = 1
    # LastResidue = Pose.n_residue()
    # time.sleep(0.5)

    # # loops through all repeats in start 
    # for RepeatStart in ConsolidatedRepeatStarts:
    #   print
    #   print RepeatStart
    #   print ConsolidatedRepeatStarts[RepeatStart]

    #   RepeatSpacings = TandemRepeats[RepeatStart]
    #   Length = RepeatSpacings[0]

    #   for FirstShiftMultipler in range(len(RepeatSpacings[:-1])):
    #     SecondShiftMultipler = i+1
    #     StartFirst = RepeatStart + (Length * FirstShiftMultipler ) 
    #     EndFirst = StartFirst + Length - 1 
    #     StartSecond = Repeat + (Length * FirstShiftMultipler )
    #     EndSecond = StartSecond + Length - 1
    #     print StartFirst, EndFirst
    #     print StartSecond, EndSecond

    # dump_pdb( (Pose)arg1, (OStream)out, (vector1_Size)residue_indices [, (str)tag='1']) -> None :
    #     for writing a specified subset of residues in pdb format
        
      # ExtendPdbOutput = subprocess.check_output([extend_pdb.py,
      #                                           '-pdb_file', '.base.pdb',
      #                                           '-start_N_cap', FirstResidue,
      #                                           '-end_N_cap', RepeatStart-1,
      #                                           '-start_repeat_1' START_REPEAT_1,
      #                                           '-end_repeat_1', END_REPEAT_1,
      #                                           '-start_repeat_2' START_REPEAT_2,
      #                                           '-end_repeat_2', END_REPEAT_2,
      #                                           '-start_C_cap' START_C_CAP,
      #                                           '-end_C_cap', LastResidue,
      #                                           '-repeat_number', REPEAT_NUMBER])

if __name__ == "__main__":
  sys.exit(main())