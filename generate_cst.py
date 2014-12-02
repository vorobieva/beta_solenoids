#! /usr/bin/env python

InfoString = ''' 
This script is to generate ATOMPAIR constraints for Rosetta,
 by default for all downstreams and oxygens
                    within 3 angstroms 
                    on non-neighbor residues 
                    within an input pose. 

'''

# uncomment just next line and copy block in multiline string for ipython mode
'''


from multiprocessing import Process
from scipy import spatial
from Bio import PDB
import numpy as np
import subprocess
import argparse
import sys
import os
import re

if '-h' not in sys.argv:
  import solenoid_tools
  import rosetta
  rosetta.init(extra_options = "-mute basic -mute core -mute protocols")

# '''

sys.argv.extend(['-pdbs', '1EZG.pdb', '-out', './' ])

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


def get_pose_constraints(Pose, MaxDist, MinPositionSeperation, SasaRadius, SasaScale, UpstreamGrep, DownstreamGrep, NeedHydrogen=True):
    '''  '''
    # AlexsSasaCalculator is from Alex's interface_fragment_matching 
    # thanks Alex!
    #
    # This is used to give buried polar contacts more weight. Thanks Alex Ford!
    try:
      from interface_fragment_matching.utility.analysis import AtomicSasaCalculator
      # make instace of Alex's sasa calculator
      AlexsSasaCalculator = AtomicSasaCalculator(probe_radius=SasaRadius)
      ResidueAtomSasa = AlexsSasaCalculator.calculate_per_atom_sasa(Pose)    
    except ImportError:
      ' Error: SASA weighting of contacts requires interface_fragment_matching from Alex Ford '

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
        AtmCoordList.append( np.array(list(Pose.residue(Res).atom(Atm).xyz())) )
      
      # add residue's lists to respective global lists
      ResAtmCoordLists.extend(AtmCoordList)
      ResAtmRecordLists.extend(AtmRecordList)

    ResidueAtomArray = np.array( ResAtmCoordLists )
    ResidueAtomKDTree = spatial.KDTree( ResidueAtomArray )

    ResidueAtomNeighbors = ResidueAtomKDTree.query_ball_point( ResidueAtomArray, MaxDist )
    # ResidueAtomNearNeighbors = ResidueAtomKDTree.query_ball_point( ResidueAtomArray, 2.0 )
    ResidueAtomHydrogens = ResidueAtomKDTree.query_ball_point( ResidueAtomArray, 1.1 )

    # holds constraints before printing
    AllConstraints = [] 
    # holds sorted cst
    AllBackboneBackboneCst = []
    AllBackboneSidechainCst = []
    AllSidechainSidechainCst = []

    # All contacts are from upstream to downstream residues to avoid double counting
    Upstream = []
    for UpIndex, UpXyzCoords in enumerate(ResAtmCoordLists):
      UpRes, UpAtm = ResAtmRecordLists[UpIndex]

      # # loop through residues storing info on oxygens
      # for UpRes in range( 1, Pose.n_residue() + 1 ):
      #   # loop through atoms
      #   for UpAtm in range( 1, Pose.residue(UpRes).natoms() + 1 ):
      UpName = Pose.residue(UpRes).atom_name(UpAtm).replace(' ', '')

      # skip virtual residues
      if Pose.residue(UpRes).is_virtual(UpAtm):
        continue

      #                                this guy 
      #                                 /
      # checks upstream name           V
      if re.match(UpstreamGrep, UpName ): 
        # print '\n'*2
        # print 'UpRes, UpName', UpRes, UpName

        # get neighbors of upstream residues
        NeighborsOfUpstream = ResidueAtomNeighbors[UpIndex]
        
        # prep for loop
        Downstreams = []

        Constraints = []
        BackboneBackboneCst = []
        BackboneSidechainCst = []
        SidechainSidechainCst = []

        # ArbitrayOrderOfAtomNames = {}
        for DownIndex in NeighborsOfUpstream:
          # name presumes downstream, checks with if imediately below
          DownRes, DownAtm = ResAtmRecordLists[DownIndex]

          # checks that downstream residue is dowstream of upstream and passes min primary sequence spacing
          if DownRes - UpRes >= MinPositionSeperation:
            DownName = Pose.residue(DownRes).atom_name(DownAtm).replace(' ', '')
            
            # skip if same atom
            if UpRes == DownRes:
              if UpName == DownName:
                continue

            # skip virtual residues
            if Pose.residue(DownRes).is_virtual(DownAtm):
              continue

            # checks downstream name
            if re.match( DownstreamGrep, DownName ):
              # print 'DownRes, DownName', DownRes, DownName

              PotentialUpstreamHydrogens = ResidueAtomHydrogens[UpIndex]
              UpstreamHydrogens = []
              # print 'PotentialUpstreamHydrogens', PotentialUpstreamHydrogens
              for UpH_I in PotentialUpstreamHydrogens:
                UpH_Res, UpH_Atm = ResAtmRecordLists[UpH_I]
                UpH_Name  = Pose.residue(UpH_Res).atom_name(UpH_Atm).replace(' ', '')
                # print 'UpH_Name', UpH_Name
                if 'H' in UpH_Name:
                  UpstreamHydrogens.append((UpH_Res, UpH_Atm, UpH_Name))
                # print 'UpstreamHydrogens', UpstreamHydrogens

              PotentialDownstreamHydrogens = ResidueAtomHydrogens[DownIndex]
              DownstreamHydrogens = []
              # print 'PotentialDownstreamHydrogens', PotentialDownstreamHydrogens
              for DownH_I in PotentialDownstreamHydrogens:
                DownH_Res, DownH_Atm = ResAtmRecordLists[DownH_I]
                DownH_Name = Pose.residue(DownH_Res).atom_name(DownH_Atm).replace(' ', '')
                # print 'DownH_Name', DownH_Name
                if 'H' in DownH_Name:
                  DownstreamHydrogens.append((DownH_Res, DownH_Atm, DownH_Name))
                # print 'DownstreamHydrogens', DownstreamHydrogens

              # check their is at least one hydrogen in system before adding constraint
              if len(UpstreamHydrogens) or len(DownstreamHydrogens) or NeedHydrogen == False:

                # these trys / excepts seperate 
                # backbone-backbone from 
                # backbone-sidechain from
                # sidechain-sidechain interactions
                # 
                # in future maybe sort into seperate lists, shouldn't rely on ResidueAtomSasa to know what is in backbone
                try:
                  UpstreamSasa = ResidueAtomSasa[UpRes][UpName]
                  DownstreamSasa = ResidueAtomSasa[DownRes][DownName]
                  AverageSasa = np.mean([UpstreamSasa, DownstreamSasa])        
                  BBBB = 1
                  BBSC = SCSC = 0
                except KeyError:                
                  # These lines handle backbone to sidechain interactions
                  # set weight equal to the most buried 
                  try:
                    UpstreamSasa = ResidueAtomSasa[UpRes][UpName]
                    AverageSasa = SasaScale.FloorSasa
                    BBSC = 1
                    BBBB = SCSC = 0
                  except KeyError:
                    try:
                      DownstreamSasa = ResidueAtomSasa[DownRes][DownName]
                      AverageSasa = SasaScale.FloorSasa 
                      BBSC = 1
                      BBBB = SCSC = 0            
                    
                    # set weight of side chain side chain equal to the most buried             
                    except KeyError:
                      AverageSasa = SasaScale.CeilingSasa 
                      SCSC = 1
                      BBSC = BBBB = 0

                # use instance of sasa_scale to calculate weight based on avg sasa of N and O
                SasaBasedWeight = SasaScale.weigh(AverageSasa)
                # print 
                # print 'AverageSasa', AverageSasa
                # print 'SasaBasedWeight', SasaBasedWeight

                # print 'found downstream neighbor %s'%DownName
                DownXyzCoords = np.array( list(Pose.residue(DownRes).atom(DownAtm).xyz()) )
                # print 'DownRes, DownName', DownRes, DownName
                # print 'DownXyzCoords', DownXyzCoords

                # ## Get neighbors for angles and torsions to use with AtomPairs

                SelectUpNeighbors = []
                # iterates through upstream atom neighbors for references for angle
                for UpNeighborIndex in NeighborsOfUpstream:
                  UpNeighborRes, UpNeighborAtm = ResAtmRecordLists[UpNeighborIndex]
                  UpNeighborName = Pose.residue(UpNeighborRes).atom_name(UpNeighborAtm).replace(' ', '')

                  # keep looking if neighbor is hyrdogen
                  if 'H' in UpNeighborName:
                    continue                

                  # skip virtual residues
                  if Pose.residue(UpNeighborRes).is_virtual(UpNeighborAtm):
                    continue

                  # keep looking if neighbor is self
                  if UpNeighborName == UpName and UpNeighborRes == UpRes:
                    continue
                  # keep looking if neighbor is downstream residue again
                  if UpNeighborName == DownName and UpNeighborRes == DownRes:
                    continue
                  UpNeighborCoords = ResAtmCoordLists[UpNeighborIndex]
                  DistanceToNeighbor = solenoid_tools.vector_magnitude( UpXyzCoords - UpNeighborCoords )
                  SelectUpNeighbors.append( (DistanceToNeighbor, UpNeighborName, UpNeighborRes, UpNeighborCoords) )

                # sort by distance to atom, nearest first
                SelectUpNeighbors.sort()                
                UpNeighbor1Tuple = SelectUpNeighbors[0]
                UpNeighbor2Tuple = SelectUpNeighbors[1]
                # print '\n'*2
                # print 'UpRes, UpName', UpRes, UpName
                # print 'UpstreamHydrogens', UpstreamHydrogens
                # print 'SelectUpNeighbors', SelectUpNeighbors

                 # get neighbors of upstream residues
                NeighborsOfDownstream = ResidueAtomNeighbors[DownIndex]
                SelectDownNeighbors = []
                # iterates through upstream atom neighbors for references for angle
                for DownNeighborIndex in NeighborsOfDownstream:
                  DownNeighborRes, DownNeighborAtm = ResAtmRecordLists[DownNeighborIndex]
                  DownNeighborName = Pose.residue(DownNeighborRes).atom_name(DownNeighborAtm).replace(' ', '')

                  # keep looking if neighbor is hyrdogen
                  if 'H' in DownNeighborName:
                    continue                

                  # skip virtual residues
                  if Pose.residue(DownNeighborRes).is_virtual(DownNeighborAtm):
                    continue

                  # keep looking if neighbor is self
                  if DownNeighborName == DownName and DownNeighborRes == DownRes:
                    continue
                  # keep looking if neighbor is upstream residue
                  if DownNeighborName == UpName and DownNeighborRes == UpRes:
                    continue

                  DownNeighborCoords = ResAtmCoordLists[DownNeighborIndex]
                  DistanceToNeighbor = solenoid_tools.vector_magnitude( DownXyzCoords - DownNeighborCoords )
                  SelectDownNeighbors.append( (DistanceToNeighbor, DownNeighborName, DownNeighborRes, DownNeighborCoords) )

                # sort by distance to atom, nearest first
                SelectDownNeighbors.sort()
                DownNeighbor1Tuple = SelectDownNeighbors[0]
                DownNeighbor2Tuple = SelectDownNeighbors[1]
                # print 'DownRes, DownName', DownRes, DownName
                # print 'DownstreamHydrogens', DownstreamHydrogens
                # print 'SelectDownNeighbors', SelectDownNeighbors

                Distance = solenoid_tools.vector_magnitude(DownXyzCoords - UpXyzCoords)
                
                DistanceCst = 'AtomPair %s %d %s %d SCALARWEIGHTEDFUNC %f HARMONIC %.2f 1.0' %( UpName, UpRes, DownName, DownRes, SasaBasedWeight, Distance )

                # Use Biopython for angle and dihedral calculations
                # here 'Vec' means PDB.Vector of atom's xyz coord
                UpstreamVec = PDB.Vector(UpXyzCoords)
                DownstreamVec = PDB.Vector(DownXyzCoords)
                
                UpNeighbor1Vec = PDB.Vector(UpNeighbor1Tuple[3])
                UpNeighbor2Vec = PDB.Vector(UpNeighbor2Tuple[3])
                DownNeighbor1Vec = PDB.Vector(DownNeighbor1Tuple[3])
                DownNeighbor2Vec = PDB.Vector(DownNeighbor2Tuple[3])

                Angle1 = PDB.calc_angle(UpNeighbor1Vec, UpstreamVec, DownstreamVec)
                AngleCst1 = 'Angle %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpNeighbor1Tuple[1], UpNeighbor1Tuple[2], UpName, UpRes, DownName, DownRes, SasaBasedWeight, Angle1 )
                Angle2 = PDB.calc_angle(UpstreamVec, DownstreamVec, DownNeighbor1Vec)
                AngleCst2 = 'Angle %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpName, UpRes, DownName, DownRes, DownNeighbor1Tuple[1], DownNeighbor1Tuple[2], SasaBasedWeight, Angle2 )

                Torsion1 = PDB.calc_dihedral(UpNeighbor2Vec, UpNeighbor1Vec, UpstreamVec, DownstreamVec)
                TorsionCst1 = 'Dihedral %s %d %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpNeighbor2Tuple[1], UpNeighbor2Tuple[2], UpNeighbor1Tuple[1], UpNeighbor1Tuple[2], UpName, UpRes, DownName, DownRes, SasaBasedWeight, Torsion1 )
                Torsion2 = PDB.calc_dihedral(UpNeighbor1Vec, UpstreamVec, DownstreamVec, DownNeighbor1Vec)
                TorsionCst2 = 'Dihedral %s %d %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpNeighbor1Tuple[1], UpNeighbor1Tuple[2], UpName, UpRes, DownName, DownRes, DownNeighbor1Tuple[1], DownNeighbor1Tuple[2], SasaBasedWeight, Torsion2 )
                Torsion3 = PDB.calc_dihedral(UpstreamVec, DownstreamVec, DownNeighbor1Vec, DownNeighbor2Vec)
                TorsionCst3 = 'Dihedral %s %d %s %d %s %d %s %d SCALARWEIGHTEDFUNC %f CIRCULARHARMONIC %.2f 0.5' %( UpName, UpRes, DownName, DownRes, DownNeighbor1Tuple[1], DownNeighbor1Tuple[2], DownNeighbor2Tuple[1], DownNeighbor2Tuple[2], SasaBasedWeight, Torsion3 )

                # adds constraint to running lists of constraints
                Constraints.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )
                if BBBB: BackboneBackboneCst.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )
                if BBSC: BackboneSidechainCst.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )
                if SCSC: SidechainSidechainCst.extend( [DistanceCst, AngleCst1, AngleCst2, TorsionCst1, TorsionCst2, TorsionCst3] )

              # else:
              #   print 'No hydrogen!'
              #   sys.exit()

        AllConstraints.extend(Constraints)
        AllBackboneBackboneCst.extend(BackboneBackboneCst)
        AllBackboneSidechainCst.extend(BackboneSidechainCst)
        AllSidechainSidechainCst.extend(SidechainSidechainCst)

    SortedConstraints = (AllBackboneBackboneCst, AllBackboneSidechainCst, AllSidechainSidechainCst)

    return AllConstraints, SortedConstraints


def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=' nc_cst_gen.py arguments ( -help ) %s'%InfoString)
  # Required arguments:
  ArgParser.add_argument('-pdbs', type=str, nargs='+', help=' input pdbs ', required=True)
  # Optional arguments:
  ArgParser.add_argument('-out', type=str, help=' output directory ', default='./')
  ArgParser.add_argument('-max_dist', type=float, default=3.4, help=' distance between the oxygens and downstreams ')
  ArgParser.add_argument('-min_seq_sep', type=int, default=3, help=' minimum seperation in primary sequece ')
  ArgParser.add_argument('-upstream_atom', type=str, default='[ON]\w?\d?', help=' grep for upstream atoms ')
  ArgParser.add_argument('-downstream_atom', type=str, default='[ON]\w?\d?', help=' grep for downstream atoms ')
  ArgParser.add_argument('-num_repeats', type=int, default=5, help=' number of repeats to extrapolate contacts for ')
  ArgParser.add_argument('-min_sasa',  type=float, default=0.0,  help=' floor for weighting downstream oxygen contacts ')
  ArgParser.add_argument('-min_sasa_weight',  type=float, default=1.0,  help=' weight of floor for downstream oxygen contacts ')
  ArgParser.add_argument('-max_sasa',  type=float, default=5.0,  help=' ceiling for cst weighting downstream oxygen contacts ')
  ArgParser.add_argument('-max_sasa_weight',  type=float, default=0.1,  help=' weight of ceiling for downstream oxygen contacts ')
  ArgParser.add_argument('-sasa_probe_radius', type=float, default=0.8,  help=' probe radius for sasa calculations ')
  ArgParser.add_argument('-renumber_pose', type=bool, default=True, help='True|False renumber pdb residues ' )
  
  ArgParser.add_argument('-disulfide', type=bool, default=True, help='True|False include disulfide constraints ' )  

  Args = ArgParser.parse_args()
  
  # if len(Args.pdbs[0]) == 1:
  #   Args.pdbs = [''.join(Args.pdbs)]

  if Args.out [-1] != '/':
    Args.out = Args.out + '/'

  import rosetta
  rosetta.init(extra_options = "-mute basic -mute core -mute protocols")

  ReportedRepeatCount = 0
  TotalPdbs = len(Args.pdbs)

  # Instance of class to convert sasas to cst weight
  SasaScale = sasa_scale( Args.min_sasa, Args.min_sasa_weight, Args.max_sasa, Args.max_sasa_weight )
  
  for iPdb, Pdb in enumerate(Args.pdbs):
    print ' Working with %s; %d of %d total pdbs '%(Pdb, iPdb+1, TotalPdbs)
    # Starting rosetta  
    Pose = rosetta.pose_from_pdb(Pdb)
    OutputPdb = Args.out+Pdb

    # Sets pdb info so residues in dumped pdbs are same as index 
    Pose.pdb_info(rosetta.core.pose.PDBInfo( Pose ))
    if Args.renumber_pose:
      rosetta.dump_pdb(Pose, OutputPdb)
    else:
      rosetta.dump_pdb(Pose, OutputPdb.replace('.pdb', '_renumbered.pdb'))

    AllConstraints, SortedConstraints = get_pose_constraints(Pose, Args.max_dist, Args.min_seq_sep, Args.sasa_probe_radius, SasaScale, Args.upstream_atom, Args.downstream_atom, True)
    
    if Args.disulfide:
      DisulfAllConstraints, DisulfSortedConstraints = get_pose_constraints(Pose, 2.3, 2, Args.sasa_probe_radius, SasaScale, 'SG', 'SG', False)
      AllConstraints.extend(DisulfAllConstraints)

    print AllConstraints, SortedConstraints
    print 
    print DisulfAllConstraints, DisulfSortedConstraints
    sys.exit()

    CstName = OutputPdb.replace('.pdb', '_All.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(AllConstraints) 

    BackboneBackboneCst, BackboneSidechainCst, SidechainSidechainCst = SortedConstraints

    CstName = OutputPdb.replace('.pdb', '_BBBB.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(BackboneBackboneCst) 
    CstName = OutputPdb.replace('.pdb', '_BBSC.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(BackboneSidechainCst) 
    CstName = OutputPdb.replace('.pdb', '_SCSC.cst')
    with open(CstName, 'w') as CstFile:
      print>>CstFile, '\n'.join(SidechainSidechainCst) 


if __name__ == "__main__":
  sys.exit(main())



# def consolidate_repeats(ListOfRepeatPositions):
#   ''' takes lists of repeat positions (based on xyz coords) and returns subsets that split the pose longest equal residue number repeat chains '''  
#   TandemIndenticalSpacings = {} #{0:[]}
  
#   for RepeatPositions in ListOfRepeatPositions:
#     RepeatPositions = RepeatPositions[:]
#     RepeatSpacings = [ RepeatPositions[i+1] - RepeatPositions[i] for i in range(len(RepeatPositions)-1) ]

#     # appending buffer value for added last repeat chain to TandemIndenticalSpacings
#     RepeatSpacings.append(0)
#     RepeatPositions.append(0)

#     LastSpacing = 0
#     RepeatChain = [RepeatSpacings[0]]
#     Start = RepeatPositions[0]
    
#     for i, Spacing in enumerate(RepeatSpacings):
#       if Spacing == LastSpacing:
#         RepeatChain.append(Spacing)
#       else:
#         TandemIndenticalSpacings[Start] = RepeatChain
#         Start = RepeatPositions[i]
#         RepeatChain = [RepeatSpacings[i]]

#       LastSpacing = Spacing
  
#   MaxNumberStart = 0
#   EqualLengthStarts = []

#   for Start in TandemIndenticalSpacings:
#     if len(TandemIndenticalSpacings[Start]) > len(TandemIndenticalSpacings[MaxNumberStart]):
#       MaxNumberStart = Start
#       EqualLengthStarts = []
#     elif len(TandemIndenticalSpacings[Start]) == len(TandemIndenticalSpacings[MaxNumberStart]):
#       EqualLengthStarts.append(Start)

#   for RepeatStart in EqualLengthStarts:
#     try:
#       assert TandemIndenticalSpacings[MaxNumberStart][0] == TandemIndenticalSpacings[RepeatStart][0], ' different repeat spacings have same max copy number ' 
#     # This is explicted (semi) silenced to prevent a large job from stopping at some later date
#     except AssertionError:
#       print '\n\n\n LOOKOUT ERROR: multiple different repeat spacings have max copy number. \n\n\n'

#   MaxNumberRepeatStarts = [MaxNumberStart] + EqualLengthStarts
#   return MaxNumberRepeatStarts, TandemIndenticalSpacings
