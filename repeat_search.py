#! /usr/bin/env python

# all hail Alex Ford
from interface_fragment_matching.utility.analysis import AtomicSasaCalculator

from multiprocessing import Process
from scipy import spatial

import numpy as np
import subprocess
import argparse
import sys
import os

import rosetta
rosetta.init()

ThreeToOne = {'GLY':'G','ALA':'A','VAL':'V','LEU':'L','ILE':'I','MET':'M','PRO':'P','PHE':'F','TRP':'W','SER':'S','THR':'T','ASN':'N','GLN':'Q','TYR':'Y','CYS':'C','CYD':'C','LYS':'K','ARG':'R','HIS':'H','ASP':'D','GLU':'E','STO':'*','UNK':'U'}
ChainAlphabetIndices = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10, 'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16, 'Q':17, 'R':18, 'S':19, 'T':20, 'U':21, 'V':22, 'W':23, 'X':24, 'Y':25, 'Z':26 }

class alpha_carbon:
	""" Calpha node for wDAG searches of proteins for repeats """
	def __init__(self, Number, CoordinateArray):
		self.CoordinateArray = CoordinateArray
		self.Number = int(Number)
		self.DownstreamNeighbors = {} # keyed with residue number, value is displacement vector to that residue's own CA



class pdb_wdag:
	""" Init method takes protein and spawns Calpha instances to populate """
	def __init__( self, Pdb, MinRepeats, DistanceTarget, DistanceFlex, AngleFlex ):
		self.Pdb = Pdb
		self.MinRepeats = MinRepeats
		self.CalphaArray = self.pdb_coord_array(Pdb)
		self.DistanceTarget = DistanceTarget
		self.DistanceFlex = DistanceFlex
		self.AngleFlex = np.radians(AngleFlex) # converts to radians

		# make instance of alpha_carbon for each row
		self.CalphaInstances = [ alpha_carbon( ResRow[0], np.array(ResRow[1:]) ) for ResRow in self.CalphaArray ]


	def find_repeat_chains(self):
		CheckRedundancyChecker = {}
		Repeats = []
		Count = 0
		for Calpha in self.CalphaInstances:
			for NeighborNumber in Calpha.DownstreamNeighbors:
				'''
				Repeat chains start with each of Calphas downstream neighbors 
				and are recursively extended so long as VecA ~= VecB
				
				   VecA      VecB
				CA- - - > CA- - - > CA
				'''
				Count += 1
				RepeatChain = [ Calpha.Number ]
				RecursiveCalpha = Calpha
				DownstreamNumber = NeighborNumber
				
				try:
					CheckRedundancyChecker['%d_%d'%(Calpha.Number, DownstreamNumber)]

				except KeyError:
					ChainExtended = 1

				while ChainExtended:
					RepeatChain.append(DownstreamNumber)
					DownstreamCalpha = self.CalphaDict[DownstreamNumber]					
					VectorA = RecursiveCalpha.DownstreamNeighbors[DownstreamNumber]
					
					ChainExtended = 0 # not extended unless valid extension found
					for NextDownstreamNumber in DownstreamCalpha.DownstreamNeighbors:
						
						# To prevent redundant calculations, upstream downstream pairs already investigated aren't considered					
						try:
							CheckRedundancyChecker['%d_%d'%(DownstreamNumber, NextDownstreamNumber)]
						# only residue pairs not in redundancy checking dictionary are considering
						except KeyError:
							VectorB = DownstreamCalpha.DownstreamNeighbors[NextDownstreamNumber]
						
							if spatial.distance.cosine(VectorA, VectorB) < self.AngleFlex:
								CheckRedundancyChecker['%d_%d'%(DownstreamNumber, NextDownstreamNumber)] = True
								DownstreamNumber = NextDownstreamNumber
								RecursiveCalpha = DownstreamCalpha
								ChainExtended = 1
								break

				if len(RepeatChain) >= self.MinRepeats:
					Repeats.append(RepeatChain)
					# print 'select repeat%d, resi %s'%(Count,'+'.join([str(Number).split('.')[0] for Number in RepeatChain]))

		return Repeats

	def find_downstream_neighbors(self):
		''' To set up DAG Calphas become nodes'''
		self.CalphaDict = {}
		for i, Instance in enumerate(self.CalphaInstances):
			# To make the graph directed, repeat chains will always be arranged from 
			for Subsequent in self.CalphaInstances[i:]:
				# Displacement from lower numbered residue to higher number
				Displacement = Subsequent.CoordinateArray - Instance.CoordinateArray
				# Magnitute of displacement vector
				Distance = (np.sum([Component**2 for Component in Displacement]))**0.5
				# check if displacement of given Instance (upstream) and Subsequent (downstream) has magnitute within target range
				if np.abs(Distance - self.DistanceTarget) <= self.DistanceFlex:
					# if in range record Subsequent as downstream neighbor 
					Instance.DownstreamNeighbors[Subsequent.Number] = Displacement
			# add residue alpha carbon instance to dictionary
			self.CalphaDict[Instance.Number] = Instance


	def pdb_coord_array(self, Pdb):
		if len(Pdb) == 5:
			Chain = Pdb[4]
			Pdb = Pdb[:4]
		else:
			assert len(Pdb) == 4, 'Pdb ID must equal 4 or 5'
			Chain = '.'

		# gets full path to pdb in lab database
		PdbfullPath = ''.join( ['/lab/databases/pdb_clean/', Pdb[1:3].lower(), '/', Pdb, '.pdb'] )

		CalphaLines = subprocess.check_output(['grep', 'ATOM.........CA......%s'%Chain, PdbfullPath]).strip('\n').split('\n')
		CalphaValues = [ [ int(Line[22:26]), float(Line[30:38]), float(Line[38:46]), float(Line[46:54]) ] for Line in CalphaLines]

		self.ResidueIndentityDict = { int(Line[22:26]) : ThreeToOne[Line[17:20]] for Line in CalphaLines}
		return np.array( CalphaValues )



def score_repeats(RepeatChains, AminoAcidIdentityDict, AminoAcidScoreMatrix):
	SortingRepeats = []
	for Repeat in RepeatChains:
		ResidueIDs = [ AminoAcidIdentityDict[Number] for Number in Repeat ]
		Score = np.sum( [ AminoAcidScoreMatrix[Residue] for Residue in ResidueIDs] )
		SortingRepeats.append( [Score, ResidueIDs, Repeat] )
	SortingRepeats.sort()
	SortingRepeats.reverse()
	return SortingRepeats



def check_sasa(Pdb, ResidueSubsets, StartingResidue, SasaProbeRadius):
	''' Uses Alex's AtomicSasaCalculator to calculate average SASA (surface area solvent accessibility) for residue sets input '''
	PdbfullPath = ''.join( ['/lab/databases/pdb_clean/', Pdb[1:3].lower(), '/', Pdb[0:4], '.pdb'] )
	
	# Load pdb into a rosetta pose object 
	PdbPose = rosetta.pose_from_pdb(PdbfullPath)

	if len(Pdb) > 4:
		TargetChainIndex = ChainAlphabetIndices[Pdb[4]]
		PdbChains = PdbPose.split_by_chain()
		# Silly loop to get pose with only the desired chain, there probably is a better way to do this
		for i, Chain in enumerate(PdbChains):
			if i == TargetChainIndex:
				PdbPose = Chain
				break

	# initalize Alex's AtomicSasaCalculator
	SasaCalculator = AtomicSasaCalculator(probe_radius=SasaProbeRadius)
	# get array of residue sasa's
	ResidueSasa = SasaCalculator.calculate_per_residue_sasa(PdbPose)

	SubsetAverageSasas = []
	# RepeatMinimumSasas = []
	# RepeatMaximumSasas = []

	count = 0
	for Residues in ResidueSubsets:
		# Converts residue number from pdb to appropriate index for sasa array
		ResidueIndices = [ResNum - StartingResidue for ResNum in Residues]
		print
		print 'select rep%d, resi'%count, '+'.join([str(R) for R in Residues])
		print [ ResidueSasa[ResIndex] for ResIndex in ResidueIndices ]
		print np.mean( [ ResidueSasa[ResIndex] for ResIndex in ResidueIndices ] )
		
		SubsetAverageSasas.append( np.mean( [ ResidueSasa[ResIndex] for ResIndex in ResidueIndices ] ) )
		# RepeatMinimumSasas.append( min( [ ResidueSasa[ResIndex] for ResIndex in ResidueIndices ] ) )
		# RepeatMaximumSasas.append( max( [ ResidueSasa[ResIndex] for ResIndex in ResidueIndices ] ) )

		count += 1
	
	return SubsetAverageSasas



def main(argv=None):
	if argv is None:
		argv = sys.argv
	
	ArgParser = argparse.ArgumentParser(description=' repeat_search.py arguments ')
	# Required arguments:
	ArgParser.add_argument('-pdbs', type=str, help=' pdb ID list to check ', required=True)
	ArgParser.add_argument('-outfolder', type=str, help=' output directory ', required=True)
	ArgParser.add_argument('-aa_matrix', type=str, help=' weight matrix for amino acid types to rank repeat rows based on desired amino acid composition ', required=True)
	# Optional arguments:
	ArgParser.add_argument('-repeat_class', type=str, default=False, help=' type of repeat a la repeatsdb.bio.unipd, e.g. "III.1" ')
	ArgParser.add_argument('-dist_targ',  type=float, default=5.0,  help=' neighboring repeat target spacing ')
	ArgParser.add_argument('-dist_flex',  type=float, default=0.5, help=' neighboring repeat spacing flexibility ')
	ArgParser.add_argument('-angle_flex', type=float, default=5.0, help=' maximum cosine (in degrees) difference between displacement vectors considered repeatative ')
	ArgParser.add_argument('-min_repeats', type=int, default=5, help=' mimium number of repeats to record ')
	ArgParser.add_argument('-min_sasa',  type=float, default=0.0,  help=' minimum average residue sasa value for reported repeats ')
	ArgParser.add_argument('-max_sasa',  type=float, default=5000.0,  help=' maximum average residue sasa value for reported repeats ')
	ArgParser.add_argument('-sasa_probe_radius', type=float, default=2.2,  help=' probe radius for sasa calculations ')
	args = ArgParser.parse_args()

	# checks for output folder and makes it if not found
	if not os.path.isdir(args.outfolder):
		subprocess.check_output(['mkdir', args.outfolder])
	
	# parses aa matrix
	with open(args.aa_matrix, 'r') as MatrixFile:
		MatrixLines = MatrixFile.readlines()
		assert len(MatrixLines) == 20, '%s contains %d lines, matrix must have 20 entries (all amino acids)'%(args.aa_matrix, len(MatrixLines) )
		MatrixLines = [Line.split('\t') for Line in MatrixLines]
		AminoAcidScoreMatrix = {Line[0]:float(Line[1]) for Line in MatrixLines}

	''' 
	Load pdb list
	'''
	with open(args.pdbs, 'r') as PdbListFile:
		PdbList = PdbListFile.readlines()

	PdbList = [ Pdb.replace('\n', '').split('\t') for Pdb in PdbList ]

	if args.repeat_class:
		PdbList = [ SplitPdb[0] for SplitPdb in PdbList if SplitPdb[2] == args.repeat_class ] 
	else:
		PdbList = [ SplitPdb[0] for SplitPdb in PdbList ]
	
	# output values for log
	print '# Beginning search of %d structures for repeats with:'%len(PdbList)
	print '# target spacing: %f'%args.dist_targ
	print '# spacing flexibility (+-): %f'%args.dist_flex
	print '# angle flexibility: %f'%args.angle_flex
	print '# minimum repeats: %d'%args.min_repeats
	print '# outputing results to: %s'%args.outfolder

	for Pdb in PdbList:
		
		PdbWDAG = pdb_wdag(Pdb, args.min_repeats, args.dist_targ, args.dist_flex, args.angle_flex )
		
		PdbWDAG.find_downstream_neighbors()
		RepeatChains = PdbWDAG.find_repeat_chains()
				
		ScoredRepeats = score_repeats(RepeatChains, PdbWDAG.ResidueIndentityDict, AminoAcidScoreMatrix)
		print '\n'*2
		print Pdb
		ChainAverageSasas = check_sasa(Pdb, RepeatChains, PdbWDAG.CalphaArray[0][0], args.sasa_probe_radius)

		# print StartingResidue
		print ScoredRepeats
		print ChainAverageSasas

		# if __name__ == '__main__':
		#     p = Process(target=check_sasa, args=(Pdb, RepeatChains, MinValue, LogFile))
		#     p.start()
		#     p.join()
	
		# for Calpha in PdbWDAG.CalphaInstances:
		# 	print Calpha.Number, 
		# 	print Calpha.DownstreamNeighbors

	# calc_total_sasa(pose, 1.5)


if __name__ == "__main__":
	sys.exit(main())