#! /usr/bin/env python

from rosetta import *
rosetta.init()
from sys import argv
import numpy as np
import solenoid_tools as st
import math as math


def main():
    # if len(sys.argv) == 3:
    pdb = sys.argv[1]
    dssp = sys.argv[2]
    abego = sys.argv[3]
#    else:
#        print "Number of arguments not valid"
#    print pdb, dssp, abego
    pose = pose_from_pdb(pdb)
    resi_dict = load_as_dict(dssp, abego, pose)
    resi_dict = strand_count(resi_dict, pose)
    print resi_dict
    tot_strands = max(data[5] for data in resi_dict.values())
    tot_turns = max(data[6] for data in resi_dict.values())
    print "Beta-hairpin twist values."
    twist_angles = sheet_twist(pose, resi_dict)
#    print vector_means
    angles_list = np.zeros([tot_strands, tot_strands, 2])
    for resi1 in range(0, pose.total_residue()):
        for resi2 in range(0, pose.total_residue()):
            if twist_angles[resi1][resi2] != 0:
                strand1 = resi_dict[resi1+1][5]
                strand2 = resi_dict[resi2+1][5]
                angles_list[strand1-1][strand2-1][0] += twist_angles[resi1][resi2]
                angles_list[strand1-1][strand2-1][1] += 1
    # print angles_list
    pairs = [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1]]
    for pair in pairs:
        avg_angle = angles_list[pair[0]-1][pair[1]-1][0]/angles_list[pair[0]-1][pair[1]-1][1]
        angle_degrees = 180.0 - (avg_angle * 180.0/math.pi)
        print "The twist angle for the pair %s is %s." %(pair, angle_degrees)
#    for strand in range(0, tot_strands - 1):
#        twist_angle = st.angle(vector_means[strand], vector_means[strand + 1])
#        print "The twist angle between strands %s and %s is %s " %(strand + 1, strand + 2, twist_angle)
    print "Twist of the individual strands:"
    twist_list = strand_twist(resi_dict, pose, tot_strands)
    for strand in range(0, tot_strands):
        print "The average twist of the strand %s is %s" %(strand + 1, twist_list[strand])
    # print "Twist of the turns (L-handed type I' turns (-52 +/- 7) are favoured in beta-sheets):"
    # turn_twist_list = turn_twist(pose, resi_dict)
    # for turn in range(0, tot_turns):
    #     print "The twist of the turn %s is %s" %(turn +1, turn_twist_list[turn])


def load_as_dict(dssp, abego, pose):
    # create a dictionary with keys corresponding to the resi number.
    resi_dict = dict()
#    resi_list = list()
    for line in open(dssp):
        vals = line.strip('\n').split(' ', 5)
        resi_dict[int(vals[0])] = list(vals[1:])
#        resi_list.append(int(vals[0]))

    # append the ABEGO types to the list in the dictionary.
    for line in open(abego):
        for i in range(0, pose.total_residue()):
            # abego[i] = line[i]
            resi_dict[i + 1].append(line[i])
    return resi_dict


def strand_twist(resi_dict, pose, strand_tot):
    # average twist of individual strands, as described by Fujiwara et al. (Proteins 2014; 82:1484, 1493)
    # BigArrayBuild = []
    # build the array with the coordinates
    # for Position in range(pose.n_residue()):
    # Position += 1
    # BigArrayBuild.append(list(pose.residue(Position).xyz("CA")))
    # sliding windows of 4 residues along the strand
    twist_list = list()
    for strand in range(1, strand_tot + 1):
        total_twist= 0.0
        nbre_residues = 0.0
        for i in range(1, pose.n_residue() + 1):
            if resi_dict[i][5] == strand and resi_dict[i][1] == 'E' and resi_dict[i+1][1] == 'E' and resi_dict[i+2][1] == 'E' and resi_dict[i+3][1] == 'E':
                Ca1 = pose.residue(i).xyz("CA")
                Ca2 = pose.residue(i+1).xyz("CA")
                Ca3 = pose.residue(i+2).xyz("CA")
                Ca4 = pose.residue(i+3).xyz("CA")
                L = (Ca1 + Ca2)/2.0
                M = (Ca2 + Ca3)/2.0
                N = (Ca3 + Ca4)/2.0
                P = (L + M)/2.0
                Q = (M + N)/2.0
                twist = rosetta.numeric.dihedral_degrees(Ca2, P, Q, Ca3)
                total_twist = total_twist + twist
                nbre_residues += 1.0
        twist_avg = total_twist/nbre_residues
        twist_list.append(twist_avg)
    return twist_list


def turn_twist(pose, resi_dict):
    # first count the num of residue per turn
    turn_twist_list = list()
#    for turn in range(1, turn_tot + 1):
    for i in range(1, pose.n_residue() + 1):
        if resi_dict[i][1] == 'E' and resi_dict[i+1][1] == 'T' and resi_dict[i+2][1] == 'T' and resi_dict[i+3][1] == 'E':
            Ca1 = pose.residue(i).xyz("CA")
            Ca2 = pose.residue(i+1).xyz("CA")
            Ca3 = pose.residue(i+2).xyz("CA")
            Ca4 = pose.residue(i+3).xyz("CA")
            twist = rosetta.numeric.dihedral_degrees(Ca1, Ca2, Ca3, Ca4)
            turn_twist_list.append(twist)
    return turn_twist_list


def strand_count(resi_dict, pose):
    # numbering the strands to keep trakking of the position in the protein
    strand_num = 0
    turn_num = 0
    turn = True
    strand = False
    for i in range(1, pose.total_residue() + 1):
        # print key
        if resi_dict[i][1] == 'E' and turn is True:
            strand = True
            strand_num += 1
            turn = False
            resi_dict[i].append(strand_num)
            resi_dict[i].append(0)
        elif resi_dict[i][1] == 'E' and turn is False:
            resi_dict[i].append(strand_num)
            resi_dict[i].append(0)
        elif resi_dict[i][1] == 'T' and strand is True:
            turn = True
            turn_num += 1
            strand = False
            resi_dict[i].append(0)
            resi_dict[i].append(turn_num)
        elif resi_dict[i][1] == 'T' and strand is False:
            resi_dict[i].append(0)
            resi_dict[i].append(turn_num)
        else:
            resi_dict[i].append(0)
            resi_dict[i].append(0)
    return resi_dict


def sheet_twist(pose, resi_dict):
    # twist (or coil) angle between two interacting beta-strands, as described by Madam et al. (Protein 2014; 82: 1721-1733)
#    twist_angles = np.empty((81, 81))
    AB_vectors = np.zeros([pose.total_residue(), 3])
    for resi in range(1, pose.total_residue() + 1):
        if resi_dict[resi][1] == 'E' and resi_dict[resi-1][1] == 'E' and resi_dict[resi + 1][1] == 'E':
            Cb0 = pose.residue(resi - 1).xyz("C")
            N1 = pose.residue(resi).xyz("N")
            Cb1 = pose.residue(resi).xyz("C")
            N2 = pose.residue(resi + 1).xyz("N")
            A = (Cb0 + N1)/2.0
            B = (N2 + Cb1)/2.0
            AB = B - A
            for i in range(0, 3):
                AB_vectors[resi-1][i] = AB[i]
        # else:
        #     AB_vectors[resi-1] = [0, 0, 0]
        # print resi, AB_vectors[resi-1]
    twist_angles = np.zeros([pose.total_residue(), pose.total_residue()])
    for entry in range(0, pose.total_residue()):
        contact1 = int(resi_dict[entry + 1][2])
        contact2 = int(resi_dict[entry + 1][3])
        if any(v != 0 for v in AB_vectors[entry]) and contact1 != 0:
            if any(v != 0 for v in AB_vectors[contact1-1]):
                # print entry +1, contact1, AB_vectors[entry], AB_vectors[contact1-1]
                angle = st.angle(AB_vectors[entry], AB_vectors[contact1-1])
                twist_angles[entry][contact1-1] = angle
        if any(v != 0 for v in AB_vectors[entry]) and contact2 != 0:
            if any(v != 0 for v in AB_vectors[contact2-1]) :
                # print entry +1, contact2, AB_vectors[entry], AB_vectors[contact2-1]
                angle = st.angle(AB_vectors[entry], AB_vectors[contact2-1])
                twist_angles[entry][contact2-1] = angle
    return twist_angles

if __name__ == "__main__":
    main()
