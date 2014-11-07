#!/usr/bin/env python

#Fabio Parmeggiani, University of Washington, Baker Lab, 2013

#script takes pdb and range of residues for first and second repeat to prepare file for calculation of helical parameters
#arg1: pdb file
#arg2: start repeat 1
#arg3: end repeat 1
#arg4: start repeat 2
#arg5: end repeat 2

#arg6: number of repeats

# the two repeats should have the same number of residues and be adjacent

#UPDATE#
#produce also capping repeats to be attached at the end
#pdbs for caps include the caps and the adjacent interbnal repeat to allow superposition afterwards 
#arg7: start cap 1
#arg8: end cap 1
#arg9: start cap 2
#arg10: end cap 2

import sys
import Bio.PDB
import numpy
import re

import argparse
import subprocess

# if len(sys.argv) <=5:
# #    print '\n'
#     print 'arguments: <pdb_file> <start_repeat_1> <end_repeat_1> <start_repeat_2> <end_repeat_2> <number of final repeats>'   
#     print 'note: the two repeats should have the same number of residues and be adjacent'
#     print 'if you want to consider also the caps, add further arguments'
#     print 'additional arguments: <start_cap_N> <end_cap_N> <start_cap_C> <end_cap_C>'
#     print 'careful that the begin and end of caps matches with internal repeats, without superpositions or gaps'
# #    print '\n '
    # sys.exit()

ArgParser = argparse.ArgumentParser(description=' extend_pdb.py arguments ')
# Required arguments:
ArgParser.add_argument('-pdb_file', type=str, help=' starting pdb ', required=True)

ArgParser.add_argument('-start_N_cap', type=int, help=' first residue N terminal cap ', required=True)
ArgParser.add_argument('-end_N_cap', type=int, help=' first residue N terminal cap ', required=True)

ArgParser.add_argument('-start_repeat_1', type=int, help=' first residue in repeat 1 ', required=True)
ArgParser.add_argument('-end_repeat_1', type=int, help=' last residue in repeat 1 ', required=True)
ArgParser.add_argument('-start_repeat_2', type=int, help=' first residue in repeat 2 ', required=True)
ArgParser.add_argument('-end_repeat_2', type=int, help=' last residue in repeat 2 ', required=True)

ArgParser.add_argument('-start_C_cap', type=int, help=' first residue C terminal cap ', required=True)
ArgParser.add_argument('-end_C_cap', type=int, help=' first residue C terminal cap ', required=True)

ArgParser.add_argument('-repeat_number', type=int, help=' number of repeats ', required=True)

ArgParser.add_argument('-relax', type=int, help=' relax boolean, 1 == call relax; 0 == do not relax ', default=0)

args = ArgParser.parse_args()

rep = args.repeat_number
lenRep = args.end_repeat_1- args.start_repeat_1+1

#Extract adjacent repeats as A and B chains in output pdb
pdb=[]
out_file = open("prep_%s" %( args.pdb_file ), 'w')
for i in open( args.pdb_file ):
   if i.split()[0]=='ATOM' and int(i.split()[5])>= args.start_repeat_1 and int(i.split()[5])<= args.end_repeat_1:
      pdb.append(i)
len_part=len(pdb)
#print pdb
for i in range(len_part):
   out_file.write(pdb[i][:21]+'A'+pdb[i][22:])
for i in open( args.pdb_file ):
   if i.split()[0]=='ATOM' and int(i.split()[5])>= args.start_repeat_2 and int(i.split()[5])<= args.end_repeat_2 :
      pdb.append(i)
#print pdb
for i in range(len_part, len(pdb)):
   out_file.write(pdb[i][:21]+'B'+pdb[i][22:])
out_file.close()

    
#extract caps
caps=1
# if len(sys.argv) == 11:
#Ncap
pdb=[]
out_file = open("capN_%s" %( args.pdb_file ), 'w')
for i in open( args.pdb_file ):
   if i.split()[0]=='ATOM' and int(i.split()[5])>= args.start_N_cap and int(i.split()[5])<= args.end_N_cap:
      pdb.append(i)
len_part=len(pdb)
#print pdb
for i in range(len_part):
   out_file.write(pdb[i][:21]+'A'+pdb[i][22:])
for i in open( args.pdb_file ):
   if i.split()[0]=='ATOM' and int(i.split()[5])> args.end_N_cap and int(i.split()[5])<= args.end_N_cap+lenRep:
      pdb.append(i)
for i in range(len_part, len(pdb)):
   out_file.write(pdb[i][:21]+'B'+pdb[i][22:])
out_file.close()


#Ccap   
pdb=[]
out_file = open("capC_%s" %( args.pdb_file ), 'w')

for i in open( args.pdb_file ):
   if i.split()[0]=='ATOM' and int(i.split()[5])>= ( args.start_C_cap - lenRep ) and int(i.split()[5])< args.start_C_cap:
      pdb.append(i)
#print pdb   
len_part=len(pdb)
for i in range(len_part):
   out_file.write(pdb[i][:21]+'B'+pdb[i][22:])
for i in open( args.pdb_file ):
   if i.split()[0]=='ATOM' and int(i.split()[5])>= args.start_C_cap and int(i.split()[5])<= args.end_C_cap :
      pdb.append(i)
for i in range(len_part, len(pdb)):
   out_file.write(pdb[i][:21]+'A'+pdb[i][22:])
out_file.close()

#superposition
#adapted from align.py http://combichem.blogspot.com/2013/08/aligning-pdb-structures-with-biopython.html

# Start the parser
pdb_parser = Bio.PDB.PDBParser(QUIET = True)

# Get the structures
ref_structure = pdb_parser.get_structure("reference", "prep_%s" %( args.pdb_file ))
sample_structure = pdb_parser.get_structure("sample", "prep_%s" %( args.pdb_file ))

# Use the first model in the pdb-files for alignment
# Change the number 0 if you want to align to another structure
ref_model = ref_structure[0]
sample_model = sample_structure[0]

#caps present 
shiftN=0
if caps==1:
   # Get the cap structures
   sample_structureN = pdb_parser.get_structure("sample", "capN_%s" %( args.pdb_file ))
   sample_structureC = pdb_parser.get_structure("sample", "capC_%s" %( args.pdb_file ))
   sample_capN = sample_structureN[0]
   sample_capC = sample_structureC[0]

   #renumber caps
   #capN
   seqN=range(1,  args.end_N_cap- args.start_N_cap+2)   
   t = 0
   for residue in sample_capN['A']:
#      print residue
      residue.id = (' ', seqN[t], ' ')
      t += 1

   seqN=range( args.end_N_cap,  args.end_N_cap+lenRep+1)   
   t = 0
   for residue in sample_capN['B']:
#      print residue
      residue.id = (' ', seqN[t], ' ')
      t += 1
   
   shiftN=len(sample_capN['A'])
   print shiftN

   #capC
   seqC=range( args.end_N_cap- args.start_N_cap+1+lenRep*rep+1,   args.end_N_cap- args.start_N_cap+1+lenRep*rep+ args.end_C_cap - args.start_C_cap+2)   
#   print seqC
   t = 0
   for residue in sample_capC['A']:
      residue.id = (' ', seqC[t], ' ')
      t += 1



models=[]
for i in range(rep) :

   # Make a list of the atoms (in the structures) you wish to align.
   # In this case we use CA atoms whose index is in the specified range
   ref_atoms = []
   sample_atoms = []
    
   # Iterate of all chains in the model in order to find all residues to align
   for ref_res in ref_model['B']:
      ref_atoms.append(ref_res['CA'])
   #print len(ref_atoms)
   #print ref_atoms
    
   # Do the same for the sample structure
   for sample_res in sample_model['A']:
      sample_atoms.append(sample_res['CA'])
   #print len(sample_atoms)
   #print sample_atoms
    
   # Now we initiate the superimposer:
   super_imposer = Bio.PDB.Superimposer()
   super_imposer.set_atoms(ref_atoms, sample_atoms)
   super_imposer.apply(sample_model.get_atoms())
    
   # Print RMSD:
   print 'rmsd_' + str(i+1) + ' ' + str(super_imposer.rms)
   
   #reset reference for next iteration of chain growth by superposition
   ref_model = sample_model

   #renumber
   seq=range((i*lenRep)+1+shiftN, (lenRep*(i+1))+1+shiftN)   
#   print seq
   t = 0
   for residue in sample_model['A']:
      residue.id = (' ', seq[t], ' ')
      t += 1

   #collect and dump as pdb files all the aligned version of chain A
   models.append(sample_model['A'])

   # Save the aligned version of repeats
   io = Bio.PDB.PDBIO()
   io.set_structure(models[i])
   io.save("rep_%s.pdb" %(i))

# kep only first repeat after alignment




# align capN
ref_structure = pdb_parser.get_structure("reference", "rep_0.pdb")
ref_model = ref_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
ref_atoms = []
sample_atoms = []
 
# Iterate of all chains in the model in order to find all residues to align
for ref_res in ref_model['A']:
   ref_atoms.append(ref_res['CA'])
#print len(ref_atoms)
#print ref_atoms
 
# Do the same for the sample structure


### caused bug
for sample_res in sample_capN['B']:
   sample_atoms.append(sample_res['CA'])



#print len(sample_atoms)
#print sample_atoms
 
# Now we initiate the superimposer:
super_imposer = Bio.PDB.Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_capN.get_atoms())
 
# Print RMSD:
print ' '
print 'rmsd_' + 'capN' + ' ' + str(super_imposer.rms)

# Save the aligned version of repeats
io = Bio.PDB.PDBIO()
io.set_structure(sample_capN['A'])
io.save("repCapN_%s.pdb" %( args.pdb_file ))


# align capC
ref_structure = pdb_parser.get_structure("reference", "rep_%s.pdb" %(rep-1))
ref_model = ref_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
ref_atoms = []
sample_atoms = []
 
# Iterate of all chains in the model in order to find all residues to align
for ref_res in ref_model['A']:
   ref_atoms.append(ref_res['CA'])
#print len(ref_atoms)
#print ref_atoms
 
# Do the same for the sample structure
for sample_res in sample_capC['B']:
   sample_atoms.append(sample_res['CA'])
#print len(sample_atoms)
#print sample_atoms
 
# Now we initiate the superimposer:
super_imposer = Bio.PDB.Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_capC.get_atoms())
 
# Print RMSD:
print ' '
print 'rmsd_' + 'capC' + ' ' + str(super_imposer.rms)

# Save the aligned version of repeats
io = Bio.PDB.PDBIO()
io.set_structure(sample_capC['A'])
io.save("repCapC_%s.pdb" %( args.pdb_file ))


full=[]
if caps==1:
   for j in open("repCapN_%s.pdb" %( args.pdb_file )):
      #print j
      if j.split()[0]=='ATOM':
         full.append(j)
for i in range(rep):
   for j in open("rep_%s.pdb" %(i)):
#      print j
      if j.split()[0]=='ATOM':
         full.append(j)
if caps==1:
   for j in open("repCapC_%s.pdb" %( args.pdb_file )):
      #print j
      if j.split()[0]=='ATOM':
         full.append(j)
full.append('END')

output_file_name = 'rep%s_%s'%(rep,  args.pdb_file )

#write full file
with open(output_file_name, 'w') as out_file:
   for i in range(len(full)):
      out_file.write(full[i])

#
if args.relax == 1:

   # relax first:
   subprocess.check_output( [ 'relax.default.linuxgccrelease', '-ignore_unrecognized_res', '-relax:constrain_relax_to_start_coords', '-ex1', '-ex2', '-use_input_sc', '-s', output_file_name] )
   relax_output_name = re.sub(r'(.*).pdb', r'\1_relax.pdb', output_file_name)
   subprocess.check_output( [ 'mv', '%s_0001.pdb'%re.sub(r'(.*).pdb', r'\1', output_file_name), relax_output_name ] )

   # then idealize:
   subprocess.check_output( [ 'idealize_jd2', '-s', relax_output_name ] )
   relax_ideal_output_name = re.sub(r'(.*).pdb', r'\1_ideal.pdb', relax_output_name)
   subprocess.check_output( [ 'mv', '%s_0001.pdb'%re.sub(r'(.*).pdb', r'\1', relax_output_name), relax_ideal_output_name ] )

   # the relax again:
   subprocess.check_output( [ 'relax.default.linuxgccrelease', '-ignore_unrecognized_res', '-relax:constrain_relax_to_start_coords', '-ex1', '-ex2', '-use_input_sc', '-s', relax_ideal_output_name] )
   relax_ideal_relax_output_name = re.sub(r'(.*).pdb', r'\1_relax.pdb', relax_ideal_output_name)
   subprocess.check_output( [ 'mv', '%s_0001.pdb'%re.sub(r'(.*).pdb', r'\1', relax_ideal_output_name), relax_ideal_relax_output_name ] )
   
   # idealize first:
   # subprocess.check_output( [ 'idealize_jd2', '-s', output_file_name ] )
   # idealized_output_name = re.sub(r'(.*).pdb', r'\1_ideal.pdb', output_file_name)
   # subprocess.check_output( [ 'mv', '%s_0001.pdb'%re.sub(r'(.*).pdb', r'\1', output_file_name), idealized_output_name ] )

   # subprocess.check_output( [ 'relax.default.linuxgccrelease', '-ignore_unrecognized_res', '-relax:constrain_relax_to_start_coords', '-ex1', '-ex2', '-use_input_sc', '-s', idealized_output_name] )
   # ideal_relax_output_name = re.sub(r'(.*).pdb', r'\1_relax.pdb', idealized_output_name)
   # subprocess.check_output( [ 'mv', '%s_0001.pdb'%re.sub(r'(.*).pdb', r'\1', idealized_output_name), ideal_relax_output_name ] )
   

