
beta_solenoids
==============

to use easily :

export PYTHONPATH=/work/pylesh/python_local.bak:$PYTHONPATH
source /work/pylesh/virtualenvs/dev/bin/activate

# Main pipeline:
Short verison:
generate_backbones.py -pdbs 3fy3.pdb -repeat 4 -max_turns_per_repeat 2

generate_cst.py -pdbs 3fy3.pdb

expand_cst.py -ref_pdb 3fy3.pdb -ref_cst 3fy3_All.cst -repeat_pdb_tag _3fy3

!!! Be careful and considerate with this script as runs things in parrallel !!!
optimize_repeat_structures.py -pdb_stem _3fy3 -thread 10
OR
!!! Now your also generating silent processes !!!
nohup optimize_repeat_structures.py -pdb_stem _3fy3 -thread 10 > log.txt &

#generate_backbones.py 
generate_backbones.py -pdbs 3fy3.pdb -repeat 4 -max_turns_per_repeat 2
 Give pdb(s) with beta solenoid repeat. Uses Pyrosetta directly and via Daniel Silva's RMSD aligner function in solenoid_tools
 
#generate_cst.py -pdbs RepeatPose.pdb
Makes cst file for each input pdb. For options see -h. Uses Pyrosetta directly and via SASA wrapper from Alex Ford's interface_fragment_matching for solvation calculation
By default, make contraints for all Nitrogen to Oxygen ( and O to O contacts with at least one hydrogen in the system ) and disulfides. 
 
 #expand_constraints.py
 expand_cst.py -ref_pdb 3fy3.pdb -ref_cst 3fy3_All.cst -repeat_pdb_tag _3fy3

 optional arguments:
  -h, --help            show this help message and exit
  -ref_pdb REF_PDB      reference pdb
  -ref_cst REF_CST      corresponding to reference pdb
  -repeat_pdb_tag       input pdb tag
  -out OUT              Output directory
 
 
# optimize_repeat_structures.py
optional arguments:
  -h, --help          show this help message and exit
  -pdb_stem PDB_STEM  pdb stem, start of globs for pdbs and csts
  -thread THREAD      number of parallel threads to run 
  
usage: score_and_select_2d.py [-h] -pdb_glob PDB_GLOB -native NATIVE -out
                           -param []
                           OUT [-score SCORE] [-plot PLOT] [-norm NORM]
                           [-name NAME]

for plotting pdb scores and selecting subsets based on absolute or per residue
scores


# score_and_select_2d.py

For -plot function (on by default, 0 for off) you need a plotly account from https://plot.ly/python/. Github accounts work. You must provide 

optional arguments:
  -h, --help       show this help message and exit
  -pdb_glob        pdb stem, start of globs for pdbs and csts
  -native_pdb      pdb to compare designs against
  -out OUT         folder to move files to
  -score SCORE     select all structures with less than this REU /
                   residue
  -plot PLOT       0|(1) plot scores with plotly
  -norm NORM       0|(1) normalize scores by residue
  -name NAME       plot tag


# Misc other scripts

# cst_to_pymol.py
cst_to_pymol.py -cst 3fy3_All.cst
use to visualize cst file in pymol (copy paste output into pymol session with pdb corresponding to cst file)

# repeat_search.py
Use to find repeat motifs involving particular residues, i.e. aromatics, charges, etc. Bias toward favored residues in a amino acid matirx like in matrix.tsv.

script needs to be in directory with repeatdb.txt (included in repo) which contains data from repeatsdb.bio.unipd


optional arguments:
  -h, --help            show this help message and exit
  -pdbs PDBS            pdb ID list to check
  -output OUTPUT        output directory
  -aa_matrix            weight matrix for amino acid types to rank repeat rows
                        based on desired amino acid composition
  -repeat_class         type of repeat a la repeatsdb.bio.unipd, e.g. "III.1"
  -dist_targ            neighboring repeat target spacing
  -dist_flex            neighboring repeat spacing flexibility
  -angle_flex           maximum cosine (in degrees) difference between
                        displacement vectors considered repeatative
  -min_repeats          mimium number of repeats to record
  -min_aa_score         minimum AA composition sum based on scores in matrix
  -min_sasa             minimum average residue sasa value for reported
                        repeats
  -max_sasa             maximum average residue sasa value for reported
                        repeats
  -sasa_probe_radius    probe radius for sasa calculations


