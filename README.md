beta_solenoids
==============

generate_cst.py -pdbs RepeatPose.pdb
# makes one cst file for each input pdb. For (some) options see -h 
# uses Pyrosetta directly and via SASA wrapper from Alex Ford's interface_fragment_matching for solvation calculation

generate_backbones.py -pdbs RepeatPose1.pdb RepeatPose2.pdb ... 
# given repeat pdbs 
# uses Pyrosetta directly and via Daniel Silva's RMSD aligner

