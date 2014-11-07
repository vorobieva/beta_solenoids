#!/bin/bash

for i in `seq 2 10` ; do
/work/fabio/git/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease -database ~/git/Rosetta/main/database -s c0_0_1.pdb @flags_run -out:prefix $i\_ -parser:protocol repeat_SD2_$i.xml > log_$i &
done
