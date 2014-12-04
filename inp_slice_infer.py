#! /usr/bin/env python

# >gi|514085005|gb|EPF65740.1| Ice nucleation protein [Pseudomonas syringae pv. syringae SM]
# MNLDKALVLRTCANNMADHCGLIWPASGTVESKYWQSTRRHENGLVGLLWGAGTSAFLSVHADARWIVCE
# VAVADIISLEEPGMVKFPRAEVVHVGDRISASHFISARQADPASTPTPTPMTAATPPPTPATANVTLPVA
# EQASHEVFDVALVSAAAPPVNTLPVTTPQNLQTA
INP_SLICES = '''TYGSTLSGDNHSRLIA
GYGSNETAGNHSDLIA
GYGSTGTAGSDSSLVA
GYGSTQTAGGDSALTA
GYGSTQTAREGSNLTA
GYGSTGTAGSDSSLIA
GYGSTQTSGEDSSLTA
GYGSTQTAQEGSNLTA
GYGSTGTAGSDSSLIA
GYGSTQTSGGDSSLTA
GYGSTQTAQEGSNLTA
GYGSTGTAGSDSSLIA
GYGSTQTSGEDSSLTA
GYGSTQTAQEGSNLTA
GYGSTGTAGSDSSLIA
GYGSTQTSGGDSSLTA
GYGSTQTAQEGSNLTS
GYGSTGTAGADSSLIA
GYGSTQTSGSDSALTA
GYGSTQTAQEGSNLTA
GYGSTGTAGSDSSLIA
GYGSTQTSGSDSSLTA
GYGSTQTAQEGSNLTA
GYGSTSTAGVDSSLIA
GYGSTQTSGSDSALTA
GYGSTQTAQEGSNLTA
GYGSTGTAGADSSLIA
GYGSTQTSGSDSALTA
GYGSTQTAQEGSNLTA
GYGSTGTAGADSSLIA
GYGSTQTSGSESSLTA
GYGSTQTAREGSTLTA
GYGSTGTAGADSSLIA
GYGSTQTSGSESSLTA
GYGSTQTAQQGSVLTS
GYGSTQTAGAASNLTT
GYGSTGTAGHESFIIA
GYGSTQTAGHKSILTA
GYGSTQTARDGSDLIA
GYGSTGTAGSGSSLIA
GYGSTQTASYRSMLTA
GYGSTQTAREHSDLVT
GYGSTSTAGSNSSLIA
GYGSTQTAGFKSILTA
GYGSTQTAQERSDLVA
GYGSTSTAGYSSSLIA
GYGSTQTAGYGSTLTT
GYGSTQTAQENSSLTT
GYGSTSTAGYSSSLIA
GYGSTQTAGYESTLTA
GYGSTQTAQERSDLVT
GYGSTSTAGYASSLIA
GYGSTQTAGYESTLTA
GYGSTQTAQENSSLTT
GYGSTSTAGFASSLIA
GYGSTQTAGYKSTLTA
GYGSTQTAEYGSSLTA
GYGSTATAGQDSSLIA
GYGSSLTSGIRSFLTA
GYGSTLIAGLRSVLIA
GYGSSLTSGIRSTLTA
GYGSNQIASYGSSLIA
GHESIQVAGNKSMLIA
GKGSSQTAGFRSTLIA
GAGSVQLAGDRSRLIA
GADSNQTAGDRSKLLA
GNNSYLTAGDRSKLTG
GHDCTLMAGDQSRLTA'''.split('\n')
# GKNSILTAGARSKLIGSEGSTLSAGEDSTLIFRLWDGKRYRQLVARTGENGVEADIPYYVNEDDDIVD
# KPDEDDDWIEVE

'''# libraries
from multiprocessing import Pool
import numpy as np
import subprocess
import argparse
import glob
import time
import sys
import os
import re

if '-h' not in sys.argv:
  import rosetta
  # rosetta.init()
  rosetta.init(extra_options = " -ex1 -ex2 -no_optH false -use_input_sc ") # -mute basic -mute core
  from rosetta.protocols.simple_moves import symmetry
  from rosetta.utility import ostream
  
  # from expand_constraints import set_all_weights_zero

'''

sys.argv = [ sys.argv[0], '-thread', '19', '-repeat', '4', '-contact', 'N-O', '-shift', '0', '-exclude', 'GSDN']

def main(argv=None):
  if argv is None:
    argv = sys.argv
  
  ArgParser = argparse.ArgumentParser(description=" args for ice nucleation slice inference ")
  ArgParser.add_argument('-thread', type=int, help=" number of threads to run simultaneously ", default=5 )     
  ArgParser.add_argument('-repeat', type=int, help=" number of repeats to produce ", default=4 )    
  ArgParser.add_argument('-contact', type=str, help="  ", default='N-O' )    
  ArgParser.add_argument('-shift', type=int, help="  ", default=0 )    
  ArgParser.add_argument('-exclude', type=str, help="  ", default='GS' )    
  Args = ArgParser.parse_args()
  
  # print CstHash
  for ThreadChunkNumber in range( (len(Pdbs)/Args.thread) + 1):
  # for ThreadChunkNumber in range( 1 ):
    Start = ThreadChunkNumber*Args.thread
    End = Start+Args.thread
    # print Start, End 
    PdbSubset = Pdbs[Start: End]
    
    SimulationTuple = []
    print 'PdbSubset:', PdbSubset 
    for i, Pdb in enumerate(PdbSubset):
      RepeatLength = int(re.sub(r'.*rep(\d+)_.*', r'\1', Pdb))
      SimulationTuple.append( (Pdb, CstHash[Pdb], RepeatLength, NativePdb) )

    # for InputTuple in SimulationTuple:
    #   # print 'Error in one of pooled multiprocessing threads; iterating sequentially for debugging '
    #   optimize_repeat_pdb(InputTuple)

    pool = Pool(processes=len(SimulationTuple))
    pool.map(optimize_repeat_pdb, SimulationTuple)


# if __name__ == "__main__":
#   sys.exit(main())
