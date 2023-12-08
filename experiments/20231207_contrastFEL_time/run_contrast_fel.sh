#!/bin/bash

# Raw run code
#time hyphy contrast-fel --alignment  data/sim.replicate.1.nex --tree  data/sim.replicate.1.nex --branch-set TEST --p-value 1.00 --q-value 1.00 
#time hyphy contrast-fel --alignment  data/HIVtransmission.nex --tree  data/HIVtransmission.nex --branch-set SOURCE --p-value 1.00 --q-value 1.00 

# Timed run code
# working
#(time hyphy contrast-fel --alignment data/sim_tree_fixed/sim.replicate.1.nex --tree data/sim_tree_fixed/sim.replicate.1.nex --branch-set TEST --p-value 1.00 --q-value 1.00) 2> time_output.txt > hyphy_output.txt
#(time hyphy contrast-fel --alignment data/sim/sim.replicate.1.nex --tree data/sim/sim.replicate.1.nex --branch-set G1 --p-value 1.00 --q-value 1.00) 2> time_output.txt > hyphy_output.txt
#(time hyphy contrast-fel --alignment data/ParvoVP.nex --tree data/ParvoVP.fasta --branch-set TEST --p-value 1.00 --q-value 1.00) 2> time_output.txt > hyphy_output.txt
#(time hyphy contrast-fel --alignment data/ParvoVP.fasta --tree data/ParvoVP_modified_no_color.nex --branch-set FELINE --branch-set CANINE --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt
#(time hyphy contrast-fel --alignment data/ParvoVP.fasta --tree data/ParvoVP_modified_no_color.nex --branch-set FELINE --branch-set CANINE --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt

(time hyphy contrast-fel --alignment data/ParvoVP.fasta  --tree data/datamonkey_tree/tree.nex --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt

