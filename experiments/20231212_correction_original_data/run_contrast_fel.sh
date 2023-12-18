#!/bin/bash

# Raw run code
#time hyphy contrast-fel --alignment  data/sim.replicate.1.nex --tree  data/sim.replicate.1.nex --branch-set TEST --p-value 1.00 --q-value 1.00 
#time hyphy contrast-fel --alignment  data/HIVtransmission.nex --tree  data/HIVtransmission.nex --branch-set SOURCE --p-value 1.00 --q-value 1.00 

# Timed run code
# working
#(time hyphy contrast-fel --alignment data/sim_tree_fixed/sim.replicate.1.nex --tree data/sim_tree_fixed/sim.replicate.1.nex --branch-set TEST --p-value 1.00 --q-value 1.00) 2> time_output.txt > hyphy_output.txt

#(time hyphy contrast-fel --alignment ParvoVP_removed_duplicates.fasta  --tree ParvoVP_datamonkey_add_group.nwk --branch-set FELINE --branch-set CANINE --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt
#(time hyphy contrast-fel --alignment ParvoVP.fasta  --tree ParvoVP.nex --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt

#(time hyphy contrast-fel --alignment ParvoVP_removed_duplicates.fasta  --tree ParvoVP_datamonkey_add_group.nwk --branch-set FELINE --branch-set CANINE --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt
(time hyphy contrast-fel --alignment ParvoVP_removed_duplicates.fasta  --tree ParvoVP_removed_duplicates.nex --branch-set G1 --branch-set G2 --p-value 1.00 --q-value 1.00) 2>&1 | tee time_output.txt hyphy_output.txt


#debug this
#time hyphy contrast-fel --alignment ParvoVP.fasta  --tree ParvoVP.nex --p-value 1.00 --q-value 1.00 # not unique sequences in ParvoVP.fasta
time hyphy contrast-fel --alignment ParvoVP_removed_duplicates.fasta   --tree ParvoVP_removed_duplicates.nex --p-value 1.00 --q-value 1.00 # now we have unique sequences, but 371 terminal nodes






