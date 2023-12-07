#!/bin/bash

# Raw run code
#time hyphy contrast-fel --alignment  data/sim.replicate.1.nex --tree  data/sim.replicate.1.nex --branch-set TEST --p-value 1.00 --q-value 1.00 
#time hyphy contrast-fel --alignment  data/HIVtransmission.nex --tree  data/HIVtransmission.nex --branch-set SOURCE --p-value 1.00 --q-value 1.00 

# Timed run code
(time hyphy contrast-fel --alignment data/sim_tree_fixed/sim.replicate.1.nex --tree data/sim_tree_fixed/sim.replicate.1.nex --branch-set TEST --p-value 1.00 --q-value 1.00) 2> time_output.txt > hyphy_output.txt
#(time hyphy contrast-fel --alignment data/ParvoVP.fasta --tree data/ParvoVP.nex --p-value 1.00 --q-value 1.00) 2> time_output.txt > hyphy_output.txt
#(time hyphy contrast-fel --alignment data/ParvoVP.fasta --tree data/ParvoVP.nex --p-value 1.00 --q-value 1.00) #2> time_output.txt > hyphy_output.txt
#(time hyphy contrast-fel --alignment data/sim_reference_test/sim.replicate.1.nex --tree data/sim_reference_test/sim.replicate.1.nex --p-value 1.00 --q-value 1.00) #2> time_output.txt > hyphy_output.txt
