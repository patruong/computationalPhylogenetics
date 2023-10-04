using Pkg
Pkg.activate("environments/MolecularEvolution")
Pkg.status()


cd("git/computationalPhylogenetics")

include("codonFrequencies.jl")  # Include the script to load the module
include("MG94.jl")  # Include the script to load the module

import .CodonFrequencies
import .CodonFrequencies.generate_nucleotide_frequencies
import .MG94
import .MG94.compute_Q
import .MG94.compute_equilibrium_rate

using FASTX, MolecularEvolution, Measures, Compose, Phylo, DataFrames
using Plots
using StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling
using MolecularEvolution, Plots
using CodonMolecularEvolution
Pkg.add("Revise")
using Revise # for code dev


# Define the parameters for the MG94 model
α = 1.0 # synonymous rate 
β = 1.0 # Nonsynonymous rate
β = 0.927456919692643
#ω = 1.0  # Nonsynonymous/synonymous rate ratio
ω = β / α
κ = 2.0  # Transition/transversion rate ratio
debug = false # debug message true or false
remove_stop_codon = true

freq_list = [0.2, 0.3, 0.3, 0.2] # Nucleotide frequencies 
nucleotide_list = ['A', 'C', 'G', 'T'] # Nucleotides


nucleotide_frequencies = generate_nucleotide_frequencies(freq_list, nucleotide_list)


codon_frequency_dictionary = generate_codon_frequencies(freq_list, nucleotide_list, remove_stop_codon, debug)

# Equilibrium codon frequencies for 61 sense codons (excluding stop codons)
# You need to specify the actual values for your specific dataset
# For example, you can calculate these frequencies from your data or use reference values.
#equilibrium_frequencies = [0.05, 0.05, ..., 0.05]  # Replace with your values

# Define the equilibrium frequencies for the 61 sense codons
equilibrium_frequencies = fill(1.0 / 61.0, 61)

# These list includes stop frequencies
codon_list, codon_freq_list = generate_codon_frequencies_key_values(freq_list, nucleotide_list, remove_stop_codon, debug)

equilibrium_frequencies = codon_freq_list



# The Q matrix now represents the rate matrix for the MG94 model
Q = compute_Q(equilibrium_frequencies, α, β, κ, false)

t = 1.0  # Time (you can specify the desired time)
equilibrium_rates = compute_equilibrium_rate(Q, t)
# Print or use the equilibrium rates as needed
println("Equilibrium Rates: ", equilibrium_rates)
