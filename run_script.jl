using Pkg
Pkg.activate("environments/MolecularEvolution")
Pkg.status()


cd("git/computationalPhylogenetics")

include("codonFrequencies.jl")  # Include the script to load the module
include("MG94.jl")  # Include the script to load the module

import .CodonFrequencies
import .MG94
import .MG94.compute_Q
import .MG94.compute_equilibrium_rate

using FASTX, MolecularEvolution, Measures, Compose, Phylo, DataFrames
using Plots
using StatsBase, Distributions, DataFrames, CSV, NLopt, ParameterHandling
using MolecularEvolution, Plots
using CodonMolecularEvolution


# Define the parameters for the MG94 model
α = 1.0 # synonymous rate 
β = 1.0 # Nonsynonymous rate
#ω = 1.0  # Nonsynonymous/synonymous rate ratio
ω = β / α
κ = 2.0  # Transition/transversion rate ratio
debug = true # debug message true or false
remove_stop_codon = true

freq_list = [0.2, 0.3, 0.3, 0.2] # Nucleotide frequencies 
nucleotide_list = ['A', 'C', 'G', 'T'] # Nucleotides

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


#### 

n(t) = exp(sin(t / 10) * 2.0 + 4)
s(t) = n(t) / 100
tree = sim_tree(500, n, s, mutation_rate=0.005)


# write as Newick string file
newick_string = newick(tree)
open("flu_sim.tre", "w") do io
    println(io, newick_string)
end

# write as nexus file
MolecularEvolution.write_nexus("flu_sim.nex", tree)

# read in sequence and fasta

#seqnames, seqs = read_fasta("Ace2_tiny_test.fasta");
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("Ace2_no_background.nex")

nuc_freqs = [0.2, 0.3, 0.3, 0.2]
nuc_rates = [1.0, 2.0, 1.0, 1.0, 1.6, 0.5]

# change nuc freqs to codon freqs and instantaneous rates
codon_list, codon_freq_list = generate_codon_frequencies_key_values(freq_list, nucleotide_list, remove_stop_codon, debug)
Q = compute_Q(equilibrium_frequencies, α, β, κ, false)

codon_freq = codon_freq_list
revQ = reversibleQ(nuc_rates, nuc_freqs)

nuc_model = DiagonalizedCTMC(reversibleQ(nuc_rates, nuc_freqs))
bm_model = BrownianMotion(0.0, 1.0)

message_template = [NucleotidePartition(nuc_freqs, 300), GaussianPartition()]
internal_message_init!(tree, message_template)

sample_down!(tree, [nuc_model, bm_model])

for n in getnodelist(tree)
    n.node_data = Dict(["mu" => n.message[2].mean])
end

MolecularEvolution.write_nexus("flu_sim_after_sample_down.nex", tree)

phylo_tree = get_phylo_tree(tree)


plot(phylo_tree, showtips=false, line_z="mu", colorbar=:none,
    linecolor=:darkrainbow, linewidth=1.0, size=(600, 600))

df = DataFrame()
df.names = [n.name for n in getleaflist(tree)]
df.seqs = [partition2obs(n.message[1]) for n in getleaflist(tree)]
df.mu = [partition2obs(n.message[2]) for n in getleaflist(tree)]
CSV.write("flu_sim_seq_and_bm.csv", df)

write_fasta("flu_sim_seq_and_bm.fasta", df.seqs, seq_names=df.names)


##### Code for molecularEvo Q_ij
#using Pkg
#using Random
#using Plots


#seed_value = 123
#Pkg.activate()
#Pkg.activate("~/.julia/environments/CodonMolecularEvolution")
#Pkg.activate("environments/CodonMolecularEvolution")
#Pkg.status()

cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()

#Pkg.add("Revise")
#using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment

using MolecularEvolution, FASTX, CodonMolecularEvolution
#using Test

# We need a fasta file
# and a .nex tree file.

analysis_name = "test_analysis"
seqnames, seqs = read_fasta("flu_sim_seq_and_bm.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("flu_sim.nex") # can this function import .tre instead of .nex file?
@time df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name);

iters = 2500
verbosity = 1
pos_thresh = 0.95

CodonMolecularEvolution.generate_tag_stripper
code = MolecularEvolution.universal_code
#tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors)
tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code)
con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4)


alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)

#df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=false)


nuc_matrix = GTRmat
F3x4 = F3x4_freqs
Q_mol = MolecularEvolution.MG94_F3x4(alpha, beta, nuc_matrix, F3x4) # MG94 model from Molecular MolecularEvolution

### MATCH Q_mol with simulated data in a new less messy script!!

