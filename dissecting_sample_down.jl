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
cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()


n(t) = exp(sin(t / 10) * 2.0 + 4)
s(t) = n(t) / 100
tree = sim_tree(500, n, s, mutation_rate=0.005)


iters = 2500
verbosity = 1
pos_thresh = 0.95

# write as Newick string file
newick_string = newick(tree)
open("flu_sim.tre", "w") do io
    println(io, newick_string)
end

# write as nexus file
MolecularEvolution.write_nexus("flu_sim.nex", tree)

nuc_freqs = [0.2, 0.3, 0.3, 0.2]
nuc_rates = [1.0, 2.0, 1.0, 1.0, 1.6, 0.5]

nuc_model = DiagonalizedCTMC(reversibleQ(nuc_rates, nuc_freqs))
bm_model = BrownianMotion(0.0, 1.0)

message_template = [NucleotidePartition(nuc_freqs, 300), GaussianPartition()]
internal_message_init!(tree, message_template)

sample_down!(tree, [nuc_model, bm_model])



