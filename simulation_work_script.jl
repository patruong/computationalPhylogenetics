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
using LinearAlgebra
using Revise # for code dev
cd("/home/ptruong/git/computationalPhylogenetics/workflow_pipeline/output")
pwd()


function nucleotide_equilibrium_matrix(diagonal_values::Vector{Float64}=[0.25, 0.25, 0.25, 0.25])
    D = Diagonal(diagonal_values)
    R = Matrix(D)
    return R
end

function nucleotide_rate_matrix(values)
    if length(values) != 6
        throw(ArgumentError("Input vector must have 6 elements"))
    end

    a, b, c, d, e, f = values

    # Create the matrix
    A = [-a-b-c a b c
        a -a-d-e d e
        b d -b-d-f f
        c e f -c-e-f]

    return A
end

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


# Nucleotide model
nuc_freqs = [0.2, 0.3, 0.3, 0.2]
nuc_rates = [1.0, 2.0, 1.0, 1.0, 1.6, 0.5]
nuc_model = DiagonalizedCTMC(reversibleQ(nuc_rates, nuc_freqs))
bm_model = BrownianMotion(0.0, 1.0)


# Codon model
# F3x4, row sums need to be equal to 1,
# Since we want to simulate the sequences we need to set these.
F3x4 = [0.20 0.30 0.30 0.20
    0.22 0.30 0.30 0.18
    0.20 0.23 0.35 0.22]
F61 = MolecularEvolution.F3x4_eq_freqs(F3x4) # codon equilibrium frequencies

# GTR matrix (instantaneous rate computations)
µ = 1
Π = nucleotide_equilibrium_matrix([0.25, 0.25, 0.25, 0.25])
R = nucleotide_rate_matrix([1, 2, 3, 4, 5, 6])
Q = µ * R * Π # GTR matrix

Q_MG94_F3x4 = MolecularEvolution.MG94_F3x4(alpha, beta, nuc_matrix, F3x4) # MG94 model from Molecular MolecularEvolution
Q_MG94_F61 = MolecularEvolution.MG94_F61(alpha, beta, nuc_matrix, F61)
Q_HB98_F61 = MolecularEvolution.HB98_F61(alpha, nuc_matrix, F61)

codon_model_MG94_F3x4 = DiagonalizedCTMC(Q_MG94_F3x4)
codon_model_MG94_F61 = DiagonalizedCTMC(Q_MG94_F61)
codon_model_HB98_F61 = DiagonalizedCTMC(Q_HB98_F61)

codon_model_MG94_F3x4 = GeneralCTMC(Q_MG94_F3x4)
codon_model_MG94_F61 = GeneralCTMC(Q_MG94_F61)
codon_model_HB98_F61 = GeneralCTMC(Q_HB98_F61)

message_template = [NucleotidePartition(nuc_freqs, 300), GaussianPartition(), CodonPartition(300), CodonPartition(300), CodonPartition(300)]
internal_message_init!(tree, message_template)
sample_down!(tree, [nuc_model, bm_model, codon_model_MG94_F3x4, codon_model_MG94_F61, codon_model_HB98_F61])

for n in getnodelist(tree)
    n.node_data = Dict(["mu" => n.message[2].mean])
end

x = getnodelist(tree)[1]
x.message[2].mean #this is from the gaussian partition


MolecularEvolution.write_nexus("flu_sim_after_sample_down.nex", tree)

phylo_tree = get_phylo_tree(tree)


plot(phylo_tree, showtips=false, line_z="mu", colorbar=:none,
    linecolor=:darkrainbow, linewidth=1.0, size=(600, 600))

df = DataFrame()
df.names = [n.name for n in getleaflist(tree)]
df.nuc_model = [partition2obs(n.message[1]) for n in getleaflist(tree)]
df.mu = [partition2obs(n.message[2]) for n in getleaflist(tree)]
df.MG94_F3x4 = [partition2obs(n.message[3]) for n in getleaflist(tree)]
df.MG94_F61 = [partition2obs(n.message[4]) for n in getleaflist(tree)]
df.HB98_F61 = [partition2obs(n.message[5]) for n in getleaflist(tree)]

CSV.write("flu_sim_codon.csv", df)
write_fasta("flu_sim_nuc_model.fasta", df.nuc_model, seq_names=df.names)
write_fasta("flu_sim_MG94_F3x4.fasta", df.MG94_F3x4, seq_names=df.names)
write_fasta("flu_sim_MG94_F61.fasta", df.MG94_F61, seq_names=df.names)
write_fasta("flu_sim_HB98_F61.fasta", df.HB98_F61, seq_names=df.names)

