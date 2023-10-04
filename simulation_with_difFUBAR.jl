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
using Revise # for code dev


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
df.seqs = [partition2obs(n.message[1]) for n in getleaflist(tree)]
df.mu = [partition2obs(n.message[2]) for n in getleaflist(tree)]

CSV.write("flu_sim_seq_and_bm.csv", df)

write_fasta("flu_sim_seq_and_bm.fasta", df.seqs, seq_names=df.names)

##### CodonMolecularEvo


seqnames, seqs = read_fasta("flu_sim_seq_and_bm.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("flu_sim_after_sample_down.nex") # can this function import .tre instead of .nex file?

code = MolecularEvolution.universal_code
code
tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code)

con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4)


nuc_matrix = GTRmat # nucelotide substitution rate 
F3x4 = F3x4_freqs
F61 = eq_freqs # 61 codon equilibrium frequencies
Q_MG94_F3x4 = MolecularEvolution.MG94_F3x4(alpha, beta, nuc_matrix, F3x4) # MG94 model from Molecular MolecularEvolution
Q_MG94_F61 = MolecularEvolution.MG94_F61(alpha, beta, nuc_matrix, F61)
Q_HB98_F61 = MolecularEvolution.HB98_F61(alpha, nuc_matrix, F61)

Q_MG94_F3x4
Q_MG94_F61
Q_HB98_F61




codon_model_MG94_F3x4 = DiagonalizedCTMC(Q_MG94_F3x4)
codon_model_MG94_F61 = DiagonalizedCTMC(Q_MG94_F61)
codon_model_HB98_F61 = DiagonalizedCTMC(Q_HB98_F61)

codon_model_MG94_F3x4 = GeneralCTMC(Q_MG94_F3x4)
codon_model_MG94_F61 = GeneralCTMC(Q_MG94_F61)
codon_model_HB98_F61 = GeneralCTMC(Q_HB98_F61)

message_template = [CodonPartition(300), CodonPartition(300), CodonPartition(300)]
internal_message_init!(tree, message_template)
sample_down!(tree, [codon_model_MG94_F3x4, codon_model_MG94_F61, codon_model_HB98_F61])

#for n in getnodelist(tree)
#n.node_data = Dict(["mu" => n.message[2].mean])
#    println(n)
#end

x = getnodelist(tree)[1]
x.message


MolecularEvolution.write_nexus("flu_sim_after_sample_down.nex", tree)

phylo_tree = get_phylo_tree(tree)


plot(phylo_tree, showtips=false, line_z="mu", colorbar=:none,
    linecolor=:darkrainbow, linewidth=1.0, size=(600, 600))

df = DataFrame()
df.names = [n.name for n in getleaflist(tree)]
df.MG94_F3x4 = [partition2obs(n.message[1]) for n in getleaflist(tree)]
df.MG94_F61 = [partition2obs(n.message[2]) for n in getleaflist(tree)]
df.HB98_F61 = [partition2obs(n.message[3]) for n in getleaflist(tree)]

CSV.write("flu_sim_codon.csv", df)
write_fasta("flu_sim_MG94_F3x4.fasta", df.MG94_F3x4, seq_names=df.names)
write_fasta("flu_sim_MG94_F61.fasta", df.MG94_F61, seq_names=df.names)
write_fasta("flu_sim_HB98_F61.fasta", df.HB98_F61, seq_names=df.names)


########



# I'm doing something wrong because the Q and Q_mol are not really the same... check the parameters
# Check sample_down!





F3x4
F3x4[p[2][1], p[2][3]]


genetic_code = MolecularEvolution.GeneticCode(MolecularEvolution.universal_genetic_code)

for p in genetic_code.syn_positions
    println("$(p): $(F3x4[p[2][1], p[2][3]])")
end

for p in genetic_code.syn_positions
    println("$(p): $(nuc_matrix[p[2][2], p[2][3]])")
end

for p in genetic_code.syn_positions
    println("$(p): $(p[2][2])")
end


genetic_code.syn_positions
F3x4[p[2][1], p[2][3]]


nuc_matrix
genetic_code.codons





