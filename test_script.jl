using Pkg
using Random
using Plots


seed_value = 123
#Pkg.activate()
#Pkg.activate("~/.julia/environments/CodonMolecularEvolution")
Pkg.activate("environments/CodonMolecularEvolution")
Pkg.status()

cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()

Pkg.add("Revise")
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment

using MolecularEvolution, FASTX, CodonMolecularEvolution
#using Test


# difFUBAR_test.jl
#analysis_name = "nobackground/Ace2"
#seqnames, seqs = read_fasta("data/larger/Ace2_tiny_tagged.fasta");
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

analysis_name = "nobackground/Ace2_test"
seqnames, seqs = read_fasta("data/larger/Ace2_tiny_tagged.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/larger/Ace2_no_background.nex")
@time df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name);

iters = 2500
verbosity = 1
pos_thresh = 0.95
# where is the model? 
CodonMolecularEvolution.difFUBAR_global_fit

CodonMolecularEvolution.generate_tag_stripper
code = MolecularEvolution.universal_code
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors)
tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code)
con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4)


alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)

#df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=false)


nuc_matrix = GTRmat
F3x4 = F3x4_freqs
MolecularEvolution.MG94_F3x4(alpha, beta, nuc_matrix, F3x4) # MG94 model from Molecular MolecularEvolution

##### THIS IS THE RUN code

nuc_freqs = [0.2, 0.3, 0.3, 0.2]
nuc_rates = [1.0, 2.0, 1.0, 1.0, 1.6, 0.5]
reversibleQ(nuc_rates, nuc_freqs)

nuc_model = DiagonalizedCTMC(reversibleQ(nuc_rates, nuc_freqs))
nuc_model = DiagonalizedCTMC(GTRmat)


codon_model =
    F3x4_freqs
bm_model = BrownianMotion(0.0, 1.0)

message_template = [NucleotidePartition(nuc_freqs, 300), GaussianPartition()]
internal_message_init!(tree, message_template)

sample_down!(tree, [nuc_model, bm_model])

for n in getnodelist(tree)
    n.node_data = Dict(["mu" => n.message[2].mean])
end

phylo_tree = get_phylo_tree(tree)
plot(phylo_tree, showtips=false, line_z="mu", colorbar=:none,
    linecolor=:darkrainbow, linewidth=1.0, size=(600, 600))

df = DataFrame()
df.names = [n.name for n in getleaflist(tree)]
df.seqs = [partition2obs(n.message[1]) for n in getleaflist(tree)]
df.mu = [partition2obs(n.message[2]) for n in getleaflist(tree)]
CSV.write("flu_sim_seq_and_bm.csv", df)

write_fasta("flu_sim_seq_and_bm.fasta", df.seqs, seq_names=df.names)


