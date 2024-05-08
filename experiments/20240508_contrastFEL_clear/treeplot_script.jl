using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse
using CSV
using DataFrames
using Tables

cd("/home/patrick/git/computationalPhylogenetics/contrastFEL_data/omnibus")

exports = true
verbosity = 1
sim = 279

fasta_file = "sims." * string(sim) * ".settings.replicate.1"
newick_file = "sims." * string(sim) * ".nwk"
output_csv = "sims." * string(sim) * ".output"
analysis_name = "output/sims." * string(sim)
seqnames, seqs = import_hyphy_simulated_FASTA(fasta_file) #for sim_tree_fixed
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)


cd("/home/patrick/git/computationalPhylogenetics/experiments/20240508_contrastFEL_clear")
comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

# Printout for information
group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]



tag_pos = [1, 2]
tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
tags, tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)

tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)



#tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), MolecularEvolution.universal_code, verbosity=verbosity)
#GTRmat # what is going on with this?

#df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)
# Save DataFrame to CSV