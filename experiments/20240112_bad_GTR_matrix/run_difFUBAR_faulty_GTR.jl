using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse
using CSV
using DataFrames
using MolecularEvolution, Plots, Phylo

cd("/home/patrick/git/computationalPhylogenetics/experiments/20240112_bad_GTR_matrix")
# params 

exports = true
verbosity = 1
fasta = "ParvoVP_removed_duplicates.fasta" #parsed_args["FASTAfile"]
#tree_file = "ParvoVP.nex" #parsed_args["treefile"]
tree_file = "ParvoVP_removed_duplicates.nwk" #parsed_args["treefile"]
analysis_name = "output" #parsed_args["output_dir"] * "/" * "ParvoVP"
seqnames, seqs = read_fasta(fasta)
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)


comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

# Printout for information
group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

# group_name printout
group_meta_data_file_name = output_csv * "group_meta_data.txt"
group_meta_data_file = open(group_meta_data_file_name, "w")
for element in group_name
    println(group_meta_data_file, "Groups")
    println(group_meta_data_file, element)
end
println(group_meta_data_file, "")

# pairwise comparison printout
for element in pairwise_comparisons_names
    println(group_meta_data_file, "Comparisons")
    println(group_meta_data_file, element)
end
close(group_meta_data_file)


tag_pos = [1, 2]
tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
tags, tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)

tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), MolecularEvolution.universal_code, verbosity=verbosity)
GTRmat # what is going on with this?



