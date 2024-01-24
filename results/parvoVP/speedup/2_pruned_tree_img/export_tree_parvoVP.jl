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

cd("/home/patrick/git/computationalPhylogenetics/results/parvoVP/speedup/2_pruned_tree_img")
fasta_file = "../../ParvoVP.fasta"
newick_file = "../1_reassign_tree/ParvoVP_regrouped.nwk"
output_csv = "output_"
analysis_name = "analysis"
exports = true
iters = 2500
verbosity = 1

seqnames, seqs = read_fasta(fasta_file)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)


function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end

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
con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid_pruned_3(tree, tags, GTRmat, F3x4_freqs, MolecularEvolution.universal_code,
    verbosity=verbosity, foreground_grid=6, background_grid=4)

strip_tags_from_name = CodonMolecularEvolution.generate_tag_stripper(tags)
color_dict = Dict(zip(getnodelist(tree), [tag_colors[model_ind(n.name, tags)] for n in getnodelist(tree)]))
label_dict = Dict(zip(getnodelist(tree), [strip_tags_from_name(n.name) for n in getnodelist(tree)]))
img = tree_draw(tree, canvas_width=30cm, canvas_height=(3 + length(getleaflist(tree)) / 5)cm,
    draw_labels=true, dot_color_dict=color_dict,
    line_color_dict=color_dict, line_width=0.3, min_dot_size=0.01,
    nodelabel_dict=label_dict)

#img |> SVG(analysis_name * "_tagged_input_tree.svg", 15cm, (3 + length(getleaflist(tree)) / 5)cm)
img |> SVG(analysis_name * "_pruned_tagged_input_tree.svg", 30cm, (3 + length(getleaflist(tree)) / 5)cm)