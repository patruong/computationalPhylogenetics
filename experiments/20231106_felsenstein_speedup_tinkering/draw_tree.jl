using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo



cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
analysis_name = "Ace2"
seqnames, seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=false, verbosity=1)

img = tree_draw(tree)
img |> PDF("imgout_.pdf", 10cm, 10cm)

############
# Sim tree #
############

using Pkg
Pkg.activate("environments/CodonMolecularEvolution")
#using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
#using MolecularEvolution, FASTX, CodonMolecularEvolution

#using Pkg
#Pkg.add("Revise")
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution



file_dir = "/home/ptruong/git/computationalPhylogenetics/contrastFEL_data/omnibus-multi/datasets/omnibus-multi/"
fasta = file_dir * "sims.0.settings.replicate.2"
tree_file = file_dir * "sims.0.nwk"
seqnames, seqs = import_hyphy_simulated_FASTA(fasta)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)

comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

# Printout for information
group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

# group_name printout
group_meta_data_file_name = "group_meta_data.txt"
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
analysis_name = "/home/ptruong/git/computationalPhylogenetics/experiments/20231027_tinker_about_efficiency/output"
analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(group_tags, tag_colors, tag_pos)
cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)

#difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name;
#    pos_thresh=0.95, iters=2500, verbosity=1,
#    exports=true, code=MolecularEvolution.universal_code)

@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, analysis_tags, analysis_tag_colors, exports=false, verbosity=1)

img = tree_draw(tree)
img |> PDF("imgout_sim.pdf", 10cm, 30cm)



typeof(tree)
tree
tree.parent
tree.children[1].name
tree.children[1].nodeindex
tree.children[1]


tree



