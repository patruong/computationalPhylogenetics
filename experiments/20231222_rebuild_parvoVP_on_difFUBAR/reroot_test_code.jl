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


dir = "/home/patrick/git/computationalPhylogenetics/experiments/20231222_rebuild_parvoVP_on_difFUBAR/"
cd(dir)
analysis_name = "output_reroot_test"
fasta_file = dir * "parvo_panleu_trans.fasta"
newick_file = dir * "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped.nwk"
#fasta_file = dir * "sim_tree_fixed/sim.replicate.1.nex"
#newick_file = dir * "sim_tree_fixed/sim_tree_fixed.nwk"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
output_csv = "test/"


seqnames, seqs = read_fasta(fasta_file)
#seqnames, seqs = import_hyphy_simulated_FASTA(fasta_file) #for sim_tree_fixed

treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)


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

# bottom feline node
getleaflist(tree)[371].name

feline_node_MW495836 = getleaflist(tree)[371]



feline_node_MW495836.parent.name
feline_node_MW495836.parent.parent.name
feline_node_MW495836.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.name #this is what we want as root here
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.branchlength

feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.branchlength
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.branchlength
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.branchlength
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.branchlength
feline_node_MW495836.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.parent.branchlength

treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init("reroot", treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
feline_node_MW495836 = getleaflist(tree)[371]
feline_node_MW495836.branchlength
feline_node_MW495836.parent.branchlength
treestring = MolecularEvolution.newick(tree)
MolecularEvolution.reroot!(feline_node_MW495836, dist_above_child=0.00001)
treestring = MolecularEvolution.newick(tree)
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init("reroot", treestring, tags, tag_colors, exports=exports, verbosity=verbosity)

# dont understand how the reroot! works....

# try work with simulated data to se what happens.... from MolecularEvolution.jl


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


tree = sim_tree(n=10)
orig_tree = deepcopy(tree)
tree_draw(orig_tree, draw_labels=true, canvas_height=10cm, canvas_width=20cm, font_size=12pt)

tree = MolecularEvolution.reroot!(orig_tree.children[2], dist_above_child=orig_tree.children[2].branchlength / 2)
tree = MolecularEvolution.reroot(orig_tree.children[2], dist_above_child=orig_tree.children[2].branchlength / 2)
tree = MolecularEvolution.recursive_reroot(orig_tree.children[2], dist_above_child=4)
tree = MolecularEvolution.recursive_reroot!(orig_tree.children[2], dist_above_child=4)

tree_draw(tree, draw_labels=true, canvas_height=10cm, canvas_width=20cm, font_size=12pt)
# reroot function is not working... remake the function

orig_tree.children[2].branchlength

# MolecularEvolution reroot functions dont work properly
# reroot! - copies and reroots
# reroot - does not se dist_above_child
# recursive_reroot! - cuts out the subtree that is in children[2]
# recursive_reroot - I think this works!
