using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Revise
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240605_MERS_COVID/results")

analysis_name = "difFUBAR/results"
fasta_file = "result.fasta"
nexus_file = "result.newick"
newick_file = "result.newick"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95

seqnames, seqs = read_fasta(fasta_file)
tree = read_newick_tree(newick_file)


tree.children[1].name
tree.children[2].name

treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_labelled_newick_tree(newick_file, true)

CodonMolecularEvolution.plot_tree(tree, original_tags, tag_colors)
node_name = "NC_007362"
node_name = "0/52"
node_name = "62.5/99"

node = CodonMolecularEvolution.get_node_by_name(tree, node_name)
dist_above_child = 0
tree = MolecularEvolution.recursive_reroot(node, dist_above_child=dist_above_child)
CodonMolecularEvolution.plot_tree(tree, original_tags, tag_colors)


tree.name
tree.children[2].children[1].name

# What is going on with the rerooting?
CodonMolecularEvolution.plot_tree(tree, original_tags, tag_colors)

tree = MolecularEvolution.recursive_reroot(node, dist_above_child=dist_above_child)

###
treestring_group_labeled, treestring, group_tags, tags, tag_colors = import_grouped_labelled_newick_tree_reroot(newick_file, "NC_007362", 0.01, true)

df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)



