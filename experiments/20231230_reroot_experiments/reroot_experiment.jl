
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


### test to see if i get it working properly
tree = sim_tree(n=10)
orig_tree = deepcopy(tree)
orig_tree.children[2].branchlength + orig_tree.children[2].children[1].branchlength
tree_draw(orig_tree, draw_labels=true, canvas_height=10cm, canvas_width=20cm, font_size=12pt)

orig_tree.children[2].branchlength + orig_tree.children[2].children[1].branchlength
tree = MolecularEvolution.recursive_reroot(orig_tree.children[2].children[1], dist_above_child=1)
tree_draw(tree, draw_labels=true, canvas_height=30cm, canvas_width=40cm, font_size=24pt)


