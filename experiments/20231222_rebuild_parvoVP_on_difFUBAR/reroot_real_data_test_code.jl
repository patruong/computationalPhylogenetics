using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse

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

treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)
tag_pos = [1, 2]
tags, tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init("root", treestring, tags, tag_colors, exports=exports, verbosity=verbosity)

tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)
found_nodes = []
for node in getnodelist(tree)
    if node.name == "0.997{G1}"
        push!(found_nodes, node)
    end
end

reroot_node = found_nodes[1]

reroot_node.branchlength
# dist_above_child=reroot_node.branchlength might cause issue, the other reroot functions in MolecularEvolution.jl does not work.
#tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=reroot_node.branchlength)
#tree.children[1].branchlength
#tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=tree.children[1].branchlength)
#tree.children[2].branchlength
#tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=tree.children[2].branchlength)
#reroot_node.parent.branchlength
#tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=reroot_node.parent.branchlength)
reroot_node.branchlength / 2
tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=reroot_node.branchlength / 2)
#reroot_node.branchlength
#tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=reroot_node.branchlength)


tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)
treestring = MolecularEvolution.newick(tree)
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init("reroot", treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)

# we need to run group internal nodes here again to correct it!
function get_group_from_name(input_string)
    group_pattern = r"\{(.+?)\}"
    match_result = match(group_pattern, input_string)
    if match_result != nothing
        return match_result.match
    else
        return ""
    end
end

function propagate_info(node)
    if isleafnode(node)
        return
    end


    left_child = propagate_info(node.children[1])
    right_child = propagate_info(node.children[2])

    if get_group_from_name(node.children[1].name) == get_group_from_name(node.children[2].name)
        group = get_group_from_name(node.children[1].name)
    else
        group = ""
    end

    # Modify the information of the current node based on its children
    node.name = node.name * group

end

propagate_info(tree)

function remove_groups(input::AbstractString)
    # Regular expression to match and remove content inside {}
    result = replace(input, r"\{[^}]*\}" => "")
    return result
end

# remove original branching
tree.children[1].name = remove_groups(tree.children[1].name)
tree.children[2].name = remove_groups(tree.children[2].name)

treestring = MolecularEvolution.newick(tree)
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init("regroup", treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)


file_path = "output.txt"

file = open(file_path, "w")
write(file, my_string)
close(file)
println("Reroot and re-assigned branch newick tree has been written to $file_path.")