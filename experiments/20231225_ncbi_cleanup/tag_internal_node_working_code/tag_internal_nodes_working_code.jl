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

cd("/home/patrick/git/computationalPhylogenetics/experiments/20231225_ncbi_cleanup/output/tree")
newick_file = "parvo_panleu_trans_tagged_fixed_aligned.nwk"
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)
tree = gettreefromnewick(treestring, FelNode)


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



getleaflist(tree)[1].parent.parent.name
getleaflist(tree)[1].parent.parent.children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2]
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[1].name
getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].name

getleaflist(tree)[1].parent.parent.children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].children[2].name


def propagate_info(node):
    # Base case: If the node is a leaf, update its information
    if node is None:
        return 0

    # Recursive step: Propagate information from children to the current node
    left_info = propagate_info(node.left)
    right_info = propagate_info(node.right)

    # Modify the information of the current node based on its children
    node.info = left_info + right_info + node.value

    return node.info






















































function add_group_to_node(node)
    if get_group_from_name(node.children[1].name) == get_group_from_name(node.children[2].name)
        group = get_group_from_name(node.children[1].name)
    else
        group = ""
    end
    return node.name * group
end







node = getleaflist(tree)[2].parent
child_groups = [get_group_from_name(child.name) for child in node.children]


felsenstein!

###################3



getleaflist(tree)
get_group_from_name(input_string)

parent_list = []
for leaf in getleaflist(tree)
    leaf_group = get_group_from_name(leaf.name)
    sibling_group = get_group_from_name(siblings(leaf)[1].name)
    if leaf_group == sibling_group && !(leaf.parent in parent_list)
        leaf.parent.name = get_group_from_name(leaf.name)
        push!(parent_list, leaf.parent)
    end
end

l = [1, 2, 3, 4, 5]
pop!(l)
push!(l, 6)
l

length(parent_list)

for node in parent_list
    node_group = get_group_from_name(node.name)
    sibling_group = get_group_from_name(siblings(node)[1].name)
    if node_group == sibling_group && !(leaf.parent in parent_list)
        leaf.parent.name = get_group_from_name(leaf.name)
        push!(parent_list, leaf.parent)
    end
end







node = getleaflist(tree)[1].parent








node.name = add_group_to_node(node)


# TEST
node = getleaflist(tree)[1]
node.parent.parent.children[1].name
node.parent.parent.children[2].children[1].name

node.parent.parent.children[2].children[1].children[1].name
node.parent.parent.children[2].children[1].children[2].name

node.parent.parent.children[2].children[2]


