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


getleaflist(tree)[371].parent.children[2].name
getleaflist(tree)[371].parent.children[1].name
getleaflist(tree)[371].parent.name

###

comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

# Printout for information
group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

# group_name printout

tag_pos = [1, 2]
tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
tags, tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)

treestring = newick(tree)


###
tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init("analysis_name", treestring, tags, tag_colors, exports=true, verbosity=1)
pwd()