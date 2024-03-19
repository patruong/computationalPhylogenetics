
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution
using CodonMolecularEvolution

cd("/home/patrick/git/computationalPhylogenetics/")

#function import_grouped_label_tree(tree_file)
# Takes a Newick tree file and return Newick tree, Newick tree with replaced tags, group tags, original tags, and randomly generated colours for each tag
#tree = read_newick_tree(tree_file)
#treestring = newick(tree)
#treestring_group_labeled, group_tags, original_tags = CodonMolecularEvolution.replace_newick_tags(treestring)
#tag_colors = generate_hex_colors(length(original_tags))
#    return treestring_group_labeled, treestring, group_tags, original_tags, tag_colors
#end


function get_unique_tags(json_file)
    json_data = JSON.parsefile(json_file)
    treestring = json_data["input"]["trees"]["0"]
    tags = json_data["tested"]["0"]
    return unique(values(tags))
end

function read_in_tree_from_hyphy_result_json(json_file)
    json_data = JSON.parsefile(json_file)
    treestring = json_data["input"]["trees"]["0"]
    tags = json_data["tested"]["0"]
    return treestring, tags
end

function tag_hyphy_tree(treestring, tags, group1, group2)
    pattern = r"[\(\),:]([^,:]+)(?=:)"
    matches = eachmatch(pattern, treestring)
    modified_treestring = treestring


    function get_last_match_from_regex(matches)
        last_match = nothing
        for match in matches
            last_match = match.match
        end
        return last_match
    end

    for match in matches
        original_string = match.captures[1]
        regex = r"[A-Za-z0-9_]+"
        pruned_match = eachmatch(regex, original_string)
        key = get_last_match_from_regex(pruned_match)
        #println(key)
        try
            tag = tags[key]
            if tag in [group1, group2]
                replacement = original_string * "{$tag}" # Note this needs to be original_string, otherwise we interpret "Node22" as "Node2"*"2", the original string contains some crap before the string which makes it kinda unique.
            else
                replacement = original_string
            end
            modified_treestring = replace(modified_treestring, original_string => replacement)
        catch
            println("No key $key")
        end
    end

    return modified_treestring
end


# function to read and tag trees without branchlenght
function tag_tree_without_branchlength(treestring, group1, group2, branchlength=0)
    """
    Do not use on a tree with branchlength as this will remove branchlength information
    """
    tree = gettreefromnewick(treestring, FelNode)
    MolecularEvolution.binarize!(tree)
    MolecularEvolution.ladderize!(tree)
    for n in getnodelist(tree)
        n.name = replace(n.name, "'" => "")
        n.branchlength = rand()  # Appends a colon and a random float to the node name for id 
    end
    treestring = newick(tree)
    ## special tagger
    # 2 group + background
    #fasta_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex"
    #json_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex.FEL.json"
    # seqnames, seqs = read_fasta(fasta_file)

    tagged_treestring = tag_hyphy_tree(treestring, tags, group1, group2)
    tree = gettreefromnewick(tagged_treestring, FelNode)
    for n in getnodelist(tree)
        n.name = replace(n.name, "'" => "")
        n.branchlength = branchlength  # set fix value for branchlength
    end
    treestring = newick(tree)
    return treestring
end

function write_string_to_file(filename, str)
    file = open(filename, "w")
    write(file, str)
    close(file)
end

function get_node_by_name(tree, nodename)
    mynode = []
    for node in getnodelist(tree)
        if node.name == nodename
            push!(mynode, node)
        end
    end
    mynode = mynode[1]
    return mynode
end

function get_group_from_name(input_string)
    group_pattern = r"\{(.+?)\}"
    match_result = match(group_pattern, input_string)
    if match_result != nothing
        return match_result.match
    else
        return ""
    end
end

function recursive_tagging(node)
    if isleafnode(node)
        return
    end


    left_child = recursive_tagging(node.children[1])
    right_child = recursive_tagging(node.children[2])

    if get_group_from_name(node.children[1].name) == get_group_from_name(node.children[2].name)
        group = get_group_from_name(node.children[1].name)
    else
        group = ""
    end

    # Modify the information of the current node based on its children
    node.name = node.name * group

end
#recursive_tagging(tree)

function remove_groups(input::AbstractString)
    # Regular expression to match and remove content inside {}
    result = replace(input, r"\{[^}]*\}" => "")
    return result
end


function plot_tree(tree, tags, tag_colors)
    p = sortperm(tags)
    tags, tag_colors = tags[p], tag_colors[p]
    push!(tag_colors, "black")
    strip_tags_from_name = CodonMolecularEvolution.generate_tag_stripper(tags)
    color_dict = Dict(zip(getnodelist(tree), [tag_colors[CodonMolecularEvolution.model_ind(n.name, tags)] for n in getnodelist(tree)]))
    label_dict = Dict(zip(getnodelist(tree), [strip_tags_from_name(n.name) for n in getnodelist(tree)]))

    tree_draw(tree, canvas_height=(3 + length(getleaflist(tree)) / 5)cm,
        draw_labels=true, dot_color_dict=color_dict,
        line_color_dict=color_dict, line_width=0.3, min_dot_size=0.01,
        nodelabel_dict=label_dict)
end

function hivRT_add_tag_to_leaf_dict(seqnames)
    mapping_dict = Dict{String,String}()
    for name in seqnames
        if startswith(name, 'T')
            mapping_dict[name] = "$name{TREATED}"
        elseif startswith(name, 'N')
            mapping_dict[name] = "$name{NAIVE}"
        end
    end
    return mapping_dict
end


function hivRT_tag_tree(tree)
    mapping_dict = hivRT_add_tag_to_leaf_dict(seqnames)
    for node in getnodelist(tree)
        if node.name in seqnames
            node.name = mapping_dict[node.name]
        end
    end
end


function gettreefromnewickstr(treestring::String;
    binarize=true,
    ladderize=true,
    strip_single_quotes=true,)
    tree = gettreefromnewick(treestring, FelNode)
    if binarize
        MolecularEvolution.binarize!(tree)
    end
    if ladderize
        MolecularEvolution.ladderize!(tree)
    end
    if strip_single_quotes
        for n in getnodelist(tree)
            n.name = replace(n.name, "'" => "")
        end
    end
    return tree
end


function hivRT_merge_trees(naive_tree, treated_tree, distance_to_root=0.01)
    naive_nwk = newick(naive_tree)
    treated_nwk = newick(treated_tree)
    naive_nwk_part = naive_nwk[1:end-1]
    treated_nwk_part = treated_nwk[1:end-1]
    #distance_to_root = 0.01
    nwk = "($(treated_nwk_part):$distance_to_root, $(naive_nwk_part):$distance_to_root)"
    treestring_group_labeled, group_tags, tags = CodonMolecularEvolution.replace_newick_tags(nwk)
    tag_colors = CodonMolecularEvolution.generate_hex_colors(length(tags))
    tree = gettreefromnewick(nwk, FelNode)
    return tree, nwk, treestring_group_labeled, group_tags, tags, tag_colors
end

# 2 group + background

fasta_file = "experiments/20240319_hivRT_rebuild/HIV_RT.fasta"
#json_file = "data/contrastFEL_empirical_data/hiv-1_reverse_transcriptase/HIV_RT.nex.FEL.json"
seqnames, seqs = read_fasta(fasta_file)


naive_tree = read_newick_tree("experiments/20240319_hivRT_rebuild/HIV_RT_NAIVE.fasta.treefile")
treated_tree = read_newick_tree("experiments/20240319_hivRT_rebuild/HIV_RT_TREATED.fasta.treefile")

hivRT_tag_tree(naive_tree)
hivRT_tag_tree(treated_tree)
recursive_tagging(naive_tree)
recursive_tagging(treated_tree)
tree, nwk, treestring_group_labeled, group_tags, tags, tag_colors = hivRT_merge_trees(naive_tree, treated_tree, 0.005)
plot_tree(tree, tags, tag_colors)

treestring = newick(tree)

write_string_to_file("experiments/20240319_hivRT_rebuild/HIV_RT_pure_clade.nwk", treestring)
