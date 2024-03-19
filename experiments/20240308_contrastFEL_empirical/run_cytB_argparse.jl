
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution
using CodonMolecularEvolution
using ArgParse

cd("/home/patrick/git/computationalPhylogenetics/")


# Set threads
using Base.Threads
ENV["JULIA_NUM_THREADS"] = "20"


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
function tag_tree_without_branchlength(treestring, tags, group1, group2, branchlength=1)
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



function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--group1", "-a"
        help = "Group1 input"
        arg_type = String
        "--group2", "-b"
        help = "Group2 input"
        arg_type = String
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    group1 = args["group1"]
    group2 = args["group2"]



    # 4 group + background # This is a bit harder to tag
    fasta_file = "data/contrastFEL_empirical_data/cytochrome_B_of_Haemosporidians/Cytb.fasta"
    json_file = "data/contrastFEL_empirical_data/cytochrome_B_of_Haemosporidians/Cytb.fasta.FEL.json"
    seqnames, seqs = read_fasta(fasta_file)
    cleaned_seqs = []
    for i in 1:length(seqs)
        seq = replace(seqs[i], " " => "")
        push!(cleaned_seqs, seq)
    end
    seqs = cleaned_seqs
    seqs = [string(x) for x in seqs if isa(x, String)]

    treestring, tags = read_in_tree_from_hyphy_result_json(json_file)
    treestring = tag_tree_without_branchlength(treestring, tags, group1, group2, 1)
    treestring_group_labeled, group_tags, tags = CodonMolecularEvolution.replace_newick_tags(treestring)
    tag_colors = CodonMolecularEvolution.generate_hex_colors(length(tags))

    ################
    # Run difFUBAR #
    ################

    exports = true
    verbosity = 1
    iters = 2500
    pos_thresh = 0.95
    analysis_name = "cytochrome_B/" * group1 * "_" * group2 * "/analysis"
    df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)

end


main()