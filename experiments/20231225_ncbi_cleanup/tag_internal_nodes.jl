using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using ArgParse

#cd("/home/patrick/git/computationalPhylogenetics/experiments/20231225_ncbi_cleanup/output/tree")
#newick_file = "parvo_panleu_trans_tagged_fixed_aligned.nwk"
#output_file = "output.nwk"

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
        help = "Path to the input newick file"
        arg_type = String
        "--output", "-o"
        arg_type = String
        help = "Output newick file"
        default = "output.nwk"
    end

    return parse_args(s)
end

function main()

    args = parse_commandline()

    newick_file = args["input"]
    output_file = args["output"]

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

    treestring = newick(tree)
    # write to file
    file = open(output_file, "w")
    write(file, treestring)
    close(file)
end

main()