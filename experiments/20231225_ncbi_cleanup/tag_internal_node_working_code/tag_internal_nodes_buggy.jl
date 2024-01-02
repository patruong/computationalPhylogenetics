using ArgParse
using MolecularEvolution
using CodonMolecularEvolution

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


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input_tree", "-i"
        help = "Path to the input newick file"
        arg_type = String
        "--output_tree", "-o"
        help = "Path the the output newick file"
        arg_type = String
        default = "output.nwk"
    end

    return parse_args(s)
end


function main()

    args = parse_commandline()
    input_file = args["input_tree"]
    output_file = args["output_tree"]

    tree = gettreefromnewick(input_file, FelNode)
    treestring = newick(tree)

    # Propagate information
    propagate_info(tree)

    # Write the modified tree to the output file
    treestring = newick(tree)

    file = open(output_file, "w")
    write(file, treestring)
    close(file)
end

# Run the script
main()
