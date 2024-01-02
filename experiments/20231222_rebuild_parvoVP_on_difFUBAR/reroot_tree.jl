using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse

#dir = "/home/patrick/git/computationalPhylogenetics/experiments/20231222_rebuild_parvoVP_on_difFUBAR/"
#cd(dir)

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--input", "-i"
        help = "Path to the input newick file"
        arg_type = String
        "--output_reroot", "-r"
        arg_type = String
        help = "Output newick file and tree .svg for rerooted tree"
        default = "output_reroot.nwk"
        "--output", "-o"
        arg_type = String
        help = "Output newick file and tree .svg for rerooted and reassigned branch tree"
        default = "output.nwk"
    end

    return parse_args(s)
end


function main()

    args = parse_commandline()

    input_file = args["input"]
    output_reroot = args["output_reroot"]
    output_reroot_reassign = args["output"]

    # Defaults 
    #newick_file = dir * "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped.nwk"
    #input = "root"
    #output_reroot = "reroot"
    #output_reroot_reassign = "reroot_regroup" 

    # placeholder variables

    analysis_name = "placeholder" # used just for convenience for difFUBAR_init
    exports = true
    verbosity = 1
    iters = 2500
    pos_thresh = 0.95
    output_csv = "test/"

    function write_tree(treestring, file_path="output.txt")
        file = open(file_path, "w")
        write(file, treestring)
        close(file)
        println("Reroot and re-assigned branch newick tree has been written to $file_path.")
    end

    println("WARNING: dist_above_child=reroot_node.branchlength might cause issue in recursive_reroot might cause a problem!")

    treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(input_file)
    tag_pos = [1, 2]
    tags, tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
    tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(input_file, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)

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
    tree = MolecularEvolution.recursive_reroot(reroot_node, dist_above_child=reroot_node.branchlength)
    #tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)
    treestring = MolecularEvolution.newick(tree)
    tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(output_reroot, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    #tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)
    write_tree(treestring, output_reroot)


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

    treestring = MolecularEvolution.newick(tree)
    tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(output_reroot_reassign, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
    #tree_draw(tree, draw_labels=true, canvas_height=200cm, canvas_width=200cm, font_size=24pt)
    write_tree(treestring, output_reroot_reassign)

end


main()