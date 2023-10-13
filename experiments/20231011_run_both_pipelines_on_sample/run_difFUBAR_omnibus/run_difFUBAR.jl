using Pkg
Pkg.activate("environments/CodonMolecularEvolution")
#using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution
using ArgParse

# setwd
#cd("/home/ptruong/git/computationalPhylogenetics/")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--FASTAfile", "-f"
        help = "A HyPhy simulation file"
        arg_type = String
        "--treefile", "-t"
        help = "Newick tree file"
        arg_type = String
        "--output_dir", "-o"
        arg_type = String
        help = "Analysis output directory"
        "--pos_thresh", "-p"
        arg_type = Float64
        default = 0.95
        help = "Posterior threshold"
        "--verbosity"
        arg_type = Int
        default = 1
        help = "Verbosity"
        "--iters"
        arg_type = Int
        default = 2500
        help = "Iters"
        "--exports"
        arg_type = Bool
        default = true
        help = "Export csv files"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg, val) in parsed_args
        println("  $arg  =>  $val")
    end
    println("working dir:" * string(pwd()))

    #println()
    seqnames, seqs = import_hyphy_simulated_FASTA(parsed_args["FASTAfile"])
    treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(parsed_args["treefile"])

    comparisons_idx = [i for i in 1:length(original_tags)]
    pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

    # Printout for information
    group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
    pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

    # group_name printout
    group_meta_data_file_name = "group_meta_data.txt"
    group_meta_data_file = open(group_meta_data_file_name, "w")
    for element in group_name
        println(group_meta_data_file, "Groups")
        println(group_meta_data_file, element)
    end
    println(group_meta_data_file, "")

    # pairwise comparison printout
    for element in pairwise_comparisons_names
        println(group_meta_data_file, "Comparisons")
        println(group_meta_data_file, element)
    end
    close(group_meta_data_file)

    for tag_pos in pairwise_comparisons
        tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
        tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
        analysis_name = parsed_args["output_dir"] * ("/" * join([tag1, tag2], "v") * "/" * "results")
        analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
        cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)

        difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name;
            pos_thresh=parsed_args["pos_thresh"], iters=parsed_args["iters"], verbosity=parsed_args["verbosity"],
            exports=parsed_args["exports"], code=MolecularEvolution.universal_code)
    end
end

main()

#fasta = "../../contrastFEL_data/omnibus-multi/datasets/omnibus-multi/sims.9.settings.replicate.1"
#tree = "../../contrastFEL_data//omnibus-multi/datasets/omnibus-multi/sims.9.nwk"
#seqnames, seqs = import_hyphy_simulated_FASTA(fasta)
#treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree)
#Parameters
#
#comparisons_idx = [i for i in 1:length(original_tags)]
#pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]
#
## Printout for information
#group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
#pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]
#
## group_name printout
#group_meta_data_file_name = "group_meta_data.txt"
#group_meta_data_file = open(group_meta_data_file_name, "w")
#for element in group_name
#    println(group_meta_data_file, "Groups")
#    println(group_meta_data_file, element)
#end
#println(group_meta_data_file, "")
#
## pairwise comparison printout
#for element in pairwise_comparisons_names
#    println(group_meta_data_file, "Comparisons")
#    println(group_meta_data_file, element)
#end
#close(group_meta_data_file)
#
#exports = true
#verbosity = 1
#iters = 2500
#pos_thresh = 0.95
#analysis_name = "experiments/working_with_contrastFEL_data/output/"
#data_dir = "contrastFEL_data/null/datasets/null/64/random-64-0.1-0.5/"
#file_name = data_dir * "replicate.replicate.1"
#tree_file = data_dir * "tree.nwk"
#
#
#
#seqnames, seqs = import_hyphy_simulated_FASTA(file_name)
#treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)
#difFUBAR(seqnames, seqs, treestring, original_tags, tag_colors, analysis_name;
#    pos_thresh=pos_thresh, iters=iters, verbosity=verbosity, exports=exports, code=MolecularEvolution.universal_code)
#
#