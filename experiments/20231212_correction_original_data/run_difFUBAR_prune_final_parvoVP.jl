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

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--fasta_file", "-f"
        help = "Path to the FASTA file"
        arg_type = String
        "--nexus_file", "-n"
        help = "Path the the NEXUS file"
        arg_type = String
        "--output", "-o"
        arg_type = String
        help = "Analysis output directory"
        default = "output"
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
    args = parse_commandline()

    analysis_name = args["output"]
    fasta_file = args["fasta_file"]
    newick_file = args["nexus_file"]
    exports = args["exports"]
    verbosity = args["verbosity"]
    iters = args["iters"]
    pos_thresh = args["pos_thresh"]
    output_csv = args["output"]

    if !isdir(output_csv)
        mkpath(output_csv)
    end

    seqnames, seqs = read_fasta(fasta_file)
    treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)

    comparisons_idx = [i for i in 1:length(original_tags)]
    pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

    # Printout for information
    group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
    pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

    # group_name printout
    group_meta_data_file_name = output_csv * "group_meta_data.txt"
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


    tag_pos = [1, 2]
    tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
    tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
    tags, tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)

    df, results = difFUBAR_prune_final(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)

    # Save DataFrame to CSV
    CSV.write(analysis_name * "df_results.csv", df)
end

main()