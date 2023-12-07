using Pkg
Pkg.activate("environments/CodonMolecularEvolution")
#using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution
using ArgParse
using CSV
using DataFrames


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

    time_rows = []
    for tag_pos in pairwise_comparisons
        tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
        tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
        analysis_name = parsed_args["output_dir"] * ("/" * join([tag1, tag2], "v") * "/" * "results")

        analysis_name_no_prune = analysis_name * "/" * "no_prune_only_grid"
        analysis_name_prune_max = analysis_name * "/" * "prune_max_only_grid"
        analysis_name_prune_patrick = analysis_name * "/" * "prune_patrick_only_grid"

        analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
        cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)


        @time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, analysis_tags, analysis_tag_colors, exports=true, verbosity=1)
        @time code = MolecularEvolution.universal_code
        @time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s


        tree_no_prune = deepcopy(tree)
        tree_prune_1 = deepcopy(tree)
        tree_prune_2 = deepcopy(tree)
        tree_prune_3 = deepcopy(tree)
        tree_prune_4 = deepcopy(tree)
        tree_prune_5 = deepcopy(tree)

        time_difFUBAR = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 50s

        time_difFUBAR_prune_max = @elapsed con_lik_matrix_prune_1, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_1, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s

        time_difFUBAR_prune_patrick = @elapsed con_lik_matrix_prune_2, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s

        time_difFUBAR_prune_patrick_max = @elapsed con_lik_matrix_prune_3, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_3(tree_prune_3, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s

        time_difFUBAR_prune_patrick_max_child = @elapsed con_lik_matrix_prune_4, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_4(tree_prune_4, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s

        time_difFUBAR_prune_final = @elapsed con_lik_matrix_prune_5, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_final(tree_prune_5, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s


        times = [parsed_args["output_dir"], tag_pos, time_difFUBAR, time_difFUBAR_prune_max, time_difFUBAR_prune_patrick, time_difFUBAR_prune_patrick_max, time_difFUBAR_prune_patrick_max_child, time_difFUBAR_prune_final]

        push!(time_rows, times)

        #time_difFUBAR = @elapsed difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name_no_prune;
        #    pos_thresh=parsed_args["pos_thresh"], iters=parsed_args["iters"], verbosity=parsed_args["verbosity"],
        #    exports=parsed_args["exports"], code=MolecularEvolution.universal_code)
        #
        #time_difFUBAR_prune_max = @elapsed difFUBAR_prune_max(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name_prune_max;
        #    pos_thresh=parsed_args["pos_thresh"], iters=parsed_args["iters"], verbosity=parsed_args["verbosity"],
        #    exports=parsed_args["exports"], code=MolecularEvolution.universal_code)
        #
        #time_difFUBAR_prune_patrick = @elapsed difFUBAR_prune_patrick(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name_prune_patrick;
        #    pos_thresh=parsed_args["pos_thresh"], iters=parsed_args["iters"], verbosity=parsed_args["verbosity"],
        #    exports=parsed_args["exports"], code=MolecularEvolution.universal_code)

        #times = [parsed_args["output_dir"], tag_pos, time_difFUBAR, time_difFUBAR_prune_max, time_difFUBAR_prune_patrick]
        #push!(time_rows, times)
    end

    print(time_rows)
    data = time_rows
    # Convert to DataFrame
    df = DataFrame(output_dir=[data[1][1]],
        tag_pos=[data[1][2]],
        time_difFUBAR=[data[1][3]],
        time_difFUBAR_prune_max=[data[1][4]],
        time_difFUBAR_prune_patrick=[data[1][5]],
        time_difFUBAR_prune_patrick_max=[data[1][6]],
        time_difFUBAR_prune_patrick_max_child=[data[1][7]],
        time_difFUBAR_prune_final=[data[1][8]])

    # Print DataFrame

    # output times
    output_filename = parsed_args["output_dir"] * "/time_results.csv"

    # make path if it does not exist
    if !isdir(parsed_args["output_dir"])
        mkdir(parsed_args["output_dir"])
    end

    CSV.write(output_filename, df)

    #open(output_filename, "w") do f
    #    CSV.write(f, time_rows)
    #end
end

main()