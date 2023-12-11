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

    fasta = parsed_args["FASTAfile"]
    tree_file = parsed_args["treefile"]
    analysis_name = parsed_args["output_dir"] * "/" * "ParvoVP"
    seqnames, seqs = read_fasta(fasta)
    treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
    time_difFUBAR_init = @elapsed tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
    @time code = MolecularEvolution.universal_code
    @time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

    tree_no_prune = deepcopy(tree)
    tree_prune_1 = deepcopy(tree)
    tree_prune_2 = deepcopy(tree)
    tree_prune_3 = deepcopy(tree)
    tree_prune_4 = deepcopy(tree)
    tree_prune_5 = deepcopy(tree)

    times = []

    time_difFUBAR_grid = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 50s
    push!(times, time_difFUBAR_grid)

    time_difFUBAR_grid_prune_max = @elapsed con_lik_matrix_prune_1, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_1, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s
    push!(times, time_difFUBAR_grid_prune_max)

    time_difFUBAR_grid_prune_patrick = @elapsed con_lik_matrix_prune_2, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
    push!(times, time_difFUBAR_grid_prune_patrick)

    time_difFUBAR_grid_prune_patrick_max = @elapsed con_lik_matrix_prune_3, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_3, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s
    push!(times, time_difFUBAR_grid_prune_patrick_max)

    time_difFUBAR_grid_prune_patrick_max_child = @elapsed con_lik_matrix_prune_4, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_4, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
    push!(times, time_difFUBAR_grid_prune_patrick_max_child)

    time_difFUBAR_grid_prune_final = @elapsed con_lik_matrix_prune_5, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_5, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
    push!(times, time_difFUBAR_grid_prune_final)


    data = times
    # Convert to DataFrame
    df = DataFrame(time_difFUBAR_grid=[data[1]],
        time_difFUBAR_grid_prune_max=[data[2]],
        time_difFUBAR_grid_prune_patrick=[data[3]],
        time_difFUBAR_grid_prune_patrick_max=[data[4]],
        time_difFUBAR_grid_prune_patrick_max_child=[data[5]],
        time_difFUBAR_grid_prune_final=[data[6]])

    # output times
    output_filename = parsed_args["output_dir"] * "/" * "ParvoVP" * "/time_results_rerun_all_2.csv"

    CSV.write(output_filename, df)

    # output equal conditional likelihood matrices
        
    matrices = [con_lik_matrix_no_prune, con_lik_matrix_prune_1, con_lik_matrix_prune_2, con_lik_matrix_prune_3, con_lik_matrix_prune_4, con_lik_matrix_prune_5]
    matrix_names = ["con_lik_matrix_no_prune", "con_lik_matrix_prune_1", "con_lik_matrix_prune_2", "con_lik_matrix_prune_3", "con_lik_matrix_prune_4", "con_lik_matrix_prune_5"]

    n = length(matrices)
    results = Matrix{Bool}(undef, n, n) # create an empty boolean matrix of size n x n
    for i in 1:n
        for j in 1:n
            results[i, j] = all(matrices[i] == matrices[j]) # compare the i-th and j-th matrices and store the result
        end
    end

    df = DataFrame(Tables.table(results))
    rename!(df, matrix_names)
    df[!, "row"] = matrix_names

    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/equality_matrix.csv", df, writeheader=true)

    # output conditional likelihod matrices
    df = DataFrame(Tables.table(con_lik_matrix_no_prune))
    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/con_lik_matrix_no_prune.csv", df, writeheader=false)
    df = DataFrame(Tables.table(con_lik_matrix_prune_1))
    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/con_lik_matrix_prune_1.csv", df, writeheader=false)
    df = DataFrame(Tables.table(con_lik_matrix_prune_2))
    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/con_lik_matrix_prune_2.csv", df, writeheader=false)
    df = DataFrame(Tables.table(con_lik_matrix_prune_3))
    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/con_lik_matrix_prune_3.csv", df, writeheader=false)
    df = DataFrame(Tables.table(con_lik_matrix_prune_4))
    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/con_lik_matrix_prune_4.csv", df, writeheader=false)
    df = DataFrame(Tables.table(con_lik_matrix_prune_5))
    CSV.write(parsed_args["output_dir"] * "/" * "ParvoVP" * "/con_lik_matrix_prune_5.csv", df, writeheader=false)


end


main()

