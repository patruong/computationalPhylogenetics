using ArgParse
using MolecularEvolution
using CodonMolecularEvolution
using Revise
using Dates

function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--analysis_name"
        help = "Name of the analysis"
        default = "difFUBAR/results"
        arg_type = String

        "--fasta_file"
        help = "Path to the fasta file"
        default = "result.fasta"
        arg_type = String

        "--newick_file"
        help = "Path to the newick file"
        default = "result.newick"
        arg_type = String

        "--verbosity"
        help = "Verbosity level"
        default = 1
        arg_type = Int64

        "--iters"
        help = "Number of iterations"
        default = 2500
        arg_type = Int64


        "--pos_thresh"
        help = "Positive selection threshold"
        default = 0.95
        arg_type = Float64
    end

    parsed_args = parse_args(s)

    analysis_name = parsed_args["analysis_name"]
    fasta_file = parsed_args["fasta_file"]
    newick_file = parsed_args["newick_file"]
    verbosity = parsed_args["verbosity"]
    iters = parsed_args["iters"]
    pos_thresh = parsed_args["pos_thresh"]

    #analysis_name = "test"
    #fasta_file = "results_aves_mammal_150/result.fasta"
    #newick_file = "results_aves_mammal_150/result.newick"
    #verbosity = 1
    #iters = 2500
    #typeof(iters)
    #pos_thresh = 0.95
    #typeof(pos_thresh)
    #cd("/home/patrick/git/computationalPhylogenetics/experiments/20240605_MERS_COVID")
    seqnames, seqs = read_fasta(fasta_file)
    treestring_group_labeled, treestring, group_tags, tags, tag_colors = import_grouped_labelled_newick_tree(newick_file)

    start_time = now()
    df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name; pos_thresh=pos_thresh, iters=iters, verbosity=verbosity, exports=true, code=MolecularEvolution.universal_code)

    for i in seqs
        println(i)
    end


    end_time = now()
    runtime = end_time - start_time

    runtime_file = joinpath(analysis_name, "runtime.txt")
    open(runtime_file, "w") do io
        println(io, "Runtime of $analysis_name: $runtime seconds")
        println(io, "Analysis completed.")
    end
end

main()
