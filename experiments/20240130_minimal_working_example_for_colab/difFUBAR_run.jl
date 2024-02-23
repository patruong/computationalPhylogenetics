using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using ArgParse
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240130_minimal_working_example_for_colab")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--analysis_name", "-a"
        help = "Analysis output folder and prefix"
        arg_type = String
        "--fasta_file", "-f"
        help = "Fasta file"
        arg_type = String
        "--nexus_file", "-n"
        arg_type = String
        help = "Nexus tree file"
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

    #analysis_name = "results/results"
    #fasta_file = "Ace2_tiny_test.fasta"
    #nexus_file = "Ace2_no_background.nex"
    #exports = true
    #verbosity = 1
    #iters = 2500
    #pos_thresh = 0.95

    analysis_name = parsed_args["analysis_name"]
    fasta_file = parsed_args["fasta_file"]
    nexus_file = parsed_args["nexus_file"]
    exports = parsed_args["exports"]
    verbosity = parsed_args["verbosity"]
    iters = parsed_args["iters"]
    pos_thresh = parsed_args["pos_thresh"]


    seqnames, seqs = read_fasta(fasta_file)
    treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nexus_file)

    df, results = difFUBAR_prune_patrick_max_child(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)
end

main()
# julia difFUBAR_run.jl --analysis_name results/results --fasta_file Ace2_tiny_test.fasta --nexus_file Ace2_no_background.nex --pos_thresh 0.95 --verbosity 1 --iters 2500 --exports true