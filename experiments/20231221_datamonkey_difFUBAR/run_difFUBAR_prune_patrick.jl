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
    nexus_file = args["nexus_file"]
    exports = args["exports"]
    verbosity = args["verbosity"]
    iters = args["iters"]
    pos_thresh = args["pos_thresh"]
    output_csv = args["output"]

    if !isdir(output_csv)
        mkpath(output_csv)
    end

    seqnames, seqs = read_fasta(fasta_file)
    treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nexus_file)

    df, results = difFUBAR_prune_patrick(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)

    # Save DataFrame to CSV
    CSV.write(analysis_name * "df_results.csv", df)
end

main()