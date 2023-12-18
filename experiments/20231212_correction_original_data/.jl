using ArgParse
using CSV

function main()
    parser = ArgParser("difFUBAR baseline")

    @add_arg_or_flag parser "analysis_name" s"A description of your analysis" default = "output/"
    @add_arg_or_flag parser "fasta_file" s"Path to the FASTA file" required = true
    @add_arg_or_flag parser "nexus_file" s"Path to the Nexus file" required = true
    @add_arg_or_flag parser "exports" b"Whether to export results" default = true
    @add_arg_or_flag parser "output_csv" s"Path to the output CSV file" default = "output.csv"

    args = parse_args(parser)

    analysis_name = args["analysis_name"]
    fasta_file = args["fasta_file"]
    nexus_file = args["nexus_file"]
    exports = args["exports"]
    output_csv = args["output_csv"]

    seqnames, seqs = read_fasta(fasta_file)
    treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nexus_file)

    df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, verbosity=0)

    # Save DataFrame to CSV
    CSV.write(output_csv, df)
end

main()