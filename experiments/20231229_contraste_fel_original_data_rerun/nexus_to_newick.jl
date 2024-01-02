using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse
using CSV
using DataFrames
using MolecularEvolution, Plots, Phylo


cd("/home/patrick/git/computationalPhylogenetics/experiments/20231212_correction_original_data/")


fasta = "ParvoVP.fasta" #parsed_args["FASTAfile"]
tree_file = "ParvoVP.nex" #parsed_args["treefile"]
analysis_name = "output" #parsed_args["output_dir"] * "/" * "ParvoVP"
seqnames, seqs = read_fasta(fasta)
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
time_difFUBAR_init = @elapsed tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code

# write Newick tree for difFUBAR run
nwk_string = MolecularEvolution.newick(tree)
file = open("ParvoVP.nwk", "w")
write(file, nwk_string)
close(file)





