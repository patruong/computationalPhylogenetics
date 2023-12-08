using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo
using ArgParse
using CSV
using DataFrames



# Parvo data
cd("/home/patrick/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time/data/")
file_dir = "/home/patrick/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time/data/" # original data
fasta = file_dir * "ParvoVP.fasta"
tree_file = file_dir * "ParvoVP.nex"
analysis_name = "ParvoVP"
seqnames, seqs = read_fasta(fasta);
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)

time_difFUBAR_init = @elapsed tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)

@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

MolecularEvolution.write_nexus(file_dir * "parvoVP/ParvoVP_molecularEvo_parsed.nex ", tree)
