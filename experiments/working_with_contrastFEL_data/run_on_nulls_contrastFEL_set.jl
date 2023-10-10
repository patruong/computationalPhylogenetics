using Pkg
Pkg.activate("environments/CodonMolecularEvolution")
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution

# setwd
cd("/home/ptruong/git/computationalPhylogenetics/")



# Parameters
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "experiments/working_with_contrastFEL_data/output/"
data_dir = "contrastFEL_data/null/datasets/null/64/random-64-0.1-0.5/"
file_name = data_dir * "replicate.replicate.1"
tree_file = data_dir * "tree.nwk"

seqnames, seqs = import_hyphy_simulated_FASTA(file_name)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)


difFUBAR(seqnames, seqs, treestring, original_tags, tag_colors, analysis_name;
    pos_thresh=pos_thresh, iters=iters, verbosity=verbosity, exports=exports, code=MolecularEvolution.universal_code)



