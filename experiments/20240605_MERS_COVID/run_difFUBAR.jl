using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Revise
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240605_MERS_COVID/results")

analysis_name = "difFUBAR/results"
fasta_file = "result.fasta"
nexus_file = "result.newick"
newick_file = "result.newick"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95



seqnames, seqs = read_fasta(fasta_file)
treestring_group_labeled, treestring, group_tags, tags, tag_colors = import_grouped_label_tree(newick_file)

df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)
