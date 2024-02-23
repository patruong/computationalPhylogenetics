using MolecularEvolution #read_fasta
using CodonMolecularEvolution
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240130_minimal_working_example_for_colab")

analysis_name = "results/results"
fasta_file = "Ace2_tiny_test.fasta"
nexus_file = "Ace2_no_background.nex"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95

seqnames, seqs = read_fasta(fasta_file)
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(nexus_file)

df, results = difFUBAR_prune_patrick_max_child(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)

# Save DataFrame to CSV
#CSV.write(analysis_name * "df_results.csv", df)
