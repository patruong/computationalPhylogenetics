using MolecularEvolution #read_fasta
using CodonMolecularEvolution
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240301_replicate_ncbi_cleanup_to_build_tree_for_covid_data")
analysis_name = "results/difFUBAR_output/results"
fasta_file = "results/mafft_alignment_removed_tags.fasta"
tree_file = "results/iqtree_fixed_internal.nwk"


exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95

seqnames, seqs = read_fasta(fasta_file)
#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
treestring_group_labeled, treestring, group_tags, tags, tag_colors = import_grouped_label_tree(tree_file)

df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)


# This works, but why does nextclade version of this not work?

# Save DataFrame to CSV
#CSV.write(analysis_name * "df_results.csv", df)


