using MolecularEvolution #read_fasta
using CodonMolecularEvolution
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/results")
analysis_name = "difFUBAR_output/results"
fasta_file = "nextclade_output/nextclade.cds_translation.S.fasta.aligned.removed.tag.postprocessed.fasta"
tree_file = "nextclade_output/nextclade.tagged.internal.tag.nwk"



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


