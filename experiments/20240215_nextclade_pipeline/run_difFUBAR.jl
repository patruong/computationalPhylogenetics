using MolecularEvolution #read_fasta
using CodonMolecularEvolution
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline")

analysis_name = "results/results"
fasta_file = "results/nextclade_output/nextclade.aligned.fasta"
newick_file = "results/nextclade_output/nextclade_fixed.nwk"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95

seqnames, seqs = read_fasta(fasta_file)
# Need to change and first label the newick file
# Need to find a nextclade newick reader

#tree = read_newick_tree(newick_file)
#treestring = read(newick_file, String)
#tree = gettreefromnewick(treestring, FelNode)
#treestring = newick(tree)

# why does tags not work?
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)

df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, original_tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)

# Save DataFrame to CSV
#CSV.write(analysis_name * "df_results.csv", df)
