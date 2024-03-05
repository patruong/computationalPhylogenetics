using MolecularEvolution #read_fasta
using CodonMolecularEvolution
cd("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline")

analysis_name = "results/results"
fasta_file = "data/combined_difFUBAR_input.fasta"
fasta_file = "results/nextclade_output/nextclade.aligned.fixed.fasta"
newick_file = "results/nextclade_output/nextclade_fixed_tagged.nwk"


# test fasta
cd("/home/patrick/git/computationalPhylogenetics/experiments/20231222_rebuild_parvoVP_on_difFUBAR")

fasta_file = "parvo_panleu_trans.fasta"
newick_file = "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped.nwk"
newick_file = "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95

seqnames, seqs = read_fasta(fasta_file)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)

df, results = difFUBAR_treesurgery_and_parallel(seqnames, seqs, treestring, original_tags, tag_colors, analysis_name, exports=exports, iters=iters, verbosity=verbosity)
