using Pkg

Pkg.activate("environments/CodonMolecularEvolution")
Pkg.status()
cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution

#cd("/home/ptruong/git/computationalPhylogenetics/contrastFEL_data/") 
cd("/home/ptruong/git/computationalPhylogenetics/")


exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "output/1v2/test_sim_data"


cd("/home/ptruong/git/computationalPhylogenetics/")
file_name = "contrastFEL_data/omnibus-multi/datasets/omnibus-multi/sims.24.settings.replicate.1"
seqnames, seqs = import_hyphy_simulated_FASTA(file_name)
tree_file = "contrastFEL_data/omnibus-multi/datasets/omnibus-multi/sims.24.nwk"
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)

comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

for tag_pos in pairwise_comparisons
    analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
    cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)
    println("/" * join(tag_pos, "v") * "/" * analysis_name)

    difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, analysis_name;
        pos_thresh=pos_thresh, iters=iters, verbosity=verbosity,
        exports=exports, code=MolecularEvolution.universal_code)
end

result_string = "/" * join([1, 2], "v") * "/"



tag_pos = [1, 3]
analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)
difFUBAR(seqnames, seqs, cleaned_tree, analysis_tags, analysis_tag_colors, parsed_args["output_dir"];
    pos_thresh=pos_thresh, iters=iters, verbosity=verbosity,
    exports=exports, code=MolecularEvolution.universal_code)

# Same as difFUBAR
#@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, cleaned_tree, analysis_tags, analysis_tag_colors, exports=exports, verbosity=verbosity)
@time code = MolecularEvolution.universal_code
@time code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code) #19s
@time con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 73 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)



