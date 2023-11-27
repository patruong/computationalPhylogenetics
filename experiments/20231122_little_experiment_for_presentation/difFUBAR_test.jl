using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment




cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231122_little_experiment_for_presentation")
analysis_name = "nobackground/Ace2"
seqnames, seqs = read_fasta("Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("Ace2_no_background.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s
con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=true)

alloc_grid

@assert size(df) == (19, 8)



