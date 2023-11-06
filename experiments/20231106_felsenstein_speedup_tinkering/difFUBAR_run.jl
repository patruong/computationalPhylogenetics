using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo



cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
analysis_name = "Ace2"
seqnames, seqs = read_fasta("data/Ace2_tiny_test.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/Ace2_no_background.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=false, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s
@time con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=false)



