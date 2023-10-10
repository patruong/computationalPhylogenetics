using Pkg

Pkg.activate("environments/CodonMolecularEvolution")
Pkg.status()
cd("/home/ptruong/git/CodonMolecularEvolution.jl/test")
pwd()
using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution

cd("/home/ptruong/git/computationalPhylogenetics/workflow_pipeline/") # how do we fix this to snakemake?

exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
analysis_name = "output/results/difFUBAR"
seqnames, seqs = read_fasta("output/data/flu_sim_MG94_F61.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("output/data/flu_sim.nex")

seqnames, seqs = read_fasta("output/data/test_data/Ace2_tiny_tagged.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("output/data/test_data/Ace2_no_background.nex")

#@time df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name);

@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=exports, verbosity=verbosity)
@time code = MolecularEvolution.universal_code
@time code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code) #19s
@time con_lik_matrix, _, codon_param_vec, alphagrid, omegagrid, _ = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 73 s
@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, iters, verbosity=verbosity)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, pos_thresh, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=verbosity, exports=exports)
