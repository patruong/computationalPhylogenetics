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

img = tree_draw(tree)
img |> PDF("imgout.pdf", 10cm, 10cm)

@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s
@time con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s

##################
# difFUBAR_grid ##
##################

verbosity = 1
foreground_grid = 6
background_grid = 4

cached_model = CodonMolecularEvolution.MG94_cacher(code)

#This is the function that assigns models to branches
#Sometimes there will be the same number of tags as omegas, but sometimes there will be one more omega.
#Depending on whether there are any background branches (UNTESTED WITH NO BACKGROUND BRANCHES)
function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
    models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
    return n::FelNode -> [models[model_ind(n.name, tags)]]
end

#Defines the grid used for inference.
function gridsetup(lb, ub, num_below_one, trin, tr)
    step = (trin(1.0) - trin(lb)) / num_below_one
    return tr.(trin(lb):step:trin(ub))
end

function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end

tr(x) = 10^x - 0.05
trinv(x) = log10(x + 0.05)
alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

num_groups = length(tags)
is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
tensor_dims = 1 + num_groups + is_background

function add_to_each_element(vec_of_vec, elems)
    return [vcat(v, [e]) for v in vec_of_vec for e in elems]
end

codon_param_vec = [[a] for a in alphagrid]
param_kinds = ["Alpha"]
for g in 1:num_groups
    push!(param_kinds, "OmegaG$(g)")
    codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
end
if is_background
    push!(param_kinds, "OmegaBackground")
    codon_param_vec = add_to_each_element(codon_param_vec, background_omega_grid)
end
codon_param_vec

num_sites = tree.message[1].sites
l = length(codon_param_vec)
log_con_lik_matrix = zeros(l, num_sites)

verbosity > 0 && println("Step 3: Calculating grid of $(length(codon_param_vec))-by-$(tree.message[1].sites) conditional likelihood values (the slowest step). Currently on:")


for (row_ind, cp) in enumerate(codon_param_vec)
    alpha = cp[1]
    omegas = cp[2:end]
    tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)

    #print("felsenstein")
    felsenstein!(tree, tagged_models)
    #This combine!() is needed because the current site_LLs function applies to a partition
    #And after a felsenstein pass, you don't have the eq freqs factored in.
    #We could make a version of log_likelihood() that returns the partitions instead of just the sum
    #print("combine")
    combine!.(tree.message, tree.parent_message)
    #println(length(tree.message[1]))
    #print("log_con_like_matrix")
    log_con_lik_matrix[row_ind, :] .= MolecularEvolution.site_LLs(tree.message[1]) #Check that these grab the scaling constants as well!
    verbosity > 0 && if mod(row_ind, 500) == 1
        print(round(100 * row_ind / length(codon_param_vec)), "% ")
        flush(stdout)
    end
end

verbosity > 0 && println()

con_lik_matrix = zeros(size(log_con_lik_matrix))
site_scalers = maximum(log_con_lik_matrix, dims=1)
for i in 1:num_sites
    con_lik_matrix[:, i] .= exp.(log_con_lik_matrix[:, i] .- site_scalers[i])
end






############
## Comment 
# Understand how the grid is built... we only need to compute one axis #
# alpha (12) x omega1 (12) x omega2 (12) x omega_background (7)
# 12*12*12*7 = 12096
# skipping the background we get
# alpha (12) x omega1 (12) x omega2 (12)
# 12*12*12 = 1728 i.e. the length of codon_param_vec
# so this is what we don't need to compute everything for... 
# 
# Look at the structure as we build up codon_param_vec

