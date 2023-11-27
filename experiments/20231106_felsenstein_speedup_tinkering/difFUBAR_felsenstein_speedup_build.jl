using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo



cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231106_felsenstein_speedup_tinkering")
analysis_name = "Ace2"
seqnames, seqs = read_fasta("data/Ace2_reallytiny.fasta");
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree("data/ACE2_reallytiny_tagged.nex")
#df, results = difFUBAR(seqnames, seqs, treestring, tags, tag_colors, analysis_name, exports=false, verbosity=0)

#difFUBAR(seqnames, seqs, treestring, tags, tag_colors, outpath; pos_thresh=0.95, iters=2500, verbosity=1, exports=true, code=MolecularEvolution.universal_code)
#analysis_name = outpath
@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code

seqnames
seqs

@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #6.9s

tree
GTRmat
F3x4_freqs
eq_freqs

#difFUBAR grid
con_lik_matrix, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 8.14 s

con_lik_matrix
log_con_lik_matrix
codon_param_vec
alphagrid
omegagrid
param_kinds

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
#Defines the grid used for inference.
function gridsetup(lb, ub, num_below_one, trin, tr)
    step = (trin(1.0) - trin(lb)) / num_below_one
    return tr.(trin(lb):step:trin(ub))
end

foreground_grid = 6
background_grid = 4
alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2
num_groups = length(tags)
[println(n.name) for n in getnodelist(tree)]
[tags for n in getnodelist(tree)]
is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups
tensor_dims = 1 + num_groups + is_background

function add_to_each_element(vec_of_vec, elems)
    return [vcat(v, [e]) for v in vec_of_vec for e in elems]
end

# check what is happening in the loop
codon_param_vec
omegagrid
a = [[1], [2], [3]]
b = [1, 2, 3]
a = add_to_each_element(a, b)
a = add_to_each_element(a, b)


codon_param_vec = [[a] for a in alphagrid]
param_kinds = ["Alpha"]

for g in 1:num_groups
    push!(param_kinds, "OmegaG$(g)")
    codon_param_vec = add_to_each_element(codon_param_vec, omegagrid)
end

param_kinds
codon_param_vec


@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
alloc_grid
theta

@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=false)


@assert size(df) == (19, 8)


alpha
beta

# rebuilding felsenstein
cached_model = CodonMolecularEvolution.MG94_cacher(code)


function N_Omegas_model_func(tags, omega_vec, alpha, nuc_mat, F3x4, code)
    models = [cached_model(alpha, alpha * o, nuc_mat, F3x4, genetic_code=code) for o in omega_vec]
    return n::FelNode -> [models[model_ind(n.name, tags)]]
end

omegas = codon_param_vec[1]


tagged_models = N_Omegas_model_func(tags, omegas, alpha, GTRmat, F3x4_freqs, code)


felsenstein!(tree)

img = tree_draw(tree)
img |> PDF("imgout_.pdf", 10cm, 10cm)



##### 
# Tabulate
########
using DataFrames
using CSV

verbosity = 1
pos_thresh = 0.95
exports = true

function collapse_counts(param_vec, count_vec; cases=nothing)
    if isnothing(cases)
        cases = sort(union(param_vec))
    end
    d = Dict(zip(cases, 1:length(cases)))
    storage = zeros(Int64, length(cases))
    for i in 1:length(count_vec)
        storage[d[param_vec[i]]] += count_vec[i]
    end
    return storage ./ sum(storage)
end

grid_size, num_sites = size(alloc_grid)

r(s) = round(s, digits=4)

detected_sites = Int64[]
group1_volumes = Vector{Float64}[]
group2_volumes = Vector{Float64}[]
alpha_volumes = Vector{Float64}[]
detections = Vector{Float64}[] #legacy name - now includes all 4 "relevant" site posteriors
param_means = Vector{Float64}[]

ω1 = [c[2] for c in codon_param_vec]
ω2 = [c[3] for c in codon_param_vec]
alphas = [c[1] for c in codon_param_vec]
ω1_greater_filt = ω1 .> ω2
ω2_greater_filt = ω2 .> ω1
ω1_pos_filt = ω1 .> 1.0
ω2_pos_filt = ω2 .> 1.0

verbosity > 0 && println("Step 5: Tabulating and plotting. Detected sites:")
for site in 1:num_sites
    ω1_greater_posterior = sum(alloc_grid[ω1_greater_filt, site]) / sum(alloc_grid[:, site])
    ω2_greater_posterior = sum(alloc_grid[ω2_greater_filt, site]) / sum(alloc_grid[:, site])
    ω1_pos_posterior = sum(alloc_grid[ω1_pos_filt, site]) / sum(alloc_grid[:, site])
    ω2_pos_posterior = sum(alloc_grid[ω2_pos_filt, site]) / sum(alloc_grid[:, site])
    detecs = [ω1_greater_posterior, ω2_greater_posterior, ω1_pos_posterior, ω2_pos_posterior]

    site_counts_ω1 = collapse_counts(ω1, alloc_grid[:, site], cases=omegagrid)
    site_counts_ω2 = collapse_counts(ω2, alloc_grid[:, site], cases=omegagrid)
    site_counts_alphas = collapse_counts(alphas, alloc_grid[:, site], cases=alphagrid)

    mean_alpha = sum(site_counts_alphas .* alphagrid)
    mean_ω1 = sum(site_counts_ω1 .* omegagrid)
    mean_ω2 = sum(site_counts_ω2 .* omegagrid)

    push!(detections, detecs)
    push!(param_means, [mean_alpha, mean_ω1, mean_ω2])
    push!(group1_volumes, site_counts_ω1)
    push!(group2_volumes, site_counts_ω2)
    push!(alpha_volumes, site_counts_alphas)

    if maximum(detecs) > pos_thresh
        verbosity > 0 && print("Site $(site) - ")
        verbosity > 0 && print("P(ω1 > ω2):", ω1_greater_posterior)
        verbosity > 0 && print("; P(ω2 > ω1):", ω2_greater_posterior)
        verbosity > 0 && print("; P(ω1 > 1):", ω1_pos_posterior)
        verbosity > 0 && println("; P(ω2 > 1):", ω2_pos_posterior)
        push!(detected_sites, site)
    end
end

#Exporting site data
df = DataFrame()
df[!, "Codon Sites"] = [1:(num_sites+1);]
df[!, "P(ω1 > ω2)"] = [d[1] for d in detections]
df[!, "P(ω2 > ω1)"] = [d[2] for d in detections]
df[!, "P(ω1 > 1)"] = [d[3] for d in detections]
df[!, "P(ω2 > 1)"] = [d[4] for d in detections]
df[!, "mean(α)"] = [d[1] for d in param_means]
df[!, "mean(ω1)"] = [d[2] for d in param_means]
df[!, "mean(ω2)"] = [d[3] for d in param_means]

verbosity > 0 && println("\nIf exports = true, writing results for all sites to CSV: " * analysis_name * "_posteriors.csv")
exports && CSV.write(analysis_name * "_posteriors.csv", df)

sites = [1:num_sites;]

#Select the sites that will get plotted, in case you want to customize this.
if isnothing(sites_to_plot)
    sites_to_plot = detected_sites
end


@time alloc_grid, theta = CodonMolecularEvolution.difFUBAR_sample(con_lik_matrix, 2500, verbosity=1)
@time df = CodonMolecularEvolution.difFUBAR_tabulate(analysis_name, 0.95, alloc_grid, codon_param_vec, alphagrid, omegagrid, tag_colors; verbosity=1, exports=true)



