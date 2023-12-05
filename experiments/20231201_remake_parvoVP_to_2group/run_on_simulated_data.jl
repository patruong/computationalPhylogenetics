using Pkg
Pkg.activate("environments/CodonMolecularEvolution")
#using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
using MolecularEvolution, FASTX, CodonMolecularEvolution
using ArgParse
using CSV
using DataFrames

cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231201_remake_parvoVP_to_2group/")
#fasta_file = "data/sim_reference_test/sim.replicate.1.nex"
#tree_file = "data/sim_reference_test/sim_tree.nwk" # 
fasta_file = "data/sim_tree_fixed/sim.replicate.1.nex"
tree_file = "data/sim_tree_fixed/sim_tree_fixed.nwk" # fixed
output_dir = "output/"

#hyphy sim
#fasta_file = "data/hyphy_sim/sims.0.settings.replicate.1"
#tree_file = "data/hyphy_sim/sims.0.nwk" # 



#println()
seqnames, seqs = import_hyphy_simulated_FASTA(fasta_file)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(tree_file)

comparisons_idx = [i for i in 1:length(original_tags)]
pairwise_comparisons = [[i, j] for i in comparisons_idx for j in comparisons_idx if i < j]

# Printout for information
group_name = ["$i = $(original_tags[i])" for i in 1:length(original_tags)]
pairwise_comparisons_names = [[original_tags[i], original_tags[j]] for i in comparisons_idx for j in comparisons_idx if i < j]

# group_name printout
group_meta_data_file_name = "group_meta_data.txt"
group_meta_data_file = open(group_meta_data_file_name, "w")
for element in group_name
    println(group_meta_data_file, "Groups")
    println(group_meta_data_file, element)
end
println(group_meta_data_file, "")

# pairwise comparison printout
for element in pairwise_comparisons_names
    println(group_meta_data_file, "Comparisons")
    println(group_meta_data_file, element)
end
close(group_meta_data_file)

#time_rows = []
#for tag_pos in pairwise_comparisons

tag_pos = pairwise_comparisons[1]
tag1 = chop(original_tags[tag_pos[1]], head=1, tail=1)
tag2 = chop(original_tags[tag_pos[2]], head=1, tail=1)
analysis_name = output_dir * ("/" * join([tag1, tag2], "v") * "/" * "results")

analysis_name_no_prune = analysis_name * "/" * "no_prune_only_grid"
analysis_name_prune_max = analysis_name * "/" * "prune_max_only_grid"
analysis_name_prune_patrick = analysis_name * "/" * "prune_patrick_only_grid"

analysis_tags, analysis_tag_colors, remove_tags = select_analysis_tags_from_newick_tree(original_tags, tag_colors, tag_pos)
cleaned_tree = remove_tags_from_newick_tree(treestring, remove_tags)

@time tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, analysis_tags, analysis_tag_colors, exports=true, verbosity=1)
@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

tree_no_prune = deepcopy(tree)
tree_prune_1 = deepcopy(tree)
tree_prune_2 = deepcopy(tree)
tree_prune_3 = deepcopy(tree)
tree_prune_4 = deepcopy(tree)
tree_prune_final = deepcopy(tree)

time_difFUBAR = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 50s

time_difFUBAR_prune_max = @elapsed con_lik_matrix_prune_1, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_1, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s

time_difFUBAR_prune_patrick = @elapsed con_lik_matrix_prune_2, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s

time_difFUBAR_prune_patrick_max = @elapsed con_lik_matrix_prune_3, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_3(tree_prune_3, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s

time_difFUBAR_prune_patrick_max_child = @elapsed con_lik_matrix_prune_4, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_4(tree_prune_4, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s

time_difFUBAR_prune_final = @elapsed con_lik_matrix_prune_4, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_final(tree_prune_final, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s


[time_difFUBAR, time_difFUBAR_prune_max, time_difFUBAR_prune_patrick, time_difFUBAR_prune_patrick_max, time_difFUBAR_prune_patrick_max_child]
# return subtree tops
prune_1_subtree_tops = CodonMolecularEvolution.difFUBAR_grid_pruned_1_return_subtree_tops(tree_prune_1, tags, GTRmat, F3x4_freqs, code; verbosity=1, foreground_grid=6, background_grid=4)

prune_3_subtree_tops = CodonMolecularEvolution.difFUBAR_grid_pruned_3_return_subtree_tops(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4)


# compare speed over many data set
time_difFUBAR_prune_max
time_difFUBAR_prune_patrick
time_difFUBAR_prune_patrick_max
length(prune_1_subtree_tops)
length(prune_3_subtree_tops)

# seems like saving things in the dict just makes the overall run-time slower.... :/

save_tree_svg(prune_1_subtree_tops[1], "prune_1_subtree_1")
save_tree_svg(prune_1_subtree_tops[2], "prune_1_subtree_2")
save_tree_svg(prune_3_subtree_tops[1], "prune_3_subtree_1")
save_tree_svg(prune_3_subtree_tops[2], "prune_3_subtree_2")

con_lik_matrix_no_prune == con_lik_matrix_prune_1
con_lik_matrix_no_prune == con_lik_matrix_prune_2
con_lik_matrix_prune_1 == con_lik_matrix_prune_2
con_lik_matrix_prune_1 == con_lik_matrix_prune_3
con_lik_matrix_prune_1 == con_lik_matrix_prune_4
con_lik_matrix_prune_3
con_lik_matrix_no_prune == con_lik_matrix_prune_3

raw_times = [time_difFUBAR, time_difFUBAR_prune_max, time_difFUBAR_prune_patrick]


times = [parsed_args["output_dir"], tag_pos, time_difFUBAR, time_difFUBAR_prune_max, time_difFUBAR_prune_patrick]
push!(time_rows, times)


### save pruned tree

function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end
function save_tree_svg(tree, analysis_name)
    strip_tags_from_name = CodonMolecularEvolution.generate_tag_stripper(tags)       #Replace with Phylo.jl based plot?

    color_dict = Dict(zip(getnodelist(tree), [tag_colors[model_ind(n.name, tags)] for n in getnodelist(tree)]))
    label_dict = Dict(zip(getnodelist(tree), [strip_tags_from_name(n.name) for n in getnodelist(tree)]))
    img = tree_draw(tree, canvas_height=(3 + length(getleaflist(tree)) / 5)cm,
        draw_labels=true, dot_color_dict=color_dict,
        line_color_dict=color_dict, line_width=0.3, min_dot_size=0.01,
        nodelabel_dict=label_dict)
    img |> SVG(analysis_name * "_tagged_input_tree.svg", 15cm, (3 + length(getleaflist(tree)) / 5)cm)
end

save_tree_svg(tree_no_prune, "no_prune")
save_tree_svg(tree_prune_1, "prune_max")
save_tree_svg(tree_prune_2, "prune_patrick")

# Codon param vec alpha beta cached_model = MG94_cacher(code)

foreground_grid = 6
background_grid = 4
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

alphagrid
omegagrid
omegagrid
background_omega_grid




### test is_background
# set is_background to flase in the no_prune method
cached_model = CodonMolecularEvolution.MG94_cacher(code)

foreground_grid = 6
background_grid = 4
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
tr(x) = 10^x - 0.05
trinv(x) = log10(x + 0.05)
alphagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
omegagrid = gridsetup(0.01, 13.0, foreground_grid, trinv, tr)
background_omega_grid = gridsetup(0.05, 6.0, background_grid, trinv, tr) #Much coarser, because this isn't a target of inference
length(background_omega_grid) * length(alphagrid) * length(omegagrid)^2

num_groups = length(tags)
is_background = maximum([model_ind(n.name, tags) for n in getnodelist(tree)]) > num_groups

num_groups
maximum([model_ind(n.name, tags) for n in getnodelist(tree)])
[model_ind(n.name, tags) for n in getnodelist(tree)]


# Change is_background to omit root

list = [model_ind(n.name, tags) for n in getnodelist(tree)]

is_background_list = []
n = getnodelist(tree)
for n in getnodelist(tree)
    if !isroot(n)
        push!(is_background_list, model_ind(n.name, tags))
    end
end



