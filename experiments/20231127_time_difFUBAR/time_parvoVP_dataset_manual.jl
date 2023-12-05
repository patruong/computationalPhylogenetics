using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo
using ArgParse
using CSV
using DataFrames



# Parvo data
cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231127_time_difFUBAR")
file_dir = "/home/ptruong/git/computationalPhylogenetics/experiments/20231127_time_difFUBAR/data/" # original data
fasta = file_dir * "ParvoVP.fasta"
tree_file = file_dir * "ParvoVP.nex"
analysis_name = "ParvoVP"
seqnames, seqs = read_fasta(fasta);
treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)


time_difFUBAR_init = @elapsed tree, tags, tag_colors, analysis_name = CodonMolecularEvolution.difFUBAR_init(analysis_name, treestring, tags, tag_colors, exports=true, verbosity=1)

@time code = MolecularEvolution.universal_code
@time tree, alpha, beta, GTRmat, F3x4_freqs, eq_freqs = CodonMolecularEvolution.difFUBAR_global_fit(seqnames, seqs, tree, CodonMolecularEvolution.generate_tag_stripper(tags), code, verbosity=1) #60s

tree_no_prune = deepcopy(tree)
tree_prune_1 = deepcopy(tree)
tree_prune_2 = deepcopy(tree)


times = []

time_difFUBAR_grid = @elapsed con_lik_matrix_no_prune, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid(tree_no_prune, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 50s
push!(times, time_difFUBAR_grid)

time_difFUBAR_grid_prune_max = @elapsed con_lik_matrix_prune_1, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_1(tree_prune_1, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 28s
push!(times, time_difFUBAR_grid_prune_max)

time_difFUBAR_grid_prune_patrick = @elapsed con_lik_matrix_prune_2, log_con_lik_matrix, codon_param_vec, alphagrid, omegagrid, param_kinds = CodonMolecularEvolution.difFUBAR_grid_pruned_2(tree_prune_2, tags, GTRmat, F3x4_freqs, code, verbosity=1, foreground_grid=6, background_grid=4) # 21s
push!(times, time_difFUBAR_grid_prune_patrick)

data = times
# Convert to DataFrame
df = DataFrame(output_dir=[data[1][1]],
    time_difFUBAR_grid=[data[1][2]],
    time_difFUBAR_grid_prune_max=[data[1][3]],
    time_difFUBAR_grid_prune_patrick=[data[1][4]])

# output times
output_filename = parsed_args["output_dir"] * "/time_results.csv"

CSV.write(output_filename, df)

################
# output trees #
################

function model_ind(str::String, tags::Vector{String})
    ind = length(tags) + 1
    for (i, t) in enumerate(tags)
        if occursin(t, str)
            ind = i
        end
    end
    return ind
end

strip_tags_from_name = CodonMolecularEvolution.generate_tag_stripper(tags)

function print_tree(tree, name)
    #Replace with Phylo.jl based plot?
    color_dict = Dict(zip(getnodelist(tree), [tag_colors[model_ind(n.name, tags)] for n in getnodelist(tree)]))
    label_dict = Dict(zip(getnodelist(tree), [strip_tags_from_name(n.name) for n in getnodelist(tree)]))
    img = tree_draw(tree, canvas_height=(3 + length(getleaflist(tree)) / 5)cm,
        draw_labels=true, dot_color_dict=color_dict,
        line_color_dict=color_dict, line_width=0.3, min_dot_size=0.01,
        nodelabel_dict=label_dict)
    img |> SVG(analysis_name * name, 15cm, (3 + length(getleaflist(tree)) / 5)cm)
end

print_tree(tree_prune_1, "_prune_1.svg")
print_tree(tree_prune_2, "_prune_2.svg")
