
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution
using CodonMolecularEvolution
using Measures


cd("/home/patrick/git/computationalPhylogenetics/")

function getpuresubclades(node::FelNode, tags::Vector{String}, pure_subclades=FelNode[])
    # Get the index of the node's tag
    tag_ind_of_node = CodonMolecularEvolution.model_ind(node.name, tags)

    # If the node is a leaf, it's pure
    if isleafnode(node)
        return pure_subclades, true, tag_ind_of_node
    end

    children_are_pure = Vector{Bool}()
    children_tag_inds = Vector{Int64}()

    for child in node.children
        pure_subclades, child_is_pure, tag_ind = getpuresubclades(child, tags, pure_subclades)
        push!(children_are_pure, child_is_pure)
        push!(children_tag_inds, tag_ind)
    end

    # Get the index of the node's first child's tag
    tag_ind_of_first_child = first(children_tag_inds)

    # This is the case where the subclade starting at node is pure
    if all(children_are_pure) && all(x == tag_ind_of_first_child for x in children_tag_inds)
        if tag_ind_of_node != tag_ind_of_first_child
            # The purity is broken at this node
            push!(pure_subclades, node)
            return pure_subclades, false, tag_ind_of_node
        end
        # The purity is not broken at this node
        return pure_subclades, true, tag_ind_of_node
    end

    # This is the case where some child has mixed tags or the children are pure with regards to different tags
    for (child_is_pure, child) in zip(children_are_pure, node.children)
        if !child_is_pure || isleafnode(child)
            # We don't want to push leaves into pure_subclades
            continue
        end
        push!(pure_subclades, child)
    end
    return pure_subclades, false, tag_ind_of_node
end

#Calculates the ratio of nodes that are in a pure clade to total nodes in the tree (1st return value)
#Calculates the number of pure omega and background omega clades (2nd and 3rd return values)
function get_purity_info(tree, tags, num_groups)
    pure_subclades, _, _ = getpuresubclades(tree, tags)
    c = 0
    for x in pure_subclades
        #The root of the pure clade is not counted
        c += length(getnodelist(x)) - 1
    end
    num_omega_clades = count(x -> CodonMolecularEvolution.model_ind(x.children[1].name, tags) <= num_groups, pure_subclades)
    num_background_omega_clades = length(pure_subclades) - num_omega_clades
    return c / length(getnodelist(tree)), num_omega_clades, num_background_omega_clades
end

function convert_string_time_to_seconds(time_str)
    split_time = split(time_str, "m")
    minutes = parse(Int, split_time[1])
    seconds_str = split(split_time[2], "s")[1]
    seconds = parse(Float64, seconds_str)
    total_seconds = minutes * 60 + seconds
    return total_seconds
end

function parse_time(time_output_file)
    data = read(time_output_file, String)
    # Split the data into lines
    lines = split(data, "\n")

    # Initialize variables to store time values
    real_time = ""
    user_time = ""
    sys_time = ""

    # Iterate through the lines to find the time data
    for line in lines
        if occursin("real", line)
            real_time = split(line)[2]
        elseif occursin("user", line)
            user_time = split(line)[2]
        elseif occursin("sys", line)
            sys_time = split(line)[2]
        end
    end
    real_time = convert_string_time_to_seconds(real_time)
    user_time = convert_string_time_to_seconds(user_time)
    sys_time = convert_string_time_to_seconds(sys_time)
    return real_time, user_time, sys_time
end

function get_purity_from_newick(newick_file)
    treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)
    tree = gettreefromnewick(treestring, FelNode)
    num_groups = length(tag_colors)
    purity, _, _ = get_purity_info(tree, original_tags, num_groups)
    return purity
end

datasets_folder = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT_branchlength_1"]
datasets_name = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT"]
nwk_files = ["epidermal_leaf_trichomes_BRT.nwk", "hiv-1_envelope.nwk", "rubisco_C3_vs_C4.nwk", "HIV_RT_branchlength_1.nwk"]

contrastFEL_time = []
difFUBAR_time = []
purity = []

for i in 1:length(datasets_folder)
    dataset = datasets_folder[i]
    contrastFEL_time_file_location = "experiments/20240308_contrastFEL_empirical/resultsfile/$dataset/contrastfel/time_output.txt"
    difFUBAR_time_file_location = "experiments/20240308_contrastFEL_empirical/resultsfile/$dataset/treesurgery_and_parallel/time_output.txt"
    contrastFEL_real_time, contrastFEL_user_time, contrastFEL_sys_time = parse_time(contrastFEL_time_file_location)
    difFUBAR_real_time, difFUBAR_user_time, difFUBAR_sys_time = parse_time(difFUBAR_time_file_location)
    tree_purity = get_purity_from_newick("experiments/20240308_contrastFEL_empirical/data/$(nwk_files[i])")

    push!(contrastFEL_time, contrastFEL_real_time)
    push!(difFUBAR_time, difFUBAR_real_time)
    push!(purity, tree_purity)

end



df = DataFrame(
    Dataset=datasets_name,
    contrastFEL_time=contrastFEL_time,
    difFUBAR_time=difFUBAR_time,
    purity=purity
)

df[!, "speedup"] = (df[!, "contrastFEL_time"] ./ df[!, "difFUBAR_time"])


df[!, :contrastFEL_time] .= round.(df[!, :contrastFEL_time], digits=1)
df[!, :difFUBAR_time] .= round.(df[!, :difFUBAR_time], digits=1)
df[!, :purity] .= round.(df[!, :purity], digits=2)
df[!, :speedup] .= round.(df[!, :speedup], digits=1)

using Latexify

# Export DataFrame to LaTeX
latex_table = pretty_table(latexify(df), header=false)

df

output_path = "results/empirical_data/contrastfel_empirical/time_benchmark.csv"
CSV.write(output_path, res)


#########
# plot ##
#########

sort!(df, "contrastFEL_time")
scatter(df[!, "difFUBAR_time"], df[!, "contrastFEL_time"], xlabel="difFUBAR time", ylabel="contrastFEL time", label="", gridlinewidth=3)

output_path = "results/empirical_data/contrastfel_empirical/time_plot.png"
plot!(margin=2mm)

savefig(output_path)

sort!(df, "speedup")
scatter(df[!, "Dataset"], df[!, "speedup"], xlabel="Dataset", ylabel="Speedup", gridlinewidth=3,
    legend=false, series_annotations=text.(round.(df[!, "purity"], digits=2), :bottom))
plot!(margin=2mm)

output_path = "results/empirical_data/contrastfel_empirical/speedup_plot.png"
savefig(output_path)




