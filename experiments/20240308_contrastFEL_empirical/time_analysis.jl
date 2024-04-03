
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution
using CodonMolecularEvolution
using Compose
cd("/home/patrick/git/computationalPhylogenetics/")

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

datasets_folder = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT_branchlength_1", "rubisco_rerooted_retagged"]

datasets_name = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT", "rubisco_reroot"]
nwk_files = ["epidermal_leaf_trichomes_BRT.nwk", "hiv-1_envelope.nwk", "rubisco_C3_vs_C4.nwk", "HIV_RT_branchlength_1.nwk", "rubisco_reroot_retagged.nwk"]

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

output_path = "results/empirical_data/contrastfel_empirical/time_benchmark.csv"
CSV.write(output_path, df)


#########
# plot ##
#########

plot()
sort!(df, "contrastFEL_time")
scatter(df[!, "difFUBAR_time"], df[!, "contrastFEL_time"], xlabel="difFUBAR time", ylabel="contrastFEL time", label="", gridlinewidth=3)

output_path = "results/empirical_data/contrastfel_empirical/time_plot.png"
savefig(output_path)

plot()
sort!(df, "speedup")
scatter(df[!, "Dataset"], df[!, "speedup"], xlabel="Dataset", ylabel="Speedup",
    gridlinewidth=3, margin=4mm,
    legend=false, series_annotations=text.(round.(df[!, "purity"], digits=2), :bottom))
output_path = "results/empirical_data/contrastfel_empirical/speedup_plot.png"
savefig(output_path)




