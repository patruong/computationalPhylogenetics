
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Revise
using DataFrames

cd("/home/patrick/git/computationalPhylogenetics/")


function read_contrastFEL_output_JSON(file_name)
    # Open the JSON file and read each line as a JSON object
    json_data = JSON.parsefile(file_name)

    content = json_data["MLE"]["content"]["0"]
    headers = [col_name[1] for col_name in json_data["MLE"]["headers"]]


    # We need to transpose the contrast-FEL JSON before being able to get it into dataframe because of the vectors of Any.
    transposed_content = [[] for i in 1:length(content[1])]
    for row_i in 1:length(content)
        row = content[row_i]
        for col_i in 1:length(row)
            push!(transposed_content[col_i], content[row_i][col_i])
        end
    end

    df_values = Dict()
    for (col, vals) in zip(headers, transposed_content)
        df_values[col] = vals
    end

    res = DataFrame(df_values)
    return res
end

function add_difference_cols_to_sim_settings(sim_settings, class1="simulator.omega.class0", class2="simulator.omega.class1")
    # sim_settings is pandas dfg
    diff_col = "difference.$class1" * ".vs." * "$class2"
    bool_diff_col = "bool.difference.$class1" * ".vs." * "$class2"
    sim_settings[!, diff_col] = sim_settings[!, class1] .- sim_settings[!, class2]
    sim_settings[!, bool_diff_col] = sim_settings[!, class1] .!= sim_settings[!, class2]
    return sim_settings, [diff_col, bool_diff_col]
end

function calculate_ω1_not_equal_ω2(df)
    #The probability that ω1 and ω2 are not equal is the complement of the probability that one is greater than the other:
    #P(ω1 ≠ ω2) = 1 - [P(ω1 > ω2) + P(ω2 > ω1)]
    #This is because there are only two possibilities: either ω1 is greater than ω2, ω2 is greater than ω1, or they are equal. So, you subtract the sum of the probabilities that ω1 is greater than ω2 and ω2 is greater than ω1 from 1 to find the probability that they are not equal.
    #df[!, "P(ω1 ≠ ω2)"] = (df[!, "P(ω1 > ω2)"] .+ df[!, "P(ω2 > ω1)"])
    df[!, "P(ω1 ≠ ω2)"] = max.(df[!, "P(ω1 > ω2)"], df[!, "P(ω2 > ω1)"])
    return df
end

function add_codon_site(df)
    df[!, :Codon_site] = axes(df, 1)
    return df
end

function cumsum_actual_positive(df)
    df[!, :cumsum_actual_positive] = cumsum(df[!, :actual_difference])
    return df
end

function cumsum_actual_negative(df) # Actual negatives
    df[!, :cumsum_actual_negative] = cumsum(df[!, :actual_difference] .== false)
    return df
end

function calculate_TRP_and_FPR(df, p_col)
    # p_col is posterior prob. or 1-(p-value)
    #  "P(ω1 ≠ ω2)" or 1 - "P-value (overall)" called ( "1-Pvalue" )
    df[!, "threshold"] = df[!, p_col]

    function calculate_TRP_and_FPR_difFUBAR_helper(df, threshold, p_col)
        predictions = df[!, p_col] .> threshold
        actual_difference = df[!, "actual_difference"]
        TP = sum(predictions .& actual_difference)
        TN = sum((predictions .== false) .& (actual_difference .== false))
        FP = sum(predictions .& (actual_difference .== false))
        FN = sum((predictions .== false) .& actual_difference)
        TPR = TP / (TP + FN)  # True Positive Rate (Sensitivity or Recall)
        FPR = FP / (FP + TN)  # False Positive Rate
        FDR = FP / (FP + TP)
        return TP, TN, FP, FN, TPR, FPR, FDR
    end

    res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR_difFUBAR_helper(df, x.threshold, p_col))) => [:TP, :TN, :FP, :FN, :TPR, :FPR, :FDR])
    return res
end

function calculate_cohen_d_effect_size(df, group_1, group_2)
    # Assuming you have a DataFrame 'df' with columns "mean(ω2)" and "mean(ω1)"
    mean_ω1 = mean(df[!, group_1])
    mean_ω2 = mean(df[!, group_2])
    n_1 = length(df[!, group_1])
    n_2 = length(df[!, group_2])
    std_ω1 = std(df[!, group_1])
    std_ω2 = std(df[!, group_2])
    pooled_std = sqrt(((n_1 - 1) * std_ω1^2 + (n_2 - 1) * std_ω2^2) / (n_1 + n_2 - 2))
    cohen_d = (mean_ω2 - mean_ω1) / pooled_std
    return cohen_d
end

function calculate_effect_size(df, group_1, group_2)
    #df[!, "actual_effect_difference"] = abs.(df[!, group_1] - df[!, group_2])
    df[!, "effect_size"] = abs.(df[!, group_1] - df[!, group_2])

    return df
end


function id_fill(df, sim, rep)
    # Create an array filled with the desired value
    values = fill("sim.$sim.rep.$rep", nrow(df))
    df[!, "id"] = values
    return df
end


# Read in functions

#difFUBAR_Res
function read_in_difFUBAR_res(pos_thresh, sim, rep)
    #difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST/"
    #difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST"
    difFUBAR_output_dir = "/home/patrick/git/computationalPhylogenetics/experiments/20240116_rerun_difFUBAR/output_difFUBAR/sim_$sim" * "_replicate_$rep/patrick_max_child/"

    #difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
    difFUBAR_res = CSV.File(difFUBAR_output_dir * "_posteriors.csv", delim=",", header=true) |> DataFrame

    difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)
    id_fill(difFUBAR_res, sim, rep)
    return difFUBAR_res
end



# contrastFEL use Q-value 
function read_in_contrastFEL_res(sim, rep)
    #contrastFEL_output_dir = "contrastFEL_data/omnibus/"
    #contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "sims.$sim.settings.replicate.$rep.FEL.json")

    contrastFEL_output_dir = "/home/patrick/git/computationalPhylogenetics/experiments/20240116_rerun_difFUBAR/output_contrastFEL"
    contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "/sim_$sim" * "_replicate_$rep" * "/contrastfel/" * "contrastfel.FEL.json")
    contrastFEL_res = add_codon_site(contrastFEL_res)
    id_fill(contrastFEL_res, sim, rep)
    return contrastFEL_res
end

# Simulator settings
function read_in_simulator_settings(sim)
    simulator_settings_dir = "contrastFEL_data/omnibus/"
    simulator_settings = CSV.File(simulator_settings_dir * "sims.$sim.settings.tsv", delim="\t", header=true) |> DataFrame
    simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")
    #sum(simulator_settings[!, settings_cols[2]])
    return simulator_settings, settings_cols
end


function filter_on(df, col_name, lower_bound, upper_bound, keep_null)
    #predictions = lower_bound .< df[!, col_name] .< upper_bound
    filtered_df = filter(row -> lower_bound <= row[col_name] <= upper_bound, df)
    if keep_null == true
        null_df = filter(row -> row[col_name] == 0, df)
        filtered_df = vcat(filtered_df, null_df)
    end
    return filtered_df
end



function calculate_TPR_and_FPR_list(df_list, p_col)
    TPR_and_FPR_list = []
    for res in df_list
        res = calculate_TRP_and_FPR(res, p_col)
        push!(TPR_and_FPR_list, res)
    end
    return TPR_and_FPR_list
end

function calculate_ROC_values_threshold(df, threshold_col, threshold)
    diff = df[df[!, threshold_col].>=threshold, :] # predicted true, no diff
    no_diff = df[df[!, threshold_col].<threshold, :] # predicted false, no diff

    TP = sum(diff[!, "actual_difference"]) #TP
    FP = sum(diff[!, "actual_difference"] .== false) #FP
    FN = sum(no_diff[!, "actual_difference"]) #FN
    TN = sum(no_diff[!, "actual_difference"] .== false) #TNcontrastFEL_res = read_in_contrastFEL_res(sim, rep)


    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    FNR = FN / (TP + FN)
    return TP, FP, FN, TN, TPR, FPR, FDR, FNR
end

function calculate_ROC_threshold(df, threshold_col)
    # Create an empty DataFrame to store the results
    result_df = DataFrame(
        Threshold=Float64[],
        TP=Int[],
        FP=Int[],
        FN=Int[],
        TN=Int[],
        TPR=Float64[],
        FPR=Float64[],
        FDR=Float64[],
        FNR=Float64[]
    )

    threshold_slices = [-Inf, collect(0:0.001:1)..., Inf]
    threshold_slices .+= rand(length(threshold_slices)) .* 0.00000001

    for i in threshold_slices
        TP, FP, FN, TN, TPR, FPR, FDR, FNR = calculate_ROC_values_threshold(df, threshold_col, i)
        push!(result_df, (i, TP, FP, FN, TN, TPR, FPR, FDR, FNR))
    end

    return result_df
end

function merge_difFUBAR_contrastFEL_res(difFUBAR_res, contrastFEL_res)
    contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]
    difFUBAR_res[!, "P(ω1 ≠ ω2)"]
    contrastFEL_res[!, "1-Pvalue"]
    difFUBAR_res[!, "Concatenated_ID"] = difFUBAR_res[!, "id"] .* "_codon.site." .* string.(difFUBAR_res[!, "Codon Sites"])
    contrastFEL_res[!, "Concatenated_ID"] = contrastFEL_res[!, "id"] .* "_codon.site." .* string.(contrastFEL_res[!, "Codon_site"])
    difFUBAR_res = rename(difFUBAR_res, names(difFUBAR_res) .=> ("difFUBAR_" .* string.(names(difFUBAR_res))))
    contrastFEL_res = rename(contrastFEL_res, names(contrastFEL_res) .=> ("contrastFEL_" .* string.(names(contrastFEL_res))))
    joined_res = innerjoin(difFUBAR_res, contrastFEL_res, on=:difFUBAR_Concatenated_ID => :contrastFEL_Concatenated_ID)

    if sum(joined_res[!, "difFUBAR_id"] .== joined_res[!, "contrastFEL_id"]) == nrow(joined_res)
        return joined_res
    else
        println("ERROR in matching row ids!")
        return false
    end
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
    n_nodes = length(getnodelist(tree))
    num_groups = length(tag_colors)
    purity, _, _ = get_purity_info(tree, original_tags, num_groups)
    return purity, n_nodes
end

function list_nwk_files(directory::AbstractString)
    files = readdir(directory)
    nwk_files = filter(file -> endswith(file, ".nwk"), files)
    return nwk_files
end

function get_purity_from_hyphy_sim(directory_path)
    nwk_files = list_nwk_files(directory_path)
    sim_list = []
    purity_list = []
    nodes_list = []
    for file in nwk_files
        purity, n_nodes = get_purity_from_newick(directory_path * "/" * file)
        match_result = match(r"\d+", file)
        integer_part = parse(Int, match_result.match)
        push!(purity_list, purity)
        push!(sim_list, integer_part)
        push!(nodes_list, n_nodes)
    end

    df = DataFrame(sim=sim_list, purity=purity_list, nodes=nodes_list)
    return df
end

directory_path = "contrastFEL_data/omnibus"

purity_df = get_purity_from_hyphy_sim(directory_path)

pos_thresh = 0.95
sims = collect(0:500)
reps = collect(1:5)


sim_list = []
rep_list = []
contrastFEL_codon_site = []
difFUBAR_codon_site = []
contrastFEL_real_time_list = []
contrastFEL_user_time_list = []
contrastFEL_sys_time_list = []
difFUBAR_real_time_list = []
difFUBAR_user_time_list = []
difFUBAR_sys_time_list = []

# Initialize variables to store aggregated results
aggregated_difFUBAR_res = DataFrame()
aggregated_contrastFEL_res = DataFrame()






for sim in sims
    for rep in reps
        try
            difFUBAR_res = read_in_difFUBAR_res(pos_thresh, sim, rep)
            contrastFEL_res = read_in_contrastFEL_res(sim, rep)

            # aggregated dataframe data
            simulator_settings, settings_cols = read_in_simulator_settings(sim)
            difFUBAR_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
            contrastFEL_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
            difFUBAR_res[!, "actual_effect_difference"] = abs.(simulator_settings[!, settings_cols[1]])
            contrastFEL_res[!, "actual_effect_difference"] = abs.(simulator_settings[!, settings_cols[1]])
            difFUBAR_res[!, "actual_directional_effect_difference"] = (simulator_settings[!, settings_cols[1]])
            contrastFEL_res[!, "actual_directional_effect_difference"] = (simulator_settings[!, settings_cols[1]])

            # time data
            difFUBAR_codons = nrow(difFUBAR_res)
            contrastFEL_codons = nrow(contrastFEL_res)
            contrastFEL_time_file_location = "experiments/20240116_rerun_difFUBAR/output_contrastFEL/" * "sim_$sim" * "_replicate_" * "$rep" * "/contrastfel/time_output.txt"
            difFUBAR_time_file_location = "experiments/20240116_rerun_difFUBAR/output_difFUBAR/" * "sim_$sim" * "_replicate_" * "$rep" * "/patrick_max_child/time_output.txt"

            contrastFEL_real_time, contrastFEL_user_time, contrastFEL_sys_time = parse_time(contrastFEL_time_file_location)
            difFUBAR_real_time, difFUBAR_user_time, difFUBAR_sys_time = parse_time(difFUBAR_time_file_location)

            push!(contrastFEL_codon_site, contrastFEL_codons)
            push!(difFUBAR_codon_site, difFUBAR_codons)
            push!(sim_list, sim)
            push!(rep_list, rep)
            push!(contrastFEL_real_time_list, contrastFEL_real_time)
            push!(contrastFEL_user_time_list, contrastFEL_user_time)
            push!(contrastFEL_sys_time_list, contrastFEL_sys_time)
            push!(difFUBAR_real_time_list, difFUBAR_real_time)
            push!(difFUBAR_user_time_list, difFUBAR_user_time)
            push!(difFUBAR_sys_time_list, difFUBAR_sys_time)
        catch
            println("Error reading files for sim $sim and rep $rep. Skipping this iteration.")
        end
    end
end

df = DataFrame(sim=sim_list, rep=rep_list, contrastFEL_codon_sites=contrastFEL_codon_site, difFUBAR_codon_sites=difFUBAR_codon_site, contrastFEL_real_time=contrastFEL_real_time_list, difFUBAR_real_time=difFUBAR_real_time_list, contrastFEL_user_time=contrastFEL_user_time_list, difFUBAR_user_time=difFUBAR_user_time_list, contrastFEL_sys_time=contrastFEL_sys_time_list, difFUBAR_sys_time=difFUBAR_sys_time_list)
df = innerjoin(df, purity_df, on=:sim)


if sum(df[!, "contrastFEL_codon_sites"] .== df[!, "difFUBAR_codon_sites"]) == nrow(df)
    df[!, "codon_sites"] = df[!, "difFUBAR_codon_sites"]
end

df = sort(df, :codon_sites)

codon_sites_list = []
contrastFEL_mean_time_list = []
contrastFEL_std_time_list = []
difFUBAR_mean_time_list = []
difFUBAR_std_time_list = []

for i in unique(df[!, "contrastFEL_codon_sites"])
    contrastFEL_mean_time = mean(df[df[!, "codon_sites"].==i, :][!, "contrastFEL_real_time"])
    contrastFEL_std_time = std(df[df[!, "codon_sites"].==i, :][!, "contrastFEL_real_time"])
    difFUBAR_mean_time = mean(df[df[!, "codon_sites"].==i, :][!, "difFUBAR_real_time"])
    difFUBAR_std_time = std(df[df[!, "codon_sites"].==i, :][!, "difFUBAR_real_time"])
    push!(codon_sites_list, i)
    push!(contrastFEL_mean_time_list, contrastFEL_mean_time)
    push!(contrastFEL_std_time_list, contrastFEL_std_time)
    push!(difFUBAR_mean_time_list, difFUBAR_mean_time)
    push!(difFUBAR_std_time_list, difFUBAR_std_time)
end

df[!, "sim"]

mean_time_df = DataFrame(codon_sites=codon_sites_list, contrastFEL_mean_time=contrastFEL_mean_time_list, contrastFEL_std_time=contrastFEL_std_time_list, difFUBAR_mean_time=difFUBAR_mean_time_list, difFUBAR_std_time=difFUBAR_std_time_list)

using Plots

# Assuming your DataFrame is called df

# Extracting data for difFUBAR
difFUBAR_mean_times = mean_time_df[!, "difFUBAR_mean_time"]
difFUBAR_std_times = mean_time_df[!, "difFUBAR_std_time"]

# Extracting data for contrastFEL
contrastFEL_mean_times = mean_time_df[!, "contrastFEL_mean_time"]
contrastFEL_std_times = mean_time_df[!, "contrastFEL_std_time"]

# Codon sites
codon_sites = mean_time_df[!, "codon_sites"]

# Create the line plot with error bars
plot(codon_sites, difFUBAR_mean_times, ribbon=difFUBAR_std_times, label="difFUBAR", yerr=difFUBAR_std_times, xlabel="Codon Sites", ylabel="Mean Time", title="difFUBAR
Mean Time with Error Bars", size=(800, 600))
savefig("results/hyphySim/timing_comparison/detailed_timing/difFUBAR_meanTime_vs_site.png")

plot(codon_sites, contrastFEL_mean_times, ribbon=contrastFEL_std_times, label="contrastFEL", yerr=contrastFEL_std_times, xlabel="Codon Sites", ylabel="Mean Time", title="ContrastFEL Mean Time with Error Bars", size=(800, 600))
savefig("results/hyphySim/timing_comparison/detailed_timing/contrastFEL_meanTime_vs_site.png")

plot(codon_sites, difFUBAR_mean_times, ribbon=difFUBAR_std_times, label="difFUBAR", yerr=difFUBAR_std_times, xlabel="Codon Sites", ylabel="Mean Time", title="Comparison of 
Mean Time with Error Bars", line=:dot, size=(800, 600))
plot!(codon_sites, contrastFEL_mean_times, ribbon=contrastFEL_std_times, label="contrastFEL", yerr=contrastFEL_std_times)
savefig("results/hyphySim/timing_comparison/detailed_timing/meanTime_vs_site.png")
################
#
###############

aggregated_df = combine(groupby(df, [:codon_sites])) do group
    DataFrame(
        mean_contrastFEL_real_time=mean(group[:, :contrastFEL_real_time]),
        std_contrastFEL_real_time=std(group[:, :contrastFEL_real_time]),
        mean_difFUBAR_real_time=mean(group[:, :difFUBAR_real_time]),
        std_difFUBAR_real_time=std(group[:, :difFUBAR_real_time]),
        # You can add more aggregation functions if needed
    )
end

sort!(aggregated_df, :codon_sites)
# Extracting columns from the DataFrame
codon_sites = aggregated_df[:, :codon_sites]
mean_difFUBAR_real_time = aggregated_df[:, :mean_difFUBAR_real_time]
std_difFUBAR_real_time = aggregated_df[:, :std_difFUBAR_real_time]
mean_contrastFEL_real_time = aggregated_df[:, :mean_contrastFEL_real_time]
std_contrastFEL_real_time = aggregated_df[:, :std_contrastFEL_real_time]

plot(codon_sites, mean_difFUBAR_real_time, ribbon=std_difFUBAR_real_time, label="DifFUBAR", yerr=std_difFUBAR_real_time, xlabel="Codon Sites", ylabel="Mean Time", title="Comparison of Mean Time with Error Bars", line=:dot, size=(800, 600))
plot!(codon_sites, mean_contrastFEL_real_time, ribbon=std_contrastFEL_real_time, label="ContrastFEL", yerr=std_contrastFEL_real_time)
savefig("results/hyphySim/timing_comparison/detailed_timing/meanTime_vs_site.png")

############
# Purity ###
############

aggregated_df = combine(groupby(df, [:purity])) do group
    DataFrame(
        mean_contrastFEL_real_time=mean(group[:, :contrastFEL_real_time]),
        std_contrastFEL_real_time=std(group[:, :contrastFEL_real_time]),
        mean_difFUBAR_real_time=mean(group[:, :difFUBAR_real_time]),
        std_difFUBAR_real_time=std(group[:, :difFUBAR_real_time]),
        # You can add more aggregation functions if needed
    )
end

sort!(aggregated_df, :purity)

aggregated_df
# Extracting columns from the DataFrame
purity = aggregated_df[:, :purity]
mean_difFUBAR_real_time = aggregated_df[:, :mean_difFUBAR_real_time]
std_difFUBAR_real_time = aggregated_df[:, :std_difFUBAR_real_time]
mean_contrastFEL_real_time = aggregated_df[:, :mean_contrastFEL_real_time]
std_contrastFEL_real_time = aggregated_df[:, :std_contrastFEL_real_time]


plot(purity, mean_difFUBAR_real_time, ribbon=std_difFUBAR_real_time, label="difFUBAR", yerr=std_difFUBAR_real_time, xlabel="Purity", ylabel="Mean Time", title="difFUBAR mean times as a function of purity", size=(800, 600))
savefig("results/hyphySim/timing_comparison/detailed_timing/difFUBAR_meanTime_vs_purity.png")

plot(purity, mean_contrastFEL_real_time, ribbon=std_contrastFEL_real_time, label="ContrastFEL", yerr=std_contrastFEL_real_time, xlabel="Purity", ylabel="Mean Time", title="difFUBAR mean times as a function of purity", size=(800, 600))
savefig("results/hyphySim/timing_comparison/detailed_timing/contrastFEL_meanTime_vs_purity.png")

plot(purity, mean_difFUBAR_real_time, ribbon=std_difFUBAR_real_time, label="DifFUBAR", yerr=std_difFUBAR_real_time, xlabel="Purity", ylabel="Mean Time", title="Mean times as a function of purity", size=(800, 600))
plot!(purity, mean_contrastFEL_real_time, ribbon=std_contrastFEL_real_time, label="ContrastFEL", yerr=std_contrastFEL_real_time)
savefig("results/hyphySim/timing_comparison/detailed_timing/meanTime_vs_purity.png")



###########
# n nodes #
###########

aggregated_df = combine(groupby(df, [:nodes])) do group
    DataFrame(
        mean_contrastFEL_real_time=mean(group[:, :contrastFEL_real_time]),
        std_contrastFEL_real_time=std(group[:, :contrastFEL_real_time]),
        mean_difFUBAR_real_time=mean(group[:, :difFUBAR_real_time]),
        std_difFUBAR_real_time=std(group[:, :difFUBAR_real_time]),
        # You can add more aggregation functions if needed
    )
end

sort!(aggregated_df, :nodes)

aggregated_df
# Extracting columns from the DataFrame
nodes = aggregated_df[:, :nodes]
mean_difFUBAR_real_time = aggregated_df[:, :mean_difFUBAR_real_time]
std_difFUBAR_real_time = aggregated_df[:, :std_difFUBAR_real_time]
mean_contrastFEL_real_time = aggregated_df[:, :mean_contrastFEL_real_time]
std_contrastFEL_real_time = aggregated_df[:, :std_contrastFEL_real_time]

plot(nodes, mean_difFUBAR_real_time, ribbon=std_difFUBAR_real_time, label="DifFUBAR", yerr=std_difFUBAR_real_time, xlabel="Number of nodes", ylabel="Mean Time", title="difFUBAR mean times as a function of nodes", size=(800, 600))
savefig("results/hyphySim/timing_comparison/detailed_timing/difFUBAR_meanTime_vs_nodes.png")

plot(nodes, mean_contrastFEL_real_time, ribbon=std_contrastFEL_real_time, label="ContrastFEL", yerr=std_contrastFEL_real_time, title="contrastFEL mean times as a function of nodes", xlabel="Number of nodes", ylabel="Mean Time", size=(800, 600))
savefig("results/hyphySim/timing_comparison/detailed_timing/contrastFEL_meanTime_vs_nodes.png")

plot(nodes, mean_difFUBAR_real_time, ribbon=std_difFUBAR_real_time, label="DifFUBAR", yerr=std_difFUBAR_real_time, xlabel="Number of nodes", ylabel="Mean Time", title="difFUBAR mean times as a function of nodes", size=(800, 600))
plot!(nodes, mean_contrastFEL_real_time, ribbon=std_contrastFEL_real_time, label="ContrastFEL", yerr=std_contrastFEL_real_time, title="Mean times as a function of nodes")
savefig("results/hyphySim/timing_comparison/detailed_timing/meanTime_vs_nodes.png")


#######################
# SURFACE PLOTS       #
# time(purity, nodes) #
#######################

aggregated_df = combine(groupby(df, [:nodes, :purity])) do group

    DataFrame(
        mean_contrastFEL_real_time=mean(group[:, :contrastFEL_real_time]),
        std_contrastFEL_real_time=std(group[:, :contrastFEL_real_time]),
        mean_difFUBAR_real_time=mean(group[:, :difFUBAR_real_time]),
        std_difFUBAR_real_time=std(group[:, :difFUBAR_real_time]),
        # You can add more aggregation functions if needed
    )
end

using Plots;
pyplot();

# Sort the DataFrame
sort!(aggregated_df, [:purity, :nodes])

# Create the surface plot
#     title="Mean difFUBAR Real Time vs. Nodes and Purity",
surface_plot = plot(
    aggregated_df.purity,
    aggregated_df.nodes,
    aggregated_df.mean_difFUBAR_real_time,
    st=:surface,
    xlabel="Purity",
    ylabel="Nodes",
    zlabel="Mean difFUBAR Real Time",
    zguidefontrotation=90,
    size=(800, 600)
)
savefig("results/hyphySim/timing_comparison/detailed_timing/difFUBAR_surface_realtime_nodes_purity.png")

# Show the plot
display(surface_plot)


# Create the surface plot
#title="Mean contrastFEL Real Time vs. Nodes and Purity",

surface_plot = plot(
    aggregated_df.purity,
    aggregated_df.nodes,
    aggregated_df.mean_contrastFEL_real_time,
    st=:surface,
    xlabel="Purity",
    ylabel="Nodes",
    zlabel="Mean contrastFEL Real Time",
    zguidefontrotation=90, size=(800, 600)
)
savefig("results/hyphySim/timing_comparison/detailed_timing/contrastFEL_surface_realtime_nodes_purity.png")


# Show the plot
display(surface_plot)


#############################
# SURFACE PLOTS             #
# time(purity, codon sites) #
#############################

aggregated_df = combine(groupby(df, [:codon_sites, :purity])) do group

    DataFrame(
        mean_contrastFEL_real_time=mean(group[:, :contrastFEL_real_time]),
        std_contrastFEL_real_time=std(group[:, :contrastFEL_real_time]),
        mean_difFUBAR_real_time=mean(group[:, :difFUBAR_real_time]),
        std_difFUBAR_real_time=std(group[:, :difFUBAR_real_time]),
        # You can add more aggregation functions if needed
    )
end

using Plots;
pyplot();

# Sort the DataFrame
sort!(aggregated_df, [:purity, :codon_sites])

# Create the surface plot
#     title="Mean difFUBAR Real Time vs. Codon sites and Purity",
surface_plot = plot(
    aggregated_df.purity,
    aggregated_df.codon_sites,
    aggregated_df.mean_difFUBAR_real_time,
    st=:surface,
    xlabel="Purity",
    ylabel="Codon sites",
    zlabel="Mean difFUBAR time",
    xguidefontrotation=345, xguidefontsize=16,
    yguidefontrotation=45, yguidefontsize=16,
    zguidefontrotation=90, zguidefontsize=16,
    size=(800, 600),
    tickfont=font(8)
)


savefig("results/hyphySim/timing_comparison/detailed_timing/difFUBAR_surface_realtime_codonSites_purity.png")

# Show the plot
display(surface_plot)


# Create the surface plot
#     title="Mean contrastFEL Real Time vs. Codon sites and Purity",
surface_plot = plot(
    aggregated_df.purity,
    aggregated_df.codon_sites,
    aggregated_df.mean_contrastFEL_real_time,
    st=:surface,
    xlabel="Purity",
    ylabel="Codon sites",
    zlabel="Mean contrast-FEL time",
    xguidefontrotation=345, xguidefontsize=16,
    yguidefontrotation=45, yguidefontsize=16,
    zguidefontrotation=90, zguidefontsize=16,
    size=(800, 600),
    tickfont=font(8)
)

savefig("results/hyphySim/timing_comparison/detailed_timing/contrastFEL_surface_realtime_codonSites_purity.png")

# Show the plot
display(surface_plot)

#############################
# SURFACE PLOTS             #
# time(nodes, codon sites) #
#############################

aggregated_df = combine(groupby(df, [:codon_sites, :nodes])) do group

    DataFrame(
        mean_contrastFEL_real_time=mean(group[:, :contrastFEL_real_time]),
        std_contrastFEL_real_time=std(group[:, :contrastFEL_real_time]),
        mean_difFUBAR_real_time=mean(group[:, :difFUBAR_real_time]),
        std_difFUBAR_real_time=std(group[:, :difFUBAR_real_time]),
        # You can add more aggregation functions if needed
    )
end

using Plots;
pyplot();

# Sort the DataFrame
sort!(aggregated_df, [:nodes, :codon_sites])

# Create the surface plot
#     title="Mean difFUBAR Real Time vs. Codon sites and Nodes",

surface_plot = plot(
    aggregated_df.codon_sites,
    aggregated_df.nodes,
    aggregated_df.mean_difFUBAR_real_time,
    st=:surface,
    xlabel="Codon sites",
    ylabel="Nodes",
    zlabel="Mean difFUBAR time",
    xguidefontrotation=345, xguidefontsize=16,
    yguidefontrotation=45, yguidefontsize=16,
    zguidefontrotation=90, zguidefontsize=16,
    size=(800, 600),
    tickfont=font(8),
    xflip=true
)
savefig("results/hyphySim/timing_comparison/detailed_timing/difFUBAR_surface_realtime_codonSites_nodes.png")

# Show the plot
display(surface_plot)


# Create the surface plot
#title="Mean contrastFEL Real Time vs. Codon sites and Nodes",

surface_plot = plot(
    aggregated_df.codon_sites,
    aggregated_df.nodes,
    aggregated_df.mean_contrastFEL_real_time,
    st=:surface,
    xlabel="Codon sites",
    ylabel="Nodes",
    zlabel="Mean contrastFEL time",
    xflip=true,
    xguidefontrotation=345, xguidefontsize=16,
    yguidefontrotation=45, yguidefontsize=16,
    zguidefontrotation=90, zguidefontsize=16,
    size=(800, 600),
    tickfont=font(8),)
savefig("results/hyphySim/timing_comparison/detailed_timing/contrastFEL_surface_realtime_codonSites_nodes.png")

# Show the plot
display(surface_plot)

n_sims = nrow(df)
nrow(df[df[!, :purity].>0.25, :])
nrow(df[df[!, :purity].>0.50, :])
nrow(df[df[!, :purity].>0.75, :])

nrow(df[df[!, :purity].>0.25, :]) / n_sims
nrow(df[df[!, :purity].>0.50, :]) / n_sims
nrow(df[df[!, :purity].>0.75, :]) / n_sims

# Make this save the plots to results in my drive