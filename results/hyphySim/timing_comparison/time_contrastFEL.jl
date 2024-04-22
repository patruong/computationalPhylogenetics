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
    contrastFEL_output_dir = "contrastFEL_data/omnibus/"
    contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "sims.$sim.settings.replicate.$rep.FEL.json")
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

function get_simulation_node_and_purity(sim)
    newick_file = "contrastFEL_data/omnibus/sims.$sim.nwk"
    treestring_group_labeled, treestring, group_tags, tags, tag_colors = import_grouped_label_tree(newick_file)
    tree = gettreefromnewick(treestring, FelNode)
    n_nodes = length(getnodelist(tree))
    num_groups = length(tag_colors)
    purity, _, _ = CodonMolecularEvolution.get_purity_info(tree, tags, num_groups)
    return n_nodes, purity
end

using Glob

function read_fasta(filename::String)
    fasta_dict = Dict{String,String}()
    current_id = ""
    current_sequence = ""

    for line in eachline(filename)
        if startswith(line, '>')
            # If a line starts with '>', it indicates a new sequence
            if !isempty(current_id)
                fasta_dict[current_id] = current_sequence
            end
            current_id = split(line[2:end])[1]  # Extracting the ID
            current_sequence = ""
        else
            current_sequence *= strip(line)  # Concatenate sequence lines
        end
    end

    # Adding the last sequence to the dictionary
    if !isempty(current_id)
        fasta_dict[current_id] = current_sequence
    end

    return fasta_dict
end

function sequence_count(sequences)
    #all_sequences = vdf_filtered = filter(row -> row["sim"] == "285", df)
    all_sequences = values(sequences)
    num_all_seq = length(all_sequences)
    unique_sequences = unique(all_sequences)
    num_unique_sequences = length(unique_sequences)
    ratio = num_unique_sequences / num_all_seq
    #println("Number of all sequences:", num_all_seq)
    #println("Number of unique sequences:", num_unique_sequences)
    #println("Ratio of unique seqeunces to all sequences:", ratio)
    return num_all_seq, num_unique_sequences, ratio
end

function hyphy_sim_n_seqs(dir, min_unique_seq)
    files = readdir(dir)
    filtered_files = filter(file -> occursin(r"[1-5]$", file), files)

    sim_list = []
    rep_list = []
    n_all_seq_list = []
    n_unique_seq_list = []
    ratio_unique_to_all_list = []
    for file in filtered_files
        fasta_file = dir * "/" * file
        parts = split(fasta_file, '.')
        sim = parts[2]
        rep = parts[5]
        sequences = read_fasta(fasta_file)
        n_all_seq, n_unique_seq, ratio_unique_to_all = sequence_count(sequences)
        push!(sim_list, sim)
        push!(rep_list, rep)
        push!(n_all_seq_list, n_all_seq)
        push!(n_unique_seq_list, n_unique_seq)
        push!(ratio_unique_to_all_list, ratio_unique_to_all)
    end

    df = DataFrame(sim=sim_list, rep=rep_list, n_all_seq=n_all_seq_list, n_unique_seq_list=n_unique_seq_list, ratio_unique_to_all=ratio_unique_to_all_list)

    #df_filtered = filter(row -> row["sim"] == "279", df)
    #df_filtered = filter(row -> row["sim"] == "194", df)
    #df_filtered = filter(row -> row["sim"] == "285", df)

    filtered_df = filter(row -> row["n_unique_seq_list"] >= min_unique_seq, df)
    filtered_df[:, "sim_rep"] = string.(filtered_df.sim, "_", filtered_df.rep)
    return filtered_df
end


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
n_nodes_list = []
purity_list = []
for sim in sims
    for rep in reps
        try
            difFUBAR_res = read_in_difFUBAR_res(pos_thresh, sim, rep)
            contrastFEL_res = read_in_contrastFEL_res(sim, rep)
            difFUBAR_codons = nrow(difFUBAR_res)
            contrastFEL_codons = nrow(contrastFEL_res)
            contrastFEL_time_file_location = "experiments/20240116_rerun_difFUBAR/output_contrastFEL/" * "sim_$sim" * "_replicate_" * "$rep" * "/contrastfel/time_output.txt"
            difFUBAR_time_file_location = "experiments/20240116_rerun_difFUBAR/output_difFUBAR/" * "sim_$sim" * "_replicate_" * "$rep" * "/patrick_max_child/time_output.txt"
            n_node, purity = get_simulation_node_and_purity(sim)
            #difFUBAR_time_file_location = "experiments/20240116_rerun_difFUBAR/output/sim_$sim" * "_replicate_" * "$rep" * "/treesurgery_and_parallel/time_output.txt"

            contrastFEL_real_time, contrastFEL_user_time, contrastFEL_sys_time = parse_time(contrastFEL_time_file_location)
            difFUBAR_real_time, difFUBAR_user_time, difFUBAR_sys_time = parse_time(difFUBAR_time_file_location)
            #difFUBAR_real_time, difFUBAR_user_time, difFUBAR_sys_time = parse_time("experiments/20240116_rerun_difFUBAR/output/sim_$sim" * "_replicate_" * "$rep" * "/treesurgery_and_parallel/time_output.txt")

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
            push!(n_nodes_list, n_node)
            push!(purity_list, purity)

        catch
            println("Error reading files for sim $sim and rep $rep. Skipping this iteration.")
        end
    end
end



df = DataFrame(sim=sim_list, rep=rep_list, contrastFEL_codon_sites=contrastFEL_codon_site, difFUBAR_codon_sites=difFUBAR_codon_site, contrastFEL_real_time=contrastFEL_real_time_list, difFUBAR_real_time=difFUBAR_real_time_list, contrastFEL_user_time=contrastFEL_user_time_list, difFUBAR_user_time=difFUBAR_user_time_list, contrastFEL_sys_time=contrastFEL_sys_time_list, difFUBAR_sys_time=difFUBAR_sys_time_list, n_nodes=n_nodes_list, purity=purity_list)
dir = "contrastFEL_data/omnibus"
hyphy_num_seqs = hyphy_sim_n_seqs(dir, 20)
df[:, "sim_rep"] = string.(df.sim, "_", df.rep)
filter_list = hyphy_num_seqs[:, "sim_rep"]
df = filter(row -> in(row[:sim_rep], filter_list), df)


if sum(df[!, "contrastFEL_codon_sites"] .== df[!, "difFUBAR_codon_sites"]) == nrow(df)
    df[!, "codon_sites"] = df[!, "difFUBAR_codon_sites"]
end

df = sort(df, :codon_sites)




df = sort(df, :difFUBAR_real_time)

# plot difFUBAR vs contrastFEL in one plot
slope = 1.0  # Change this to adjust the steepness of the line
offset = 0.0  # Change this to adjust the y-intercept of the line
x_min, x_max = extrema(df[:, :difFUBAR_real_time])
y_diag = slope * (x_min:x_max) .+ offset
plot(df[!, :difFUBAR_real_time], df[!, :contrastFEL_real_time],
    xlim=(0, 3000), ylim=(0, 80000),  # Set specific ranges for x and y axes
    xlabel="difFUBAR time", ylabel="contrastFEL Time", legend=false) #title="difFUBAR vs Contrast-FEL Time"
plot!(x_min:x_max, y_diag, linestyle=:dot, color=:lightgray)  # Adjust color name for different shades of grey

savefig("results/hyphySim/timing_comparison/difFUBAR_vs_contrastFEL_relative_time.png")  # Saves as PNG by default




index_values = axes(df, 1)
df = sort(df, :codon_sites)
df = sort(df, :n_nodes)
df = sort(df, :purity)

df = sort(df, :difFUBAR_real_time)
plot(index_values=axes(df, 1), df[!, :difFUBAR_real_time],
    xlabel="difFUBAR time", ylabel="contrastFEL Time", legend=false, linecolor=:blue, linewidth=1.5) # title="difFUBAR vs Contrast-FEL Time"

df = sort(df, :contrastFEL_real_time)
plot(index_values=axes(df, 1), df[!, :contrastFEL_real_time],
    xlabel="difFUBAR time", ylabel="contrastFEL Time", legend=false, linecolor=:red, linewidth=1.5) # title="difFUBAR vs Contrast-FEL Time"


# Panel plot version
df = sort(df, :difFUBAR_real_time)
plot(index_values=axes(df, 1), df[!, :difFUBAR_real_time],
    xlabel="Simulations sorted by increasing contrastFEL Time", ylabel="Time (s)", legend=false, linecolor=:blue, linewidth=1.5) # title="difFUBAR vs Contrast-FEL Time"

df = sort(df, :contrastFEL_real_time)

#df
deleteat!(df, 1) #remove first time (sim 1, rep 1) contrastFEL crashes here
df = filter(row -> row["sim"] != 279, df)
plot(index_values=axes(df, 1), df[!, :difFUBAR_real_time],
    xlabel="", ylabel="log\$_{10}\$(s)", legend=false, linecolor=:blue, linewidth=2, size=(700, 300), margin=4mm) #title="difFUBAR vs Contrast-FEL Time"
plot!(index_values=axes(df, 1), df[!, :contrastFEL_real_time],
    xlabel="", ylabel="log\$_{10}\$(s)", legend=false, linecolor=:red, linewidth=2, yscale=:log10)
savefig("results/hyphySim/timing_comparison/difFUBAR_vs_contrastFEL_time_log.png")  # Saves as PNG by default


df = sort(df, :contrastFEL_real_time)
plot(index_values=axes(df, 1), df[!, :difFUBAR_real_time],
    xlabel="Simulations sorted by increasing contrastFEL Time", ylabel="s", legend=false, linecolor=:blue, linewidth=2, size=(700, 300), margin=3mm) #title="difFUBAR vs Contrast-FEL Time"
plot!(index_values=axes(df, 1), df[!, :contrastFEL_real_time],
    xlabel="Simulations sorted by increasing contrastFEL Time", ylabel="s", legend=false, linecolor=:red, linewidth=2)
savefig("results/hyphySim/timing_comparison/difFUBAR_vs_contrastFEL_time.png")  # Saves as PNG by default



# speedup
df[!, :speed_up] = df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time]
maximum(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])
minimum(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])
mean(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])




minimum(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])
min_index = argmin(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])
max_index = argmax(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])


(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])
min_row = df[min_index, :]
max_row = df[max_index, :]


max_row = df[max_index, :]


# Find the maximum index
max_index = argmax(df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time])

# Remove the maximum element from the array
filtered_array = df[!, :contrastFEL_real_time] ./ df[!, :difFUBAR_real_time]
deleteat!(filtered_array, max_index)

# Find the second maximum index
second_max_index = argmax(filtered_array)

# Remove the second maximum element from the array
deleteat!(filtered_array, second_max_index)

# Find the third maximum index
third_max_index = argmax(filtered_array)

# Display the indices
println("Max Index: ", max_index)
println("Second Max Index: ", second_max_index)
println("Third Max Index: ", third_max_index)


df[10, :]
df[10, :]
df[10, :]

#################################
# Mean Time vs Codon Sites plot #
#################################


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
plot(codon_sites, difFUBAR_mean_times, ribbon=difFUBAR_std_times, label="difFUBAR", yerr=difFUBAR_std_times, xlabel="Codon Sites", ylabel="Mean Time", title="Comparison of Mean Time with Error Bars", line=:dot)
plot!(codon_sites, contrastFEL_mean_times, ribbon=contrastFEL_std_times, label="contrastFEL", yerr=contrastFEL_std_times)





# work with this and make some Statistics
# speedup difFUBAR_real_time vs contrastFEL_real_time
# plot average time per codon sites (x-axis codon_sites, y-axis average time)

# Set up a basic difFUBAR 
# Set up a basic triqler run



