
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using KernelDensity

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


pos_thresh = 0.95
sims = collect(0:500)
reps = collect(1:5)

# Initialize variables to store aggregated results
aggregated_difFUBAR_res = DataFrame()
aggregated_contrastFEL_res = DataFrame()

for sim in sims
    for rep in reps
        try
            difFUBAR_res = read_in_difFUBAR_res(pos_thresh, sim, rep)
            contrastFEL_res = read_in_contrastFEL_res(sim, rep)
            simulator_settings, settings_cols = read_in_simulator_settings(sim)
            difFUBAR_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
            contrastFEL_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
            difFUBAR_res[!, "actual_effect_difference"] = abs.(simulator_settings[!, settings_cols[1]])
            contrastFEL_res[!, "actual_effect_difference"] = abs.(simulator_settings[!, settings_cols[1]])
            difFUBAR_res[!, "actual_directional_effect_difference"] = (simulator_settings[!, settings_cols[1]])
            contrastFEL_res[!, "actual_directional_effect_difference"] = (simulator_settings[!, settings_cols[1]])
            # Append the current batch of data to the aggregated DataFrames
            append!(aggregated_difFUBAR_res, difFUBAR_res)
            append!(aggregated_contrastFEL_res, contrastFEL_res)
        catch
            println("Error reading files for sim $sim and rep $rep. Skipping this iteration.")
        end
    end
end


difFUBAR_res = aggregated_difFUBAR_res
contrastFEL_res = aggregated_contrastFEL_res
contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]
contrastFEL_res[!, "Codon Sites"] = contrastFEL_res[!, "Codon_site"]

joined_res = merge_difFUBAR_contrastFEL_res(difFUBAR_res, contrastFEL_res)


function get_variance_from_sim(df, col_name, sim)
    filtered_df = df[contains.(df.id, "sim.$sim"), :]
    filtered_df = filtered_df[filtered_df[!, "Codon Sites"].==1, :]
    variance = var(filtered_df[!, col_name])
    return variance
end


function get_sim_variances(df, col_name, sims)
    variances = []
    for sim in sims
        sim = string(sim)
        variance = get_variance_from_sim(df, col_name, sim)
        push!(variances, variance)
    end
    return variances
end


# plot
difFUBAR_var = get_sim_variances(difFUBAR_res, "P(ω1 ≠ ω2)", sims)
contrastFEL_var = get_sim_variances(contrastFEL_res, "1-Pvalue", sims)

variance_df = DataFrame(difFUBAR_var=difFUBAR_var, contrastFEL_var=contrastFEL_var)

variance_df.difFUBAR_var = convert(Vector{Float64}, variance_df.difFUBAR_var)
variance_df.contrastFEL_var = convert(Vector{Float64}, variance_df.contrastFEL_var)

variance_df = variance_df[.!isnan.(variance_df.difFUBAR_var), :]
variance_df = variance_df[.!isnan.(variance_df.contrastFEL_var), :]

difFUBAR_density = kde(variance_df.difFUBAR_var)
contrastFEL_density = kde(variance_df.contrastFEL_var)

# Plot the distributions using KDEs with truncated x-axis at 0
plot(difFUBAR_density.x, difFUBAR_density.density .* 0.01, label="difFUBAR P(ω1 ≠ ω2)", xlims=(0, Inf), linewidth=2)
plot!(contrastFEL_density.x, contrastFEL_density.density .* 0.01, label="contrastFEL 1-Pvalue", xlims=(0, Inf), linewidth=2)

# Add labels and title
xlabel!("Variance of P(ω1 ≠ ω2) / 1-Pvalue")
ylabel!("Density")
title!("Variance of P(ω1 ≠ ω2) / 1-Pvalue")

# Display the plot
plot!(legend=:topright)  # Adjust legend position if needed

plot_dir = "/home/patrick/git/computationalPhylogenetics/results/hyphySim/distribution_of_variation/"
savefig(plot_dir * "distribution_of_variation.png")

########################
# p-value distribution #
########################

histogram(difFUBAR_res[!, "P(ω1 ≠ ω2)"], label="P(ω1 ≠ ω2)", alpha=0.5, density=true)
xlabel!("P(ω1 ≠ ω2)")
ylabel!("Count")
title!("difFUBAR P(ω1 ≠ ω2) distribution")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "difFUBAR_posterior_distribution.png")

histogram(contrastFEL_res[!, "1-Pvalue"], label="1-Pvalue", alpha=0.5, density=true)
xlabel!("1-Pvalue")
ylabel!("Count")
title!("contrastFEL 1-Pvalue distribution")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "contrastFEL_pvalue_distribution.png")

histogram(difFUBAR_res[!, "P(ω1 ≠ ω2)"], label="P(ω1 ≠ ω2)", alpha=0.5, density=true)
histogram!(contrastFEL_res[!, "1-Pvalue"], label="1-Pvalue", alpha=0.5, density=true)
xlabel!("P(ω1 ≠ ω2) / 1-Pvalue")
ylabel!("Count")
title!("P(ω1 ≠ ω2) / 1-Pvalue - distribution")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "contrastFEL_pvalue_distribution.png")