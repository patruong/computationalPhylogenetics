#cd("/home/ptruong/git/computationalPhylogenetics/")
cd("/home/patrick/git/computationalPhylogenetics/")

using Pkg
using CSV
using JSON
using DataFrames
using Plots
using Statistics


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
    TN = sum(no_diff[!, "actual_difference"] .== false) #TN

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

# sim.0.replicate.1 

# This is a good work setting
#pos_thresh = 0.75
#sim = 421
#rep = 3

pos_thresh = 0.75
sim = 141
rep = 3
difFUBAR_res = read_in_difFUBAR_res(pos_thresh, sim, rep)
contrastFEL_res = read_in_contrastFEL_res(sim, rep)
simulator_settings, settings_cols = read_in_simulator_settings(sim)


pos_thresh = 0.95
#sims = [182, 250, 421, 466, 500]
#sims = [182, 250, 421, 466] # sim 500 totally breaks
#reps = [1, 2]
sims = collect(0:500) # sim 500 totally breaks contrast-FEL
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

#############################################################
# DO NOT USE BOTH FILTERS, IT WILL FILTER AWAY EVERYTHING   #
# - MAYBE I SHOULD SPLIT UP THIS CODE TO AVOID THIS MISTAKE #
#############################################################


# Filter on w1 > w2
#difFUBAR_res = difFUBAR_res[difFUBAR_res[:, "actual_directional_effect_difference"].>=0, :]
#contrastFEL_res = contrastFEL_res[contrastFEL_res[:, "actual_directional_effect_difference"].>=0, :]

# Filter on w2 > w1
#difFUBAR_res = difFUBAR_res[difFUBAR_res[:, "actual_directional_effect_difference"].<=0, :]
#contrastFEL_res = contrastFEL_res[contrastFEL_res[:, "actual_directional_effect_difference"].<=0, :]



#maximum(difFUBAR_res[!, "P(ω1 ≠ ω2)"])
#minimum(difFUBAR_res[!, "P(ω1 ≠ ω2)"])
#difFUBAR_res[400, :]
#contrastFEL_res[400, :]
# Filtering, should we filter????
#difFUBAR_res = filter(row -> row."P(ω1 ≠ ω2)" < 0.95, difFUBAR_res)
#contrastFEL_res = filter(row -> row."P-value (overall)" < 0.5, contrastFEL_res)

# sorting comes after aggregating batches
difFUBAR_res = sort(difFUBAR_res, "P(ω1 ≠ ω2)", rev=true)
contrastFEL_res = sort!(contrastFEL_res, "P-value (overall)")
contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]

# We can filter on effect size to get the curves Ben has.
#calculate_effect_size(difFUBAR_res, "mean(ω1)", "mean(ω2)")
#calculate_effect_size(contrastFEL_res, "beta (background)", "beta (TEST)")

#lower_bound = 0
#upper_bound = 0.25
#difFUBAR_plot = calculate_ROC_threshold(filter_on(difFUBAR_res, "actual_effect_difference", lower_bound, upper_bound, true), "P(ω1 ≠ ω2)")
#contrastFEL_plot = calculate_ROC_threshold(filter_on(contrastFEL_res, "actual_effect_difference", lower_bound, upper_bound, true), "1-Pvalue")


function find_closest_value(df::DataFrame, column::Symbol, target::Float64)
    values = df[!, column]
    closest_value = argmin(abs.(values .- target))
    return values[closest_value]
end

function filter_df_closest_value(df::DataFrame, target_column::Symbol, target_value::Float64, tolerance::Float64)
    return filter(row -> abs(row.Threshold - find_closest_value(df, :Threshold, target_value)) < tolerance, df)
end


#tolerance = 1e-5
#filter_df_closest_value(difFUBAR_plot, :Threshold, 0.95, tolerance)



#difFUBAR_res = difFUBAR_res[difFUBAR_res[!, "P(ω1 ≠ ω2)"].>0.75, :]
#contrastFEL_res = contrastFEL_res[contrastFEL_res[!, "1-Pvalue"].>0.95, :]


bounds = [(0, 0.25), (0.25, 0.5), (0.5, 3.0), (3.0, Inf)]
plot_colors = [:red, :blue, :green, :black]

labels_added = false

difFUBAR_plot_data = []
contrastFEL_plot_data = []
for i in 1:length(plot_colors)
    lower_bound = bounds[i][1]
    upper_bound = bounds[i][2]
    difFUBAR_plot = calculate_ROC_threshold(filter_on(difFUBAR_res, "actual_effect_difference", lower_bound, upper_bound, true), "P(ω1 ≠ ω2)")
    contrastFEL_plot = calculate_ROC_threshold(filter_on(contrastFEL_res, "actual_effect_difference", lower_bound, upper_bound, true), "1-Pvalue")
    push!(difFUBAR_plot_data, difFUBAR_plot)
    push!(contrastFEL_plot_data, contrastFEL_plot)
end

difFUBAR_dot_threshold = 0.80
contrastFEL_dot_threshold = 0.95

plot(size=(450, 550))

for i in 1:length(plot_colors)
    lower_bound = bounds[i][1]
    upper_bound = bounds[i][2]
    difFUBAR_plot = difFUBAR_plot_data[i]
    contrastFEL_plot = contrastFEL_plot_data[i]

    p = plot!(difFUBAR_plot.FPR, difFUBAR_plot.TPR, xlabel="FPR", ylabel="TPR", label="", linecolor=plot_colors[i], linewidth=1.5)
    p = plot!(contrastFEL_plot.FPR, contrastFEL_plot.TPR, label="", line=(:dash, 1.5, plot_colors[i]))

    display(p)

end

labels_added = false
difFUBAR_dots = []
contrastFEL_dots = []
for i in 1:length(plot_colors)
    difFUBAR_plot = difFUBAR_plot_data[i]
    contrastFEL_plot = contrastFEL_plot_data[i]

    tolerance = 1e-5
    dot_difFUBAR_plot = filter_df_closest_value(difFUBAR_plot, :Threshold, difFUBAR_dot_threshold, tolerance)
    dot_contrasteFEL_plot = filter_df_closest_value(contrastFEL_plot, :Threshold, contrastFEL_dot_threshold, tolerance)
    push!(difFUBAR_dots, dot_difFUBAR_plot)
    push!(contrastFEL_dots, dot_contrasteFEL_plot)

    if !labels_added
        scatter!([dot_difFUBAR_plot.FPR], [dot_difFUBAR_plot.TPR], label="", markershape=:circle, markercolor=:grey, markersize=5, markerstrokecolor=:grey, markeralpha=0.85)
        scatter!([dot_contrasteFEL_plot.FPR], [dot_contrasteFEL_plot.TPR], label="", markershape=:circle, markercolor=:white, markersize=5, markerstrokecolor=:grey, markeralpha=0.85)
        labels_added = true  # Set the variable to true once labels are added
    else
        scatter!([dot_difFUBAR_plot.FPR], [dot_difFUBAR_plot.TPR], label="", markershape=:circle, markercolor=:grey, markersize=5, markerstrokecolor=:grey, markeralpha=0.85)
        scatter!([dot_contrasteFEL_plot.FPR], [dot_contrasteFEL_plot.TPR], label="", markershape=:circle, markercolor=:white, markersize=5, markerstrokecolor=:grey, markeralpha=0.85)
    end
end

#difFUBAR_dots[1][["FPR", "TPR"]]
#difFUBAR_dots.TPR
#contrastFEL_dots.FPR
#contrastFEL_dots.TPR

legend_added = false  # Initialize a boolean variable to track whether legend entries have been added
contrastfel_dot_label = round(1 - contrastFEL_dot_threshold, digits=2)

for i in 1:length(plot_colors)
    if !legend_added
        plot!([], [], line=:solid, linecolor=:black, label="difFUBAR")
        plot!([], [], line=:dash, linecolor=:black, label="contrastFEL")
        plot!([], [], line=:solid, linecolor=:red, label="E = 0.0 to 0.25")
        plot!([], [], line=:solid, linecolor=:blue, label="E = 0.25 to 0.5")
        plot!([], [], line=:solid, linecolor=:green, label="E = 0.5 to 3.0")
        plot!([], [], line=:solid, linecolor=:black, label="E = 3.0 to Inf")

        scatter!([], [], label="P(ω1 ≠ ω2) > $difFUBAR_dot_threshold", markershape=:circle, markercolor=:grey, markersize=4, markerstrokecolor=:grey, markeralpha=0.85)
        scatter!([], [], label="P-value < $contrastfel_dot_label", markershape=:circle, markercolor=:white, markersize=4, markerstrokecolor=:grey, markeralpha=0.85)

        #scatter!([], [], label="P(ω1 ≠ ω2) > $difFUBAR_dot_threshold", markershape=:circle, markercolor=:grey, markersize=2, markerstrokecolor=:grey, markeralpha=0.3)
        #scatter!([], [], label="P-value < $contrastfel_dot_label", markershape=:circle, markercolor=:white, markersize=2, markerstrokecolor=:grey, markeralpha=0.3)
        legend_added = true  # Set the variable to true once legend entries are added
    end
end


slope = 1
plot!(x -> slope * x, c=:grey, line=:dash, label="", legend=:bottomright)
#legend(:bottomright, title="Legend Title", framealpha=0.7)


savefig("results/hyphySim/ROC/ROC.svg")

# FPR plot
plot()
#for i in 1:length(plot_colors)
#    lower_bound = bounds[i][1]
#    upper_bound = bounds[i][2]
#difFUBAR_plot = calculate_TRP_and_FPR(filter_on(difFUBAR_res, "actual_effect_difference", lower_bound, upper_bound), "P(ω1 ≠ ω2)")
#contrastFEL_plot = calculate_TRP_and_FPR(filter_on(contrastFEL_res, "actual_effect_difference", lower_bound, upper_bound), "1-Pvalue")

#difFUBAR_plot = filter_on_effect_size(difFUBAR_res, lower_bound, upper_bound)
#contrastFEL_plot = filter_on_effect_size(contrastFEL_res, lower_bound, upper_bound)
difFUBAR_plot = calculate_ROC_threshold(difFUBAR_res, "P(ω1 ≠ ω2)")
#contrastFEL_plot = calculate_ROC_threshold(filter_on(contrastFEL_res, "actual_effect_difference", lower_bound, upper_bound, true), "1-Pvalue")

plot(size=(550, 300))
p = plot!(difFUBAR_plot.Threshold, difFUBAR_plot.FPR, xlabel="Posterior Probability Threshold", ylabel="FPR", label="difFUBAR False Positive", linecolor=:black, linewidth=1.5)
#p = plot!(contrastFEL_plot.Threshold, contrastFEL_plot.FPR, label="contrastFEL, $lower_bound to $upper_bound", line=(:dash, 1.5, plot_colors[i]))

display(p)
#end
slope = 1
plot!(x -> 1 - slope * x, c=:grey, line=:dash, label="", legend=:topright)
#legend(:bottomright, title="Legend Title", framealpha=0.7)

savefig("results/hyphySim/ROC/TPR_threshold.png")
savefig("results/hyphySim/ROC/TPR_threshold.svg")


difFUBAR_plot

contrastFEL_plot = calculate_ROC_threshold(contrastFEL_res, "1-Pvalue")
p = plot!(contrastFEL_plot.Threshold, contrastFEL_plot.FPR, xlabel="P-value Threshold", ylabel="FPR", label="difFUBAR False Positive", linecolor=:black, linewidth=1.5)


#### investigate FPR plot 
lower_bound = 3.0
upper_bound = Inf
difFUBAR_plot = calculate_ROC_threshold(filter_on(difFUBAR_res, "actual_effect_difference", lower_bound, upper_bound, true), "P(ω1 ≠ ω2)")
contrastFEL_plot = calculate_ROC_threshold(filter_on(contrastFEL_res, "actual_effect_difference", lower_bound, upper_bound, true), "1-Pvalue")

difFUBAR_plot

contrastFEL_plot


# Statistical power plot
# We need to bin the data into 0.5 step intervals to calculate the binomial confidence intervals
using CategoricalArrays
using Distributions
using StatsBase
using Random
using HypothesisTests
using StatsAPI


df = difFUBAR_res
threshold_col = "P(ω1 ≠ ω2)"
threshold = 0.75

function calculate_power(df, threshold_col="P(ω1 ≠ ω2)", threshold=0.75)
    TP, FP, FN, TN, TPR, FPR, FDR, FNR = calculate_ROC_values_threshold(df, threshold_col, threshold)
    power = 1 - FNR
    return power
end



calculate_power(df, "P(ω1 ≠ ω2)", 0.75)
#bin the data into range based on effect_size

bins = cut(df.actual_effect_difference, -0.25:0.5:7, extend=true)


df.effect_size_bins = bins
df = sort(df, "actual_effect_difference")


# Binomial confidence interval as defined by google BARD.
function binomial_ci(s, n, level)
    #lower, upper = binomial_ci(5000000, 10000000, 0.75)
    #println("Lower bound:", lower)
    #println("Upper bound:", upper)

    # Check that s and n are non-negative integers
    if s < 0 || n < 0 || floor(s) != s || floor(n) != n
        error("Invalid input: s and n must be non-negative integers")
    end

    # Compute the maximum likelihood estimate of p
    p_hat = s / n

    # Calculate the lower and upper bounds of the confidence interval
    z = quantile(Normal(), (1 - level) / 2)
    q = sqrt((p_hat * (1 - p_hat)) / n)
    lower = p_hat - z * q
    upper = p_hat + z * q

    return lower, upper
end
#success, a correctly classified actual difference i.e if actual difference is 0 we want to the be false, else we want to the be true
# we check what is classified as actual difference by thresholding on posterior error prob.
function success_condition(df::DataFrame, threshold_col="P(ω1 ≠ ω2)", threshold=0.75)
    df[!, "classified_difference"] = df[:, threshold_col] .> threshold
    df[!, :success] = [row[:actual_difference] == row[:classified_difference] for row in eachrow(df)]
    return df
end

function extract_values_from_bin(bin)
    matches = match(r"\[(-?\d+(?:\.\d+)?),\s*(-?\d+(?:\.\d+)?)\)", bin)
    return (parse(Float64, matches.captures[1]), parse(Float64, matches.captures[2]))
end

function extract_values_from_bin(bin)
    try
        matches = match(r"\[\s*(-?\d+(?:\.\d+)?),?\s*(-?\d+(?:\.\d+)?)\s*\)", bin)
        return (parse(Float64, matches.captures[1]), parse(Float64, matches.captures[2]))
    catch
        matches = match(r"\[\s*(-?\d+(?:\.\d+)?),?\s*(-?\d+(?:\.\d+)?)\s*\]", bin)
        return (parse(Float64, matches.captures[1]), parse(Float64, matches.captures[2]))
    end
end



bins = unique(df.effect_size_bins)
bins[1]
df_bin = df[df.effect_size_bins.==bins[1], :]
power = calculate_power(df_bin, "P(ω1 ≠ ω2)", 0.75)
df_bin = success_condition(df_bin, "P(ω1 ≠ ω2)", 0.75)


# this seems about right
success = sum(df_bin.success)
n = length(df_bin.success)
bin_ci = binomial_ci(success, n, 0.75)
mu_bin_ci = mean(bin_ci)



#power = calculate_power(df_bin, "P(ω1 ≠ ω2)", 0.75)


# the actual plot


plot_x = []
plot_y = []
plot_y_ci = []
for bin in bins
    df_bin = df[df.effect_size_bins.==bin, :]
    df_bin = success_condition(df_bin, "P(ω1 ≠ ω2)", 0.75)
    success = sum(df_bin.success)
    n = length(df_bin.success)

    bin = extract_values_from_bin(bin)
    bin_ci = binomial_ci(success, n, 0.75)
    mean_bin_ci = mean(bin_ci)
    bin_ci = (bin_ci[1] - bin_ci[2]) / 2
    push!(plot_x, mean(bin))
    push!(plot_y, mean_bin_ci)
    push!(plot_y_ci, bin_ci)
end


# Create a new plot with error bars and a line
scatter(plot_x, plot_y, yerr=plot_y_ci, label="data", markercolor=:red)
plot!(plot_x, plot_y, label="fit", linewidth=2)


df_bin = df[df.effect_size_bins.==bins[15], :]
df_bin = success_condition(df_bin, "P(ω1 ≠ ω2)", 0.75)
success = sum(df_bin.success)
n = length(df_bin.success)

### calculate power

plot_x = []
plot_y = []
for bin in bins
    df_bin = df[df.effect_size_bins.==bin, :]
    df_bin = success_condition(df_bin, "P(ω1 ≠ ω2)", 0.75)
    success = sum(df_bin.success)
    n = length(df_bin.success)
    power = calculate_power(df_bin)
    bin = extract_values_from_bin(bin)
    bin_ci = binomial_ci(success, n, 0.75)
    push!(plot_x, mean(bin))
    push!(plot_y, power)
end


#similar
#, title="Clopper-Pearson Intervals"
plot(size=(550, 300))
p = plot!(plot_x, plot_y, label="", xlabel="Effect Size", ylabel="Power", ylim=[0, 1], color=:red)


#### Clopper-Pearson

plot_x = []
plot_y = []
plot_y_lower_bound = []
plot_y_upper_bound = []
for bin in bins
    df_bin = df[df.effect_size_bins.==bin, :]
    df_bin = success_condition(df_bin, "P(ω1 ≠ ω2)", 0.75)
    success = sum(df_bin.success)
    n = length(df_bin.success)
    point_estimate = success / n
    binTest = BinomialTest(success, n, 0.75) # threshold at 0.75
    bin_ci = confint(binTest; level=0.95, tail=:both, method=:clopper_pearson) # confidence interval at 0.95
    bin = extract_values_from_bin(bin)
    push!(plot_x, mean(bin))
    push!(plot_y, point_estimate)
    push!(plot_y_lower_bound, bin_ci[1])
    push!(plot_y_upper_bound, bin_ci[2])
end

using Plots

#p = plot(plot_x, plot_y, label="Point Estimate", xlabel="Effect Size", ylabel="Power", title="Clopper-Pearson Intervals", ylim=[0, 1])
# Plot the vertical line
for i in 1:length(plot_x)
    plot!([plot_x[i], plot_x[i]], [plot_y_lower_bound[i], plot_y_upper_bound[i]], color=:red, label="")
end

display(p)

savefig("results/hyphySim/ROC/clopper_pearson_interval.png")


