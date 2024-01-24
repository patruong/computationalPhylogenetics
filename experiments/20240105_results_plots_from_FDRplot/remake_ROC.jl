cd("/home/patrick/git/computationalPhylogenetics/")

using Pkg
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using FreqTables

function id_fill(df, sim, rep)
    # Create an array filled with the desired value
    values = fill("sim.$sim.rep.$rep", nrow(df))
    df[!, "id"] = values
    return df
end
function calculate_ω1_not_equal_ω2(df)
    #The probability that ω1 and ω2 are not equal is the complement of the probability that one is greater than the other:
    #P(ω1 ≠ ω2) = 1 - [P(ω1 > ω2) + P(ω2 > ω1)]
    #This is because there are only two possibilities: either ω1 is greater than ω2, ω2 is greater than ω1, or they are equal. So, you subtract the sum of the probabilities that ω1 is greater than ω2 and ω2 is greater than ω1 from 1 to find the probability that they are not equal.
    df[!, "P(ω1 ≠ ω2)"] = (df[!, "P(ω1 > ω2)"] .+ df[!, "P(ω2 > ω1)"])
    return df
end

#difFUBAR_Res
function read_in_difFUBAR_res(pos_thresh, sim, rep)
    #difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST/"
    #difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST"
    difFUBAR_output_dir = "/home/patrick/git/computationalPhylogenetics/experiments/20240116_rerun_difFUBAR/output/sim_$sim" * "_replicate_$rep/patrick_max_child/"

    #difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
    difFUBAR_res = CSV.File(difFUBAR_output_dir * "_posteriors.csv", delim=",", header=true) |> DataFrame

    difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)
    id_fill(difFUBAR_res, sim, rep)
    return difFUBAR_res
end



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


function add_codon_site(df)
    df[!, :Codon_site] = axes(df, 1)
    return df
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
function add_difference_cols_to_sim_settings(sim_settings, class1="simulator.omega.class0", class2="simulator.omega.class1")
    # sim_settings is pandas dfg
    diff_col = "difference.$class1" * ".vs." * "$class2"
    bool_diff_col = "bool.difference.$class1" * ".vs." * "$class2"
    sim_settings[!, diff_col] = sim_settings[!, class1] .- sim_settings[!, class2]
    sim_settings[!, bool_diff_col] = sim_settings[!, class1] .!= sim_settings[!, class2]
    return sim_settings, [diff_col, bool_diff_col]
end
function read_in_simulator_settings(sim)
    simulator_settings_dir = "contrastFEL_data/omnibus/"
    simulator_settings = CSV.File(simulator_settings_dir * "sims.$sim.settings.tsv", delim="\t", header=true) |> DataFrame
    simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")
    #sum(simulator_settings[!, settings_cols[2]])
    return simulator_settings, settings_cols
end

function calculate_confusion_matrix_vals(df, prob_col, threshold)

    df[!, "predicted_difference"] = df[!, prob_col] .> threshold
    confusion_matrix = freqtable(df, :actual_difference, :predicted_difference)

    #function result = trycatch(expression, default)
    #    try
    #        result = expression;
    #    catch
    #        result = default;
    #    end
    #end

    TN = nothing
    FP = nothing
    FN = nothing
    TP = nothing

    try
        TN = confusion_matrix[1, 1]
    catch
        TN = 0
    end

    try
        FP = confusion_matrix[1, 2]
    catch
        FP = 0
    end

    try
        FN = confusion_matrix[2, 1]
    catch
        FN = 0
    end

    try
        TP = confusion_matrix[2, 2]
    catch
        TP = 0
    end

    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    FNR = FN / (TP + FN)

    return TP, FP, TN, FN, TPR, FPR, FDR, FNR
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
        TP, FP, TN, FN, TPR, FPR, FDR, FNR = calculate_confusion_matrix_vals(df, threshold_col, i)
        push!(result_df, (i, TP, FP, FN, TN, TPR, FPR, FDR, FNR))
    end

    return result_df
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

function find_closest_value(df::DataFrame, column::Symbol, target::Float64)
    values = df[!, column]
    closest_value = argmin(abs.(values .- target))
    return values[closest_value]
end

function filter_df_closest_value(df::DataFrame, target_column::Symbol, target_value::Float64, tolerance::Float64)
    return filter(row -> abs(row.Threshold - find_closest_value(df, :Threshold, target_value)) < tolerance, df)
end



pos_thresh = 0.75
sims = collect(0:550) # sim 500 totally breaks contrast-FEL
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

difFUBAR_res = sort(difFUBAR_res, "P(ω1 ≠ ω2)", rev=true)
contrastFEL_res = sort!(contrastFEL_res, "P-value (overall)")
contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]



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


plot()
for i in 1:length(plot_colors)
    lower_bound = bounds[i][1]
    upper_bound = bounds[i][2]
    difFUBAR_plot = difFUBAR_plot_data[i]
    contrastFEL_plot = contrastFEL_plot_data[i]

    p = plot!(difFUBAR_plot.FPR, difFUBAR_plot.TPR, xlabel="FPR", ylabel="TPR", label="", linecolor=plot_colors[i], linewidth=1.5)
    p = plot!(contrastFEL_plot.FPR, contrastFEL_plot.TPR, label="", line=(:dash, 1.5, plot_colors[i]))

    display(p)

end

difFUBAR_dot_threshold = 0.75
contrastFEL_dot_threshold = 0.95

labels_added = false
for i in 1:length(plot_colors)
    difFUBAR_plot = difFUBAR_plot_data[i]
    contrastFEL_plot = contrastFEL_plot_data[i]

    tolerance = 1e-5
    dot_difFUBAR_plot = filter_df_closest_value(difFUBAR_plot, :Threshold, difFUBAR_dot_threshold, tolerance)
    dot_contrasteFEL_plot = filter_df_closest_value(contrastFEL_plot, :Threshold, contrastFEL_dot_threshold, tolerance)

    if !labels_added
        scatter!([dot_difFUBAR_plot.FPR], [dot_difFUBAR_plot.TPR], label="", markershape=:circle, markercolor=:grey, markersize=8, markerstrokecolor=:grey)
        scatter!([dot_contrasteFEL_plot.FPR], [dot_contrasteFEL_plot.TPR], label="", markershape=:circle, markercolor=:white, markersize=8, markerstrokecolor=:grey)
        labels_added = true  # Set the variable to true once labels are added
    else
        scatter!([dot_difFUBAR_plot.FPR], [dot_difFUBAR_plot.TPR], label="", markershape=:circle, markercolor=:grey, markersize=8, markerstrokecolor=:grey)
        scatter!([dot_contrasteFEL_plot.FPR], [dot_contrasteFEL_plot.TPR], label="", markershape=:circle, markercolor=:white, markersize=8, markerstrokecolor=:grey)
    end
end

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

        scatter!([], [], label="P(ω1 ≠ ω2) > $difFUBAR_dot_threshold", markershape=:circle, markercolor=:grey, markersize=8, markerstrokecolor=:grey)
        scatter!([], [], label="P-value < $contrastfel_dot_label", markershape=:circle, markercolor=:white, markersize=8, markerstrokecolor=:grey)
        legend_added = true  # Set the variable to true once legend entries are added
    end
end


slope = 1
plot!(x -> slope * x, c=:grey, line=:dash, label="", legend=:bottomright)




##############
# FPR PLOT ###
##############


# FPR plot
plot()

fdr_plot_difFUBAR = difFUBAR_res[difFUBAR_res[!, "P(ω1 ≠ ω2)"].>0.75, :]
fdr_plot_contrastFEL = contrastFEL_res[contrastFEL_res[!, "1-Pvalue"].>0.95, :]

#difFUBAR_plot = calculate_ROC_threshold(difFUBAR_res, "P(ω1 ≠ ω2)")
#contrastFEL_plot = calculate_ROC_threshold(contrastFEL_res, "1-Pvalue")

difFUBAR_plot = calculate_ROC_threshold(fdr_plot_difFUBAR, "P(ω1 ≠ ω2)")
contrastFEL_plot = calculate_ROC_threshold(fdr_plot_contrastFEL, "1-Pvalue")


p = plot!(difFUBAR_plot.Threshold, difFUBAR_plot.FPR, xlabel="Posterior Probability Threshold / (1-Pvalue threshold)", ylabel="FPR", label="difFUBAR", linewidth=1.5)
p = plot!(contrastFEL_plot.Threshold, contrastFEL_plot.FPR, label="contrastFEL", line=(:dash, 1.5))

display(p)
slope = 1
plot!(x -> 1 - slope * x, c=:grey, line=:dash, label="diagonal line", legend=:bottomright)
#legend(:bottomright, title="Legend Title", framealpha=0.7)









function noisy(label; λ=0.0)
    if label
        return 1 - λ * rand()
    else
        return λ * rand()
    end
end


labels = rand(Bool, 200);

scores(λ) =
    map(labels) do label
        noisy(label, λ=λ)
    end

using ROC

scores(0.6)
roc_good = roc(scores(0.6), labels, true);
roc_bad = roc(scores(1.0), labels, true);

area_good = AUC(roc_good)
area_bad = AUC(roc_bad)


using Plots
plot(roc_good, label="good");
plot!(roc_bad, label="bad")




