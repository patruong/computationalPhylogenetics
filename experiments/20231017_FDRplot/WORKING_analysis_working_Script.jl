cd("/home/ptruong/git/computationalPhylogenetics/")

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

#function read_hyphy_simulator_settings(file_name)
#    df_python = extract_simulator_site_settings_to_tsv.get_simulation_settings_from_hyphy(file_name)
#    df_julia = df_python |> Pandas.DataFrame #|> DataFrames.DataFrame
#    return df_julia
#end

function read_hyphy_simulator_settings(file_name)
    pass
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
    df[!, "P(ω1 ≠ ω2)"] = (df[!, "P(ω1 > ω2)"] .+ df[!, "P(ω2 > ω1)"])
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
        return TP, TN, FP, FN, TPR, FPR
    end

    res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR_difFUBAR_helper(df, x.threshold, p_col))) => [:TP, :TN, :FP, :FN, :TPR, :FPR])
    return res
end


#
#function calculate_TRP_and_FPR_difFUBAR(df)
#    df[!, "threshold"] = df[!, "P(ω1 ≠ ω2)"]
#
#    function calculate_TRP_and_FPR_difFUBAR_helper(df, threshold)
#        predictions = df[!, "P(ω1 ≠ ω2)"] .> threshold
#        actual_difference = df[!, "actual_difference"]
#        TP = sum(predictions .& actual_difference)
#        TN = sum((predictions .== false) .& (actual_difference .== false))
#        FP = sum(predictions .& (actual_difference .== false))
#        FN = sum((predictions .== false) .& actual_difference)
#        TPR = TP / (TP + FN)  # True Positive Rate (Sensitivity or Recall)
#        FPR = FP / (FP + TN)  # False Positive Rate
#        return TP, TN, FP, FN, TPR, FPR
#    end
#
#    res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR_difFUBAR_helper(df, x.threshold))) => [:TP, :TN, :FP, :FN, :TPR, :FPR])
#    return res
#end
#

#function calculate_TPR_and_FPR_contrastFEL(df)
#    df[!, "threshold"] = df[!, "P-value (overall)"]
#
#    function calculate_TRP_and_FPR_contrastFEL_helper(df, threshold)
#        predictions = df[!, "P-value (overall)"] .< threshold
#        actual_difference = df[!, "actual_difference"]
#        TP = sum(predictions .& actual_difference)
#        TN = sum((predictions .== false) .& (actual_difference .== false))
#        FP = sum(predictions .& (actual_difference .== false))
#        FN = sum((predictions .== false) .& actual_difference)
#        TPR = TP / (TP + FN)  # True Positive Rate (Sensitivity or Recall)
#        FPR = FP / (FP + TN)  # False Positive Rate
#        return TP, TN, FP, FN, TPR, FPR
#    end
#
#    res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR_contrastFEL_helper(df, x.threshold))) => [:TP, :TN, :FP, :FN, :TPR, :FPR])
#    return res
#end

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
    difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST/"
    difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
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


function filter_on_effect_size(df, lower_bound, upper_bound)
    predictions = lower_bound .< df[!, "effect_size"] .< upper_bound
    filtered_df = filter(row -> lower_bound < row.effect_size < upper_bound, df)
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



# sim.0.replicate.1 

# This is a good work setting
#pos_thresh = 0.75
#sim = 421
#rep = 3

pos_thresh = 0.75
sim = 421
rep = 3
difFUBAR_res = read_in_difFUBAR_res(pos_thresh, sim, rep)
contrastFEL_res = read_in_contrastFEL_res(sim, rep)
simulator_settings, settings_cols = read_in_simulator_settings(sim)


pos_thresh = 0.95
#sims = [182, 250, 421, 466, 500]
#sims = [182, 250, 421, 466] # sim 500 totally breaks
#reps = [1, 2]
sims = collect(150:250) # sim 500 totally breaks contrast-FEL
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
            difFUBAR_res[!, "actual_effect_difference"] = simulator_settings[!, settings_cols[1]]
            contrastFEL_res[!, "actual_effect_difference"] = simulator_settings[!, settings_cols[1]]
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




# Filtering, should we filter????
#difFUBAR_res = filter(row -> row."P(ω1 ≠ ω2)" < 0.95, difFUBAR_res)
#contrastFEL_res = filter(row -> row."P-value (overall)" < 0.5, contrastFEL_res)

# sorting comes after aggregating batches
difFUBAR_res = sort(difFUBAR_res, "P(ω1 ≠ ω2)", rev=true)
contrastFEL_res = sort!(contrastFEL_res, "P-value (overall)")
contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]

# We can filter on effect size to get the curves Ben has.
calculate_effect_size(difFUBAR_res, "mean(ω1)", "mean(ω2)")
calculate_effect_size(contrastFEL_res, "beta (background)", "beta (TEST)")

# Last step
#difFUBAR_res = calculate_TRP_and_FPR(difFUBAR_res, "P(ω1 ≠ ω2)")
#contrastFEL_res = calculate_TRP_and_FPR(contrastFEL_res, "1-Pvalue")

#difFUBAR_res_list = [filter_on_effect_size(difFUBAR_res, 0, 0.25),
#    filter_on_effect_size(difFUBAR_res, 0.25, 0.5),
#    filter_on_effect_size(difFUBAR_res, 0.5, 3.0),
#    filter_on_effect_size(difFUBAR_res, 3.0, Inf)]
#
#contrastFEL_res_list = [filter_on_effect_size(contrastFEL_res, 0, 0.25),
#    filter_on_effect_size(contrastFEL_res, 0.25, 0.5),
#    filter_on_effect_size(contrastFEL_res, 0.5, 3.0),
#    filter_on_effect_size(contrastFEL_res, 3.0, Inf)]
#    
#difFUBAR_res_list = calculate_TPR_and_FPR_list(difFUBAR_res_list, "P(ω1 ≠ ω2)")
#contrastFEL_res_list = calculate_TPR_and_FPR_list(contrastFEL_res_list, "1-Pvalue")



#calculate_TRP_and_FPR(difFUBAR_res, "P(ω1 ≠ ω2)")
#calculate_TRP_and_FPR(contrastFEL_res, "1-Pvalue")





bounds = [(0, 0.25), (0.25, 0.5), (0.5, 3.0), (3.0, Inf)]
plot_colors = [:red, :blue, :green, :black]

plot()
for i in 1:length(plot_colors)
    lower_bound = bounds[i][1]
    upper_bound = bounds[i][2]
    difFUBAR_plot = calculate_TRP_and_FPR(filter_on_effect_size(difFUBAR_res, lower_bound, upper_bound), "P(ω1 ≠ ω2)")
    contrastFEL_plot = calculate_TRP_and_FPR(filter_on_effect_size(contrastFEL_res, lower_bound, upper_bound), "1-Pvalue")

    #difFUBAR_plot = filter_on_effect_size(difFUBAR_res, lower_bound, upper_bound)
    #contrastFEL_plot = filter_on_effect_size(contrastFEL_res, lower_bound, upper_bound)

    p = plot!(difFUBAR_plot.FPR, difFUBAR_plot.TPR, xlabel="FPR", ylabel="TPR", label="difFUBAR, $lower_bound to $upper_bound", linecolor=plot_colors[i], linewidth=1.5)
    p = plot!(contrastFEL_plot.FPR, contrastFEL_plot.TPR, label="contrastFEL, $lower_bound to $upper_bound", line=(:dash, 1.5, plot_colors[i]))
    display(p)
end
slope = 1
plot!(x -> slope * x, c=:grey, line=:dash, label="diagonal line", legend=:bottomright)
#legend(:bottomright, title="Legend Title", framealpha=0.7)


filter_on_effect_size(contrastFEL_res, 3.0, Inf)

#CSV.write("check_difFUBAR.csv", difFUBAR_res)
#CSV.write("check_contrastFEL.csv", contrastFEL_res)











plot!(filter_on_effect_size(difFUBAR_res, 0, 0.25).FPR, filter_on_effect_size(difFUBAR_res, 0, 0.25).TPR, xlabel="FPR", ylabel="TPR", label="difFUBAR", linecolor=:red, linewidth=2)
plot!(filter_on_effect_size(difFUBAR_res, 0.25, 0.5).FPR, filter_on_effect_size(difFUBAR_res, 0.25, 0.5).TPR, xlabel="FPR", ylabel="TPR", label="difFUBAR", linecolor=:blue, linewidth=2)
plot!(filter_on_effect_size(difFUBAR_res, 0.5, 3.0).FPR, filter_on_effect_size(difFUBAR_res, 0.5, 3.0).TPR, xlabel="FPR", ylabel="TPR", label="difFUBAR", linecolor=:green, linewidth=2)
plot!(filter_on_effect_size(difFUBAR_res, 3.0, Inf).FPR, filter_on_effect_size(difFUBAR_res, 3.0, Inf).TPR, xlabel="FPR", ylabel="TPR", label="difFUBAR", linecolor=:black, linewidth=2)


contrastFEL_res

# Write the code so we are aggregating all the simulated data

########### Cohen d effect size.... but in the plot he uses absolut difference between beta1 and beta2
df = difFUBAR_res
group_1 = "mean(ω1)"
group_2 = "mean(ω2)"
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
calculate_cohen_d_effect_size(df, group_1, group_2)


df = contrastFEL_res
group_1 = "beta (background)"
group_2 = "beta (TEST)"
calculate_cohen_d_effect_size(df, group_1, group_2)


#println("Cohen's d Effect Size: $cohen_d")

# 0.75 effect size 

df =
    difFUBAR_res
contrastFEL_res
# incorporate different effect size filtering.

