cd("/home/ptruong/git/computationalPhylogenetics/")

using Pkg
using CSV
using JSON
using DataFrames
using Plots

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
    df[!, "P(ω1 ≠ ω2)"] = (difFUBAR_res[!, "P(ω1 > ω2)"] .+ difFUBAR_res[!, "P(ω2 > ω1)"])
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

function calculate_TRP_and_FPR_difFUBAR(df)
    df[!, "threshold"] = df[!, "P(ω1 ≠ ω2)"]

    function calculate_TRP_and_FPR_difFUBAR_helper(df, threshold)
        predictions = df[!, "P(ω1 ≠ ω2)"] .> threshold
        actual_difference = df[!, "actual_difference"]
        TP = sum(predictions .& actual_difference)
        TN = sum((predictions .== false) .& (actual_difference .== false))
        FP = sum(predictions .& (actual_difference .== false))
        FN = sum((predictions .== false) .& actual_difference)
        TPR = TP / (TP + FN)  # True Positive Rate (Sensitivity or Recall)
        FPR = FP / (FP + TN)  # False Positive Rate
        return TP, TN, FP, FN, TPR, FPR
    end

    res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR_difFUBAR_helper(df, x.threshold))) => [:TP, :TN, :FP, :FN, :TPR, :FPR])
    return res
end

function calculate_TPR_and_FPR_contrastFEL(df)
    df[!, "threshold"] = df[!, "Q-value (overall)"]
    function calculate_TRP_and_FPR_contrastFEL_helper(df, threshold)
        predictions = df[!, "Q-value (overall)"] .< threshold
        actual_difference = df[!, "actual_difference"]
        TP = sum(predictions .& actual_difference)
        TN = sum((predictions .== false) .& (actual_difference .== false))
        FP = sum(predictions .& (actual_difference .== false))
        FN = sum((predictions .== false) .& actual_difference)
        TPR = TP / (TP + FN)  # True Positive Rate (Sensitivity or Recall)
        FPR = FP / (FP + TN)  # False Positive Rate
        return TP, TN, FP, FN, TPR, FPR
    end
    res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR_contrastFEL_helper(df, x.threshold))) => [:TP, :TN, :FP, :FN, :TPR, :FPR])
    return res
end

# sim.0.replicate.1 

#difFUBAR_Res
difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_0.95/sim.421.replicate.1/REFERENCEvTEST/"
difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)

# contrastFEL use Q-value 
contrastFEL_output_dir = "contrastFEL_data/omnibus/"
contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "sims.421.settings.replicate.1.FEL.json")
contrastFEL_res = add_codon_site(contrastFEL_res)

# Simulator settings
simulator_settings_dir = "contrastFEL_data/omnibus/"
simulator_settings = CSV.File(simulator_settings_dir * "sims.421.settings.tsv", delim="\t", header=true) |> DataFrame
simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")
sum(simulator_settings[!, settings_cols[2]])

difFUBAR_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
contrastFEL_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
difFUBAR_res = sort(difFUBAR_res, "P(ω1 ≠ ω2)", rev=true)
contrastFEL_res = sort!(contrastFEL_res, "Q-value (overall)")

difFUBAR_res = calculate_TRP_and_FPR_difFUBAR(difFUBAR_res)
contrastFEL_res = calculate_TPR_and_FPR_contrastFEL(contrastFEL_res)

plot(res.FPR, res.TPR, xlabel="FPR", ylabel="TPR", label="difFUBAR") # This should be roughly equivalent to FPR?
plot(res.FPR, res.TPR, xlabel="FPR", ylabel="TPR", label="contrastFEL") # This should be roughly equivalent to FPR?



