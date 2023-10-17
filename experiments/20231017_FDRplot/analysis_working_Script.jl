cd("/home/ptruong/git/computationalPhylogenetics/")

using Pkg
using CSV
using JSON
using DataFrames

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



# sim.0.replicate.1 

#difFUBAR_Res
difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_0.75/sim.0.replicate.1/REFERENCEvTEST/"
difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)

# contrastFEL use Q-value 
contrastFEL_output_dir = "contrastFEL_data/omnibus/"
contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "sims.0.settings.replicate.1.FEL.json")
contrastFEL_res = add_codon_site(contrastFEL_res)

# Simulator settings
simulator_settings_dir = "contrastFEL_data/omnibus/"
simulator_settings = CSV.File(simulator_settings_dir * "sims.0.settings.tsv", delim="\t", header=true) |> DataFrame
simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")
sum(simulator_settings[!, settings_cols[2]])

difFUBAR_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
contrastFEL_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
difFUBAR_res = sort(difFUBAR_res, "P(ω1 ≠ ω2)", rev=true)
contrastFEL_res = sort!(contrastFEL_res, "Q-value (overall)")




###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
###################################################################


df = difFUBAR_res
# Assuming you have a DataFrame 'df' with the P(ω1 ≠ ω2) column
# and a boolean list 'actual_difference'

# Calculate TPR and FPR
threshold = 0.5  # Set your desired threshold

function calculate_TRP_and_FPR(df, threshold)
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


using DataFrames

# Define your function to accept a DataFrame, a column name for probabilities, and a threshold
function calculate_TRP_and_FPR(df, prob_column, threshold)
    predictions = df[!, prob_column] .> threshold
    actual_difference = df[!, "actual_difference"]
    TP = sum(predictions .& actual_difference)
    TN = sum((predictions .== false) .& (actual_difference .== false))
    FP = sum(predictions .& (actual_difference .== false))
    FN = sum((predictions .== false) .& actual_difference)
    TPR = TP / (TP + FN)  # True Positive Rate (Sensitivity or Recall)
    FPR = FP / (FP + TN)  # False Positive Rate
    return TP, TN, FP, FN, TPR, FPR
end

# Assuming you have a DataFrame 'df', a column 'prob_column', and a threshold
df = DataFrame(
    prob_column=[0.6, 0.7, 0.4, 0.8],
    actual_difference=[true, false, true, false]
)

threshold_value = 0.5  # Set your desired threshold

# Apply the function and store the results in variables
TP, TN, FP, FN, TPR, FPR = calculate_TRP_and_FPR(df, :prob_column, threshold_value)

# Display the results
println("True Positives (TP): $TP")
println("True Negatives (TN): $TN")
println("False Positives (FP): $FP")
println("False Negatives (FN): $FN")
println("True Positive Rate (TPR): $TPR")
println("False Positive Rate (FPR): $FPR")

# Output TPR and FPR
#println("True Positive Rate (TPR): $TPR")
#println("False Positive Rate (FPR): $FPR")



using DataFrames
df = DataFrame(a=1:7, b=4:10, c=7:13)
df[!, "test test"] = df[!, :a]
df[!, "test"] = df[!, :a]

function calculate_result(x, y, df)
    a_sum = sum(df[1:x, :a])
    return a_sum * y
end

function test2(x, y)
    return 2x + y
end
transform(df, AsTable([:a, :c]) => ByRow(x -> (calculate_result(x.a, x.c, df), 4x.c)) => ["ae what", "f"])


df = difFUBAR_res
df[!, "threshold"] = df[!, "P(ω1 ≠ ω2)"]

function calculate_TRP_and_FPR(df, threshold)
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

res = transform(df, AsTable([:threshold]) => ByRow(x -> (calculate_TRP_and_FPR(df, x.threshold))) => [:TP, :TN, :FP, :FN, :TPR, :FPR])
res
using Plots
plot(res.FPR, res.TPR, xlabel="FPR", ylabel="TPR", label="Data") # This should be roughly equivalent to FPR?


#Nice work on this
