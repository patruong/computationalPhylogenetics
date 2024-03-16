
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
using MolecularEvolution
using CodonMolecularEvolution

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
function read_in_difFUBAR_res(filename)
    difFUBAR_res = CSV.File(filename, delim=",", header=true) |> DataFrame
    difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)
    #id_fill(difFUBAR_res, sim, rep)
    return difFUBAR_res
end



# contrastFEL use Q-value 
function read_in_contrastFEL_res(filename)
    contrastFEL_res = read_contrastFEL_output_JSON(filename)
    contrastFEL_res = add_codon_site(contrastFEL_res)
    #id_fill(contrastFEL_res, sim, rep)
    return contrastFEL_res
end

datasets_folder = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT_branchlength_1"]
datasets_name = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT"]

range = (0.01:0.01:1)
results = []
for i in 1:length(datasets_folder)
    dataset = datasets_folder[i]
    col_name = datasets_name[i]
    difFUBAR_location = "experiments/20240308_contrastFEL_empirical/results/$dataset/treesurgery_and_parallel/analysis_posteriors.csv"
    contrastFEL_location = "experiments/20240308_contrastFEL_empirical/results/$dataset/contrastfel/contrastfel.FEL.json"
    difFUBAR_res = read_in_difFUBAR_res(difFUBAR_location)
    contrastFEL_res = read_in_contrastFEL_res(contrastFEL_location)
    contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]

    difFUBAR_count = []
    contrastFEL_count = []
    for i in range
        n_difFUBAR = nrow(difFUBAR_res[difFUBAR_res[!, "P(ω1 ≠ ω2)"].>i, :])
        n_contrastFEL = nrow(contrastFEL_res[contrastFEL_res[!, "1-Pvalue"].>i, :])
        push!(difFUBAR_count, n_difFUBAR)
        push!(contrastFEL_count, n_contrastFEL)
    end

    res = DataFrame(threshold=range)
    res[!, "difFUBAR_"*col_name] = difFUBAR_count
    res[!, "contrastFEL_"*col_name] = contrastFEL_count
    push!(results, res)
end


results = reduce((df1, df2) -> innerjoin(df1, df2, on=:threshold), results)

results[!, "P(ω1 ≠ ω2) / (1-Pvalue)"] = (results[!, "threshold"] .- 1) .* -1


colors = Dict(
    "epidermal_leaf" => :blue,
    "hiv_envelope" => :green,
    "rubisco_C3_vs_C4" => :red,
    "hivRT" => :orange
)

line_styles = Dict(
    "difFUBAR" => :solid,
    "contrastFEL" => :dash
)

datasets_name = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT"]
methods = ["difFUBAR", "contrastFEL"]

plot()

for (i, dataset) in enumerate(datasets_name)
    for (j, method) in enumerate(methods)
        col_name = Symbol("$(method)_$(dataset)")
        data = results[!, col_name]
        #plot!(results[!,  "P(ω1 ≠ ω2) / (1-Pvalue)"], data, label = i == 1 ? "$dataset $method" : "", color = colors[dataset],
        #      linestyle = line_styles[method], xlabel = "Threshold", ylabel = "Value")
        plot!(results[!, "P(ω1 ≠ ω2) / (1-Pvalue)"], label="", data, color=colors[dataset],
            linestyle=line_styles[method], xlabel="P(ω1 ≠ ω2) / (1-Pvalue)", ylabel="Differential sites")

    end
end

plot!([], [], line=:solid, linecolor=:black, label="difFUBAR")
plot!([], [], line=:dash, linecolor=:black, label="contrastFEL")
plot!([], [], line=:solid, linecolor=:green, label="hiv_envelope")
plot!([], [], line=:solid, linecolor=:blue, label="epidermal_leaf")
plot!([], [], line=:solid, linecolor=:red, label="rubisco")
plot!([], [], line=:solid, linecolor=:orange, label="hivRT")

plot!(legend=:topleft)


savefig("results/empirical_data/contrastfel_empirical/differential.png")
results
# Set x-axis limits
#xlims!(0, 0.2)
#ylims!(0, 100)


datasets_folder = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT_branchlength_1"]
datasets_name = ["epidermal_leaf", "hiv_envelope", "rubisco_C3_vs_C4", "hivRT"]

difFUBAR_results = []
contrastFEL_results = []
for i in 1:length(datasets_folder)
    dataset = datasets_folder[i]
    col_name = datasets_name[i]
    difFUBAR_location = "experiments/20240308_contrastFEL_empirical/results/$dataset/treesurgery_and_parallel/analysis_posteriors.csv"
    contrastFEL_location = "experiments/20240308_contrastFEL_empirical/results/$dataset/contrastfel/contrastfel.FEL.json"
    difFUBAR_res = read_in_difFUBAR_res(difFUBAR_location)
    contrastFEL_res = read_in_contrastFEL_res(contrastFEL_location)
    #contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]

    n_difFUBAR = nrow(difFUBAR_res[difFUBAR_res[!, "P(ω1 ≠ ω2)"].>0.75, :])
    n_contrastFEL = nrow(contrastFEL_res[contrastFEL_res[!, "P-value (overall)"].<0.05, :])

    push!(difFUBAR_results, n_difFUBAR)
    push!(contrastFEL_results, n_contrastFEL)
end

# threshold at difFUBAR 0.75 and p-value 0.05
res = DataFrame(dataset=datasets_name)
res[!, "difFUBAR_0.75"] = difFUBAR_results
res[!, "contrastFEL_0.95"] = contrastFEL_results
res[!, "purity"] = purity
# purity

res
output_path = "results/empirical_data/contrastfel_empirical/benchmark.csv"
CSV.write(output_path, res)



