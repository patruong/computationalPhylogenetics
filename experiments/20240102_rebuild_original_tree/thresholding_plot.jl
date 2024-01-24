cd("/home/patrick/git/computationalPhylogenetics/experiments/20240102_rebuild_original_tree/")

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

function read_in_difFUBAR_res(pos_thresh, sim, rep)
    difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST/"
    difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
    difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)
    id_fill(difFUBAR_res, sim, rep)
    return difFUBAR_res
end


function calculate_ω1_not_equal_ω2(df)
    #The probability that ω1 and ω2 are not equal is the complement of the probability that one is greater than the other:
    #P(ω1 ≠ ω2) = 1 - [P(ω1 > ω2) + P(ω2 > ω1)]
    #This is because there are only two possibilities: either ω1 is greater than ω2, ω2 is greater than ω1, or they are equal. So, you subtract the sum of the probabilities that ω1 is greater than ω2 and ω2 is greater than ω1 from 1 to find the probability that they are not equal.
    df[!, "P(ω1 ≠ ω2)"] = (df[!, "P(ω1 > ω2)"] .+ df[!, "P(ω2 > ω1)"])
    return df
end

contrastFEL_output_dir = "output/contrastfel/"
contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "contrastfel.FEL.json")

difFUBAR_output_dir = "output/baseline/"
difFUBAR_res = CSV.File(difFUBAR_output_dir * "_posteriors.csv", delim=",", header=true) |> DataFrame
difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)

# plot this x-axis threshold, y axis number significant
contrastFEL_res[contrastFEL_res[!, "P-value (overall)"].<0.05, :]
contrastFEL_res[contrastFEL_res[!, "Permutation p-value"].<0.05, :]
contrastFEL_res[contrastFEL_res[!, "Q-value (overall)"].<0.05, :]
difFUBAR_res[difFUBAR_res[!, "P(ω1 ≠ ω2)"].>0.95, :]


using Plots

# Assuming contrastFEL_res and difFUBAR_res are your DataFrames
# Define a range of threshold values
threshold_values = 0.00:0.001:0.1

# Initialize arrays to store the number of significant hits for each threshold
hits_p_value = zeros(Int, length(threshold_values))
hits_permutation_p_value = zeros(Int, length(threshold_values))
hits_q_value = zeros(Int, length(threshold_values))
hits_p_omega = zeros(Int, length(threshold_values))

# Iterate over each threshold value
for (i, threshold) in enumerate(threshold_values)
    # Count the number of significant hits for each condition
    hits_p_value[i] = sum(contrastFEL_res[!, "P-value (overall)"] .< threshold)
    hits_permutation_p_value[i] = sum(contrastFEL_res[!, "Permutation p-value"] .< threshold)
    hits_q_value[i] = sum(contrastFEL_res[!, "Q-value (overall)"] .< threshold)
    #hits_p_omega[i] = sum(difFUBAR_res[!, "P(ω1 ≠ ω2)"] .> 1 - threshold)
    hits_p_omega[i] = sum((1 .- difFUBAR_res[!, "P(ω1 ≠ ω2)"]) .< threshold)

end



# Plot the results
plot(threshold_values, hits_p_value, label="P-value (overall)", xlabel="P-value / Q-value / 1-Posterior", ylabel="Number of Significant Hits", legend=:topleft)
plot!(threshold_values, hits_permutation_p_value, label="Permutation p-value")
plot!(threshold_values, hits_q_value, label="Q-value (overall)")
plot!(threshold_values, hits_p_omega, label="1 - P(ω1 ≠ ω2)")
