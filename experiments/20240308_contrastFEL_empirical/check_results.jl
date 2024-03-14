
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


dir = "experiments/20240308_contrastFEL_empirical/output/rubisco_C3_vs_C4/contrastfel/"
df = read_contrastFEL_output_JSON(dir * "contrastfel.FEL.json")
filter(row -> row."P-value (overall)" < 0.05, df) # this is ok as in the paper 15 codons are significant

dir = "experiments/20240308_contrastFEL_empirical/output/epidermal_leaf/contrastfel/"
df = read_contrastFEL_output_JSON(dir * "contrastfel.FEL.json")
filter(row -> row."P-value (overall)" < 0.05, df) # the oaoer has 28 codons, and I have 29 codons. This is ok!


dir = "experiments/20240308_contrastFEL_empirical/output/hiv_envelope/contrastfel/"
df = read_contrastFEL_output_JSON(dir * "contrastfel.FEL.json")
filter(row -> row."P-value (overall)" < 0.05, df) # the paper has 31 codons, and I have 39 codons. This is ok? 




