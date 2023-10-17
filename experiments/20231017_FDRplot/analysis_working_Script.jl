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
    sim_settings[!, "difference.$class1"*".vs."*"$class2"] = sim_settings[!, class1] .- sim_settings[!, class2]
    sim_settings[!, "bool.difference.$class1"*".vs."*"$class2"] = sim_settings[!, class1] .!= sim_settings[!, class2]
    return sim_settings
end

function calculate_ω1_not_equal_ω2(df)
    #The probability that ω1 and ω2 are not equal is the complement of the probability that one is greater than the other:
    #P(ω1 ≠ ω2) = 1 - [P(ω1 > ω2) + P(ω2 > ω1)]
    #This is because there are only two possibilities: either ω1 is greater than ω2, ω2 is greater than ω1, or they are equal. So, you subtract the sum of the probabilities that ω1 is greater than ω2 and ω2 is greater than ω1 from 1 to find the probability that they are not equal.
    df[!, "P(ω1 ≠ ω2)"] = 1 .- (difFUBAR_res[!, "P(ω1 > ω2)"] .+ difFUBAR_res[!, "P(ω2 > ω1)"])
    return df
end


# sim.0.replicate.1 

#difFUBAR_Res
difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_0.75/sim.0.replicate.1/REFERENCEvTEST/"
difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)


sort(difFUBAR_res)

sorted_df = sort(difFUBAR_res, order=(:"P(ω1 ≠ ω2)",))

typeof(difFUBAR_res)
sort(difFUBAR_res, "P(ω1 ≠ ω2)")

# contrastFEL use Q-value 
contrastFEL_output_dir = "contrastFEL_data/omnibus/"
contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "sims.0.settings.replicate.1.FEL.json")


# Simulator settings
simulator_settings_dir = "contrastFEL_data/omnibus/"
simulator_settings = CSV.File(simulator_settings_dir * "sims.0.settings.tsv", delim="\t", header=true) |> DataFrame
add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")



difFUBAR_res
contrastFEL_res
sim_settings

difFUBAR_res



