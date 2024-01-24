cd("/home/ptruong/git/computationalPhylogenetics/")

using Pkg
using CSV
using JSON
using DataFrames
using Plots

function add_difference_cols_to_sim_settings(sim_settings, class1="simulator.omega.class0", class2="simulator.omega.class1")
    # sim_settings is pandas dfg
    diff_col = "difference.$class1" * ".vs." * "$class2"
    bool_diff_col = "bool.difference.$class1" * ".vs." * "$class2"
    sim_settings[!, diff_col] = sim_settings[!, class1] .- sim_settings[!, class2]
    sim_settings[!, bool_diff_col] = sim_settings[!, class1] .!= sim_settings[!, class2]
    return sim_settings, [diff_col, bool_diff_col]
end

# Simulator settings
simulator_settings_dir = "contrastFEL_data/omnibus/"
simulator_settings = CSV.File(simulator_settings_dir * "sims.1.settings.tsv", delim="\t", header=true) |> DataFrame
simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")
sum(simulator_settings[!, settings_cols[2]])

using DataFrames, CSV

function process_simulator_settings()
    results = DataFrame(filename=String[], sum_value=Float64[])
    simulator_settings_dir = "contrastFEL_data/omnibus/"  # Replace with your directory path

    for i in 1:500
        # Construct the file path
        file_path = joinpath(simulator_settings_dir, "sims.$i.settings.tsv")

        try
            # Load the CSV file into a DataFrame
            simulator_settings = CSV.File(file_path, delim="\t", header=true) |> DataFrame
            simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")

            # Calculate the sum for the desired column (adjust the column name as needed)
            sum_value = sum(simulator_settings[!, settings_cols[2]])

            # Append the results to the DataFrame
            push!(results, (filename=file_path, sum_value=sum_value))
        catch e
            # Handle any potential errors (e.g., file not found)
            println("Error processing file $file_path: $e")
        end
    end

    return results
end

# Call the function to process the simulator settings
simulator_results = process_simulator_settings()
simulator_results[!, "index"] = axis()


maximum(simulator_results[!, "sum_value"])
# Display the results DataFrame
println(simulator_results)
sort(simulator_results, "sum_value")[!, "filename"]