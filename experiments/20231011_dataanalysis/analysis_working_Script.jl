
cd("/home/ptruong/git/computationalPhylogenetics/experiments/20231011_run_both_pipelines_on_sample")
basedir = "/home/ptruong/git/computationalPhylogenetics/experiments/20231011_run_both_pipelines_on_sample/"

# null 
contrastFEL_null_data_dir = basedir * "/run_contrastFEL_null/run_contrastFEL_null/random/"
difFUBAR_null_data_dir = basedir * "/run_difFUBAR_null/output/null/pos_thresh_0.75/random-64-0.1-0.5/replicate.100/1v2/"


#cd(data_dir)
#readdir(data_dir; join=true)
using Plots
using Pkg
using CSV
using DataFrames
using JSON

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


########################
# diFFUBAR_null_data ###
########################


readdir(difFUBAR_null_data_dir, join=true)
# work on difFUBAR null
df = CSV.File(difFUBAR_null_data_dir * "results_posteriors.csv") |> DataFrame

findmax(df."P(ω1 > ω2)")
#df[313, "P(ω1 > ω2)"]

# ω1 ≠ ω2, don't how to do this propertly. but the one below should be two different version of the same thing.
#test, actually there should never be more than 400 trues, because we have only 400 sites
df = CSV.File(difFUBAR_null_data_dir * "results_posteriors.csv") |> DataFrame
data = sort(vcat(df."P(ω1 > ω2)", sort(df."P(ω1 > ω2)"))) # Is this really true though? we need to hard check this
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω1 ≠ ω2)", ylabel="False Positives Rate", label="Data") # This should be roughly equivalent to FPR?
#slope = length(data) 
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)


# Null set contains only false positives
# (ω1 > ω2)
df = CSV.File(difFUBAR_null_data_dir * "results_posteriors.csv") |> DataFrame
data = sort(df."P(ω1 > ω2)")
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω1 > ω2)", ylabel="False Positives Rate", label="Data") # This should be roughly equivalent to FPR?
#slope = length(data) 
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)

# P(ω2 > ω1)
df = CSV.File(difFUBAR_null_data_dir * "results_posteriors.csv") |> DataFrame
data = sort(df."P(ω2 > ω1)")
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω2 > ω1)", ylabel="False Positives Rate", label="Data")
#slope = length(data)
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)

# Am i wrong when i say p-value does not say anything on wether there is not difference, so TN is automatically incorprated. 



##########################
# difFUBAR omnibusMulti ##
##########################

contrastFEL_omnibusMulti_data_dir = basedir * "run_contrastFEL_omnibus_multi"
difFUBAR_omnibusMulti_data_dir = basedir * "/run_difFUBAR_omnibusMulti/output/omnibus_multi/pos_thresh_0.75/sim.0.replicate.1"
readdir(contrastFEL_omnibusMulti_data_dir, join=true)

df = CSV.File(difFUBAR_omnibusMulti_data_dir * "/1v2/results_posteriors.csv") |> DataFrame
# ω1 ≠ ω2, don't how to do this propertly. but the one below should be two different version of the same thing.
#test, actually there should never be more than 400 trues, because we have only 400 sites
df = CSV.File(difFUBAR_omnibusMulti_data_dir * "/1v2/results_posteriors.csv") |> DataFrame
data = sort(vcat(df."P(ω1 > ω2)", sort(df."P(ω1 > ω2)"))) # Is this really true though? we need to hard check this
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω1 ≠ ω2)", ylabel="What are the actual rates of this simulation?", label="Data") # This should be roughly equivalent to FPR?
#slope = length(data) 
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)


# Null set contains only false positives
# (ω1 > ω2)
df = CSV.File(difFUBAR_omnibusMulti_data_dir * "/1v2/results_posteriors.csv") |> DataFrame
data = sort(df."P(ω1 > ω2)")
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω1 > ω2)", ylabel="What are the actual rates of this simulation?", label="Data") # This should be roughly equivalent to FPR?
#slope = length(data) 
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)

# P(ω2 > ω1)
df = CSV.File(difFUBAR_omnibusMulti_data_dir * "/1v2/results_posteriors.csv") |> DataFrame
data = sort(df."P(ω2 > ω1)")
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω2 > ω1)", ylabel="What are the actual rates of this simulation?", label="Data")
#slope = length(data)
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)

#############################
# contrastFEL omnibusMulti ##
#############################

contrastFEL_omnibusMulti_data_dir = basedir * "run_contrastFEL_omnibus_multi/"
readdir(contrastFEL_omnibusMulti_data_dir, join=true)

file_name = contrastFEL_omnibusMulti_data_dir * "sims.0.settings.replicate.1.FEL.json"

df = read_contrastFEL_output_JSON(file_name)



