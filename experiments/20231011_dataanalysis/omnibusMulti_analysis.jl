using Plots
using Pkg
using CSV
using DataFrames
using JSON


####################
# contrast-FEL part#
####################

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

basedir = "/home/ptruong/git/computationalPhylogenetics/experiments/20231011_run_both_pipelines_on_sample/"
difFUBAR_omnibusMulti_data_dir = basedir * "/run_difFUBAR_omnibusMulti/output/omnibus_multi/pos_thresh_0.75/sim.0.replicate.1"
contrastFEL_omnibusMulti_data_dir = basedir * "run_contrastFEL_omnibus_multi/"

#

# Read in contrastFEL
file_name = contrastFEL_omnibusMulti_data_dir * "sims.0.settings.replicate.1.FEL.json"
df = read_contrastFEL_output_JSON(file_name)
CSV.write("output_file.csv", df)


pairwise_comparisons = []
for i in names(df)
    if occursin(r"^P-value", i)
        # Extract "GROUP_X" and "GROUP_Y" from the column name
        match_result = match(r"P-value for (\w+) vs (\w+)", i)

        if match_result !== nothing
            group_x = match_result.captures[1]
            group_y = match_result.captures[2]
            push!(pairwise_comparisons, [group_x, group_y])
        end
    end
end


group1 = pairwise_comparisons[3][1]
group2 = pairwise_comparisons[3][2]

group1 = "GROUP_3"
group2 = "GROUP_2"

df = read_contrastFEL_output_JSON(file_name)

CSV.write("output_file.csv", df)

data = sort(df[!, "P-value for $group1 vs $group2"])
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)

plot(df.value, df.index, xlabel="P(ω1 ≠ ω2)", ylabel="What are the actual rates of this simulation?", label="Data") # This should be roughly equivalent to FPR?
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)



################
# difFUBAR part#
################

comparisons_pairs = []
for s in readdir(difFUBAR_omnibusMulti_data_dir, join=true)
    groups = split(s, "v") # split s on "v" and get a list of two elements
    group1 = groups[1][end-6:end] # get the first group name by slicing the string from the end
    group2 = groups[2] # get the second group name
    push!(comparisons_pairs, [group1, group2]) # print the group names
end


group1 = comparisons_pairs[1][1]
group2 = comparisons_pairs[1][2]

df = CSV.File(difFUBAR_omnibusMulti_data_dir * "/$group1" * "v" * "$group2/results_posteriors.csv") |> DataFrame
data = sort(vcat(df."P(ω1 > ω2)", sort(df."P(ω2 > ω1)"))) # Is this really true though? we need to hard check this
# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
df.index = df.index ./ length(data)
plot(df.value, df.index, xlabel="P(ω1 ≠ ω2)", ylabel="What are the actual rates of this simulation?", label="Data") # This should be roughly equivalent to FPR?
#slope = length(data) 
slope = 1
plot!(x -> slope * x, c=:red, line=:dash)
