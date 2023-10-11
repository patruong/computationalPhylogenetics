
data_dir = "/home/ptruong/git/computationalPhylogenetics/experiments/20231009_pipeline_for_contrastFEL/output/pos_thresh_0.75/1v2"
cd(data_dir)
readdir(data_dir; join=true)

using Pkg
Pkg.add("CSV")
using CSV
using DataFrames

df = CSV.File("results_posteriors.csv") |> DataFrame
df

p_value_threshold = 0.05

df[df."P(ω1 > ω2)".<p_value_threshold, :]
df[df."P(ω2 > ω1)".<p_value_threshold, :]
df[df."P(ω1 > 1)".<p_value_threshold, :]
df[df."P(ω2 > 1)".<p_value_threshold, :]

df."P(ω1 > ω2)" .< p_value_threshold
df."P(ω2 > ω1)" .< p_value_threshold
df."P(ω1 > 1)" .< p_value_threshold
df."P(ω2 > 1)" .< p_value_threshold


# we do this over replicates... to get a distribution

# Do we for p-value loop like we did when benchmarking triqler
# Rank the df by pvalue and the and then plot, that should give us the curve.

using Plots
df = CSV.File("results_posteriors.csv") |> DataFrame

data = sort(df."P(ω1 > ω2)")

# Create a dataframe with one column named "value"
df = DataFrame(value=data)
insertcols!(df, :index => eachindex(data))
plot(df.value, df.index, xlabel="P(ω1 > ω2)", ylabel="Discovery", label="Data")
slope = length(data)
plot!(x -> slope * x, c=:red, line=:dash)


# we could also have this plot as percentage of data discovered.

# Something like this is good...
#
# Now can we run contrast-FEL 