using DataFrames
using CSV

cd("/home/patrick/git/computationalPhylogenetics/experiments/20231222_rebuild_parvoVP_on_difFUBAR")
file_baseline = "output/baseline/_con_lik_matrix.csv"
file_max = "output/max/_con_lik_matrix.csv"
file_patrick = "output/patrick/_con_lik_matrix.csv"
file_patrick_max = "output/patrick_max/_con_lik_matrix.csv"
file_final = "output/final/_con_lik_matrix.csv"

con1 = CSV.File(file_baseline) |> DataFrame
con2 = CSV.File(file_max) |> DataFrame
con3 = CSV.File(file_patrick) |> DataFrame
con4 = CSV.File(file_patrick_max) |> DataFrame
con5 = CSV.File(file_final) |> DataFrame


con1 == con2
con1 == con3
con1 == con4
con1 == con5












