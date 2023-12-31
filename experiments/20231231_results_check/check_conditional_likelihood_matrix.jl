using DataFrames
using CSV

cd("/home/patrick/git/computationalPhylogenetics/experiments/20231231_results_check/difFUBAR_output_orig/")
file_baseline = "output/baseline/_con_lik_matrix.csv"
file_max = "output/max/_con_lik_matrix.csv"
file_patrick = "output/patrick/_con_lik_matrix.csv"
file_patrick_max = "output/patrick_max/_con_lik_matrix.csv"
file_patrick_max_child = "output/patrick_max_child/_con_lik_matrix.csv"
file_final = "output/final/_con_lik_matrix.csv"

function try_read_df(file)
    try
        return CSV.File(file, header=false) |> DataFrame
    catch
        return missing
    end
end


con1 = try_read_df(file_baseline)
con2 = try_read_df(file_max)
con3 = try_read_df(file_patrick)
con4 = try_read_df(file_patrick_max)
con5 = try_read_df(file_patrick_max_child)
con6 = try_read_df(file_final)




function df_comparison(df1, df2)
    try
        return df1 == df2
    catch
        return missing
    end
end

function df_max_residual(df1, df2)
    try
        error = abs.(df1 .- df2)
        return maximum(maximum.(eachcol(error)))
    catch
        return missing
    end
end

function df_sum_residual(df1, df2)
    try
        error = abs.(df1 .- df2)
        return sum(sum.(eachcol(error)))
    catch
        return missing
    end
end

df_comparison(con1, con2)
df_comparison(con1, con3)
df_comparison(con1, con4)
df_comparison(con1, con5)
df_comparison(con1, con6)


df_max_residual(con1, con2)
df_max_residual(con1, con3)
df_max_residual(con1, con4)
df_max_residual(con1, con5)
df_max_residual(con1, con6)

df_sum_residual(con1, con2)
df_sum_residual(con1, con3)
df_sum_residual(con1, con4)
df_sum_residual(con1, con5)
df_sum_residual(con1, con6)


# Create a DataFrame with summary results
summary_df = DataFrame(
    Comparison=["con1 == con2", "con1 == con3", "con1 == con4", "con1 == con5", "con1 == con6"],
    Comparison_res=[df_comparison(con1, con2),
        df_comparison(con1, con3),
        df_comparison(con1, con4),
        df_comparison(con1, con5),
        df_comparison(con1, con6)],
    Max_Error=[
        df_max_residual(con1, con2),
        df_max_residual(con1, con3),
        df_max_residual(con1, con4),
        df_max_residual(con1, con5),
        df_max_residual(con1, con6)
    ],
    Total_Error=[df_sum_residual(con1, con2),
        df_sum_residual(con1, con3),
        df_sum_residual(con1, con4),
        df_sum_residual(con1, con5),
        df_sum_residual(con1, con6)]
)

# Display the resulting DataFrame
println(summary_df)


