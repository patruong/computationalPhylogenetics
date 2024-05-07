
using Pkg
using Revise
using CSV
using JSON
using DataFrames
using Plots
using Statistics
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
function read_in_difFUBAR_res(pos_thresh, sim, rep)
    #difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST/"
    #difFUBAR_output_dir = "experiments/20231011_run_both_pipelines_on_sample/run_difFUBAR_omnibus/output/omnibus/pos_thresh_$pos_thresh/sim.$sim.replicate.$rep/REFERENCEvTEST"
    difFUBAR_output_dir = "/home/patrick/git/computationalPhylogenetics/experiments/20240116_rerun_difFUBAR/output_difFUBAR/sim_$sim" * "_replicate_$rep/patrick_max_child/"

    #difFUBAR_res = CSV.File(difFUBAR_output_dir * "results_posteriors.csv", delim=",", header=true) |> DataFrame
    difFUBAR_res = CSV.File(difFUBAR_output_dir * "_posteriors.csv", delim=",", header=true) |> DataFrame

    difFUBAR_res = calculate_ω1_not_equal_ω2(difFUBAR_res)
    id_fill(difFUBAR_res, sim, rep)
    return difFUBAR_res
end

# contrastFEL use Q-value 
function read_in_contrastFEL_res(sim, rep)
    contrastFEL_output_dir = "contrastFEL_data/omnibus/"
    contrastFEL_res = read_contrastFEL_output_JSON(contrastFEL_output_dir * "sims.$sim.settings.replicate.$rep.FEL.json")
    contrastFEL_res = add_codon_site(contrastFEL_res)
    id_fill(contrastFEL_res, sim, rep)
    return contrastFEL_res
end

# Simulator settings
function read_in_simulator_settings(sim)
    simulator_settings_dir = "contrastFEL_data/omnibus/"
    simulator_settings = CSV.File(simulator_settings_dir * "sims.$sim.settings.tsv", delim="\t", header=true) |> DataFrame
    simulator_settings, settings_cols = add_difference_cols_to_sim_settings(simulator_settings, "simulator.omega.class0", "simulator.omega.class1")
    #sum(simulator_settings[!, settings_cols[2]])
    return simulator_settings, settings_cols
end


function filter_on(df, col_name, lower_bound, upper_bound, keep_null)
    #predictions = lower_bound .< df[!, col_name] .< upper_bound
    filtered_df = filter(row -> lower_bound <= row[col_name] <= upper_bound, df)
    if keep_null == true
        null_df = filter(row -> row[col_name] == 0, df)
        filtered_df = vcat(filtered_df, null_df)
    end
    return filtered_df
end



function calculate_TPR_and_FPR_list(df_list, p_col)
    TPR_and_FPR_list = []
    for res in df_list
        res = calculate_TRP_and_FPR(res, p_col)
        push!(TPR_and_FPR_list, res)
    end
    return TPR_and_FPR_list
end

function calculate_ROC_values_threshold(df, threshold_col, threshold)
    diff = df[df[!, threshold_col].>=threshold, :] # predicted true, no diff
    no_diff = df[df[!, threshold_col].<threshold, :] # predicted false, no diff

    TP = sum(diff[!, "actual_difference"]) #TP
    FP = sum(diff[!, "actual_difference"] .== false) #FP
    FN = sum(no_diff[!, "actual_difference"]) #FN
    TN = sum(no_diff[!, "actual_difference"] .== false) #TNcontrastFEL_res = read_in_contrastFEL_res(sim, rep)


    TPR = TP / (TP + FN)
    FPR = FP / (FP + TN)
    FDR = FP / (FP + TP)
    FNR = FN / (TP + FN)
    return TP, FP, FN, TN, TPR, FPR, FDR, FNR
end

function calculate_ROC_threshold(df, threshold_col)
    # Create an empty DataFrame to store the results
    result_df = DataFrame(
        Threshold=Float64[],
        TP=Int[],
        FP=Int[],
        FN=Int[],
        TN=Int[],
        TPR=Float64[],
        FPR=Float64[],
        FDR=Float64[],
        FNR=Float64[]
    )

    threshold_slices = [-Inf, collect(0:0.001:1)..., Inf]
    threshold_slices .+= rand(length(threshold_slices)) .* 0.00000001

    for i in threshold_slices
        TP, FP, FN, TN, TPR, FPR, FDR, FNR = calculate_ROC_values_threshold(df, threshold_col, i)
        push!(result_df, (i, TP, FP, FN, TN, TPR, FPR, FDR, FNR))
    end

    return result_df
end

function merge_difFUBAR_contrastFEL_res(difFUBAR_res, contrastFEL_res)
    contrastFEL_res[!, "1-Pvalue"] = 1 .- contrastFEL_res[!, "P-value (overall)"]
    difFUBAR_res[!, "P(ω1 ≠ ω2)"]
    contrastFEL_res[!, "1-Pvalue"]
    difFUBAR_res[!, "Concatenated_ID"] = difFUBAR_res[!, "id"] .* "_codon.site." .* string.(difFUBAR_res[!, "Codon Sites"])
    contrastFEL_res[!, "Concatenated_ID"] = contrastFEL_res[!, "id"] .* "_codon.site." .* string.(contrastFEL_res[!, "Codon_site"])
    difFUBAR_res = rename(difFUBAR_res, names(difFUBAR_res) .=> ("difFUBAR_" .* string.(names(difFUBAR_res))))
    contrastFEL_res = rename(contrastFEL_res, names(contrastFEL_res) .=> ("contrastFEL_" .* string.(names(contrastFEL_res))))
    joined_res = innerjoin(difFUBAR_res, contrastFEL_res, on=:difFUBAR_Concatenated_ID => :contrastFEL_Concatenated_ID)

    if sum(joined_res[!, "difFUBAR_id"] .== joined_res[!, "contrastFEL_id"]) == nrow(joined_res)
        return joined_res
    else
        println("ERROR in matching row ids!")
        return false
    end
end


pos_thresh = 0.95
sims = collect(1:500)
reps = collect(1)

# Initialize variables to store aggregated results
aggregated_difFUBAR_res = DataFrame()
aggregated_contrastFEL_res = DataFrame()

for sim in sims
    for rep in reps
        try
            difFUBAR_res = read_in_difFUBAR_res(pos_thresh, sim, rep)
            contrastFEL_res = read_in_contrastFEL_res(sim, rep)
            simulator_settings, settings_cols = read_in_simulator_settings(sim)
            difFUBAR_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
            contrastFEL_res[!, "actual_difference"] = simulator_settings[!, settings_cols[2]]
            difFUBAR_res[!, "actual_effect_difference"] = abs.(simulator_settings[!, settings_cols[1]])
            contrastFEL_res[!, "actual_effect_difference"] = abs.(simulator_settings[!, settings_cols[1]])
            difFUBAR_res[!, "actual_directional_effect_difference"] = (simulator_settings[!, settings_cols[1]])
            contrastFEL_res[!, "actual_directional_effect_difference"] = (simulator_settings[!, settings_cols[1]])
            # Append the current batch of data to the aggregated DataFrames
            append!(aggregated_difFUBAR_res, difFUBAR_res)
            append!(aggregated_contrastFEL_res, contrastFEL_res)
        catch
            println("Error reading files for sim $sim and rep $rep. Skipping this iteration.")
        end
    end
end


difFUBAR_res = aggregated_difFUBAR_res
contrastFEL_res = aggregated_contrastFEL_res
joined_res = merge_difFUBAR_contrastFEL_res(difFUBAR_res, contrastFEL_res)

names(joined_res)
joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"]
joined_res[!, "contrastFEL_1-Pvalue"]

# Assuming joined_res is your DataFrame
#title="Scatter Plot of difFUBAR vs contrastFEL",

scatter(joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"], joined_res[!, "contrastFEL_P-value (overall)"],
    xlabel="difFUBAR P(ω1 ≠ ω2)", ylabel="contrastFEL P-value",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, joined_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)

#cd("/home/patrick/git/computationalPhylogenetics/results/hyphySim/coherence_between_difFUBAR_and_contrastFEL")

plot_dir = "/home/patrick/git/computationalPhylogenetics/results/hyphySim/coherence_between_difFUBAR_and_contrastFEL/"
savefig(plot_dir * "1_coherence_scatterplot.svg")


actual_true_res = filter(row -> row["difFUBAR_actual_difference"] != false, joined_res)

# Assuming joined_res is your DataFrame
#title="Scatter Plot of difFUBAR vs contrastFEL",

scatter(actual_true_res[!, "difFUBAR_P(ω1 ≠ ω2)"], actual_true_res[!, "contrastFEL_P-value (overall)"],
    xlabel="difFUBAR P(ω1 ≠ ω2)", ylabel="contrastFEL P-value",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, actual_true_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)
plot_dir = "/home/patrick/git/computationalPhylogenetics/results/hyphySim/coherence_between_difFUBAR_and_contrastFEL/"
savefig(plot_dir * "1_coherence_scatterplot_actual_true.svg")


actual_false_res = filter(row -> row["difFUBAR_actual_difference"] .== false, joined_res)

# Assuming joined_res is your DataFrame
#title="Scatter Plot of difFUBAR vs contrastFEL",

scatter(actual_false_res[!, "difFUBAR_P(ω1 ≠ ω2)"], actual_false_res[!, "contrastFEL_P-value (overall)"],
    xlabel="difFUBAR P(ω1 ≠ ω2)", ylabel="contrastFE P-value",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, actual_false_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)
plot_dir = "/home/patrick/git/computationalPhylogenetics/results/hyphySim/coherence_between_difFUBAR_and_contrastFEL/"
savefig(plot_dir * "1_coherence_scatterplot_actual_false.png")


#########################################
# Deeper analysis of contrastFEL issue  
# Filter 
# "difFUBAR_P(ω1 ≠ ω2)"].>0.9
# "contrastFEL_1-Pvalue"].<0.5
#########################################
joined_res
filtered_joined_res = joined_res[joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"].>0.99, :]
filtered_joined_res = filtered_joined_res[filtered_joined_res[!, "contrastFEL_1-Pvalue"].<0.01, :]

scatter(filtered_joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"], filtered_joined_res[!, "contrastFEL_1-Pvalue"],
    xlabel="difFUBAR_P(ω1 ≠ ω2)", ylabel="contrastFEL_1-Pvalue",
    title="P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, filtered_joined_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)
savefig(plot_dir * "2_coherence_scatterplot_difFUBAR_certain_contrastFEL_uncertain.png")
# High confidence from difFUBAR and low confidence from contrastFEL are actual hits that are correctly classified by difFUBAR, but missed for contrastFEL


csv_dir = plot_dir * "2_high_difFUBAR_certainty_low_contrastFEL_certainty_codon_sites_df.csv"
CSV.write(csv_dir, filtered_joined_res)


#######################################################
# investigate if same codon is faulty in every sim ####
#######################################################

filtered_joined_res[!, "difFUBAR_id"]

#unique_elements = unique(filtered_joined_res[!, "difFUBAR_id"])
# go through all sim on the filtered_joined_res

sims = collect(1:500)
reps = collect(1:5)

sim_dict = Dict()
for sim in sims
    sim_array = []
    for rep in reps
        simulation_name = "sim." * string(sim) * ".rep." * string(rep)
        if !isempty(filtered_joined_res[filtered_joined_res[!, "difFUBAR_id"].==simulation_name, :][!, "difFUBAR_Codon Sites"])
            push!(sim_array, filtered_joined_res[filtered_joined_res[!, "difFUBAR_id"].==simulation_name, :][!, "difFUBAR_Codon Sites"])
        end
    end
    #if sim_array != Any[]
    if length(sim_array) > 1
        sim_dict[sim] = sim_array
    end
    #end
end


sim_dict






sims = []
values = []
for (key, value) in sim_dict
    push!(sims, key)
    push!(values, value)
end

high_difFUBAR_certainty_low_contrastFEL = DataFrame(sim=sims, codon_sites=values)
csv_dir = plot_dir * "2_high_difFUBAR_certainty_low_contrastFEL_certainty_codon_sites.csv"
CSV.write(csv_dir, high_difFUBAR_certainty_low_contrastFEL)


sim_dict

recurring_codon_errors = Dict()
for (key, value) in sim_dict
    if intersect(value...) != Int64[]
        recurring_codon_errors[key] = intersect(value...)
    end
end


###############################################
# Check these manually one-by-one
#
# THESE ERROR ARE CONSISTENT AMONG REPLICATES #
recurring_codon_errors #list of simulations with codon sites that are faulty in many simulations.

sims = []
values = []
for (key, value) in recurring_codon_errors
    push!(sims, key)
    push!(values, value)
end
recurring_faulty_sites = DataFrame(sim=sims, codon_sites=values)
csv_dir = plot_dir * "2_recurring_error_codon_sites.csv"
CSV.write(csv_dir, recurring_faulty_sites)



joined_res


function select_sim_and_codon_site(df, sim_n, codon_sites)
    sim_n = string(sim_n)
    target_codon_sites = codon_sites
    target_ids = ["sim.$sim_n.rep.1", "sim.$sim_n.rep.2", "sim.$sim_n.rep.3", "sim.$sim_n.rep.4", "sim.$sim_n.rep.5"]
    filtered_joined_res = filter(row -> row["difFUBAR_id"] ∈ target_ids && row["difFUBAR_Codon Sites"] ∈ target_codon_sites, df)
    return filtered_joined_res
end

x1 = select_sim_and_codon_site(joined_res, 83, [38, 87, 206, 591, 658, 758])
x2 = select_sim_and_codon_site(joined_res, 56, [705])  # only 3 samples
x3 = select_sim_and_codon_site(joined_res, 175, [634])  # only 3 samples
x4 = select_sim_and_codon_site(joined_res, 26, [324])
x5 = select_sim_and_codon_site(joined_res, 347, [3, 315, 672])
recurrent_error_df = vcat(x1, x2, x3, x4, x5)

x1 = select_sim_and_codon_site(joined_res, 56, [705]) # only 3 samples
x2 = select_sim_and_codon_site(joined_res, 175, [634])  # only 3 samples
x3 = select_sim_and_codon_site(joined_res, 60, [238])
x4 = select_sim_and_codon_site(joined_res, 26, [324])
recurrent_error_df = vcat(x1, x2, x3, x4)


csv_dir = plot_dir * "2_recurring_error_codon_sites_df.csv"
CSV.write(csv_dir, recurrent_error_df)


# Fix the check of codon_sites  in bens chat 
names(diagnose_df)



function df_col_prefix(df, prefix)
    selected_columns_values = Dict{String,Vector}()
    for col_name in names(df)
        if occursin(prefix, col_name)
            selected_columns_values[col_name] = df[!, col_name]
        end
    end
    selected_df = DataFrame(selected_columns_values)
    return selected_df
end

df_col_prefix(joined_res, "contrastFEL")
df_col_prefix(joined_res, "alpha") # ERROR
df_col_prefix(joined_res, "beta") # ERROR
df_col_prefix(diagnose_df, "contrastFEL")
df_col_prefix(diagnose_df, "alpha") # ERROR 
df_col_prefix(diagnose_df, "beta") # ERRPR
names(df_col_prefix(diagnose_df, "contrastFEL"))

names(df_col_prefix(diagnose_df, "difFUBAR"))

df_col_prefix(diagnose_df, "difFUBAR")
df_col_prefix(diagnose_df, "difFUBAR_mean")
df_col_prefix(joined_res, "difFUBAR_mean(α)")



# The contrastFEL errors when we have errors in setting alpha and beta
df_col_prefix(filtered_joined_res, "alpha")


####################################
# test the code for sim 82, codon 133 - alpha and beta error
# sim 82, codon 601 - beta error
# sim 82, codon 604 - beta error




################################
# plot up distribution of beta #
################################
df_col_prefix(joined_res, "alpha")
df_col_prefix(joined_res, "beta")

x1 = df_col_prefix(joined_res, "alpha")[!, "contrastFEL_alpha"]
histogram(x1, label="Alpha (no filter)", alpha=0.5, density=true)
xlabel!("Alpha")
ylabel!("Count")
title!("contrastFEL alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_contrastFEL_alpha.png")

x2 = df_col_prefix(joined_res, "beta")[!, "contrastFEL_beta (TEST)"]
histogram(x2, label="Beta (test) (no filter)", alpha=0.5, density=true)
xlabel!("Beta (test)")
ylabel!("Count")
title!("contrastFEL beta (test)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_contrastFEL_beta_test.png")


x3 = df_col_prefix(joined_res, "beta")[!, "contrastFEL_beta (background)"]
histogram(x3, label="Beta (background) (no filter)", alpha=0.5, density=true)
xlabel!("Beta (background)")
ylabel!("Count")
title!("contrastFEL beta (background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_contrastFEL_beta_background.png")


x = joined_res[:, ["contrastFEL_alpha", "contrastFEL_beta (TEST)", "contrastFEL_beta (background)"]]
x4 = maximum.(eachrow(x))

histogram(x4, label="max(alpha, beta(test), beta(background)) (no filter)", alpha=0.5, density=true)
xlabel!("Max of alpha, beta(test), beta(background)")
ylabel!("Count")
title!("Max of alpha, beta(test), beta(background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_contrastFEL_max_of_alpha_betas.png")



names(difFUBAR_res)
x5 = df_col_prefix(joined_res, "mean(α)")[!, "difFUBAR_mean(α)"]
histogram(x5, label="Alpha (no filter)", alpha=0.5, density=true)
xlabel!("mean(α)")
ylabel!("Count")
title!("difFUBAR alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_difFUBAR_alpha.png")


x6 = df_col_prefix(joined_res, "mean(ω1)")[!, "difFUBAR_mean(ω1)"]
histogram(x6, label="mean(ω1) (no filter)", alpha=0.5, density=true)
xlabel!("mean(ω1)")
ylabel!("Count")
title!("difFUBAR ω1")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_difFUBAR_omega_1.png")


x7 = df_col_prefix(joined_res, "mean(ω2)")[!, "difFUBAR_mean(ω2)"]
histogram(x7, label="mean(ω2) (no filter)", alpha=0.5, density=true)
xlabel!("mean(ω2)")
ylabel!("Count")
title!("difFUBAR ω2")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "1_histogram_all_difFUBAR_omega_2.png")


#######################################
# Plot of filtered data (faulty data) #
#######################################

df_col_prefix(filtered_joined_res, "Pvalue")

maximum(df_col_prefix(filtered_joined_res, "Pvalue")[!, "contrastFEL_1-Pvalue"])

df_col_prefix(filtered_joined_res, "alpha")
df_col_prefix(filtered_joined_res, "beta")

x1 = df_col_prefix(filtered_joined_res, "alpha")[!, "contrastFEL_alpha"]
histogram(x1, label="Alpha P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true, nbins=50)
xlabel!("Alpha")
ylabel!("Count")
title!("contrastFEL alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_contrastFEL_alpha.png")


x2 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (TEST)"]
histogram(x2, label="Beta(test) P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (test)")
ylabel!("Count")
title!("contrastFEL beta (test)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_contrastFEL_beta_test.png")


x3 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (background)"]
histogram(x3, label="Beta (background) P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (background)")
ylabel!("Count")
title!("contrastFEL beta (background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_contrastFEL_beta_background.png")


# Make analysis of what max alpha, beta (background), beta (test)

x = filtered_joined_res[:, ["contrastFEL_alpha", "contrastFEL_beta (TEST)", "contrastFEL_beta (background)"]]
x4 = maximum.(eachrow(x))

histogram(x4, label="max(alpha, beta(test), beta(background)) P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true, nbins=50)
xlabel!("Max of alpha, beta(test), beta(background)")
ylabel!("Count")
title!("Max of alpha, beta(test), beta(background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_contrastFEL_max_of_alpha_betas.png")


x5 = df_col_prefix(filtered_joined_res, "mean(α)")[!, "difFUBAR_mean(α)"]
histogram(x5, label="Alpha (no filter) P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true)
xlabel!("mean(α)")
ylabel!("Count")
title!("difFUBAR alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_difFUBAR_alpha.png")


x6 = df_col_prefix(filtered_joined_res, "mean(ω1)")[!, "difFUBAR_mean(ω1)"]
histogram(x6, label="mean(ω1) (no filter) P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true)
xlabel!("mean(ω1)")
ylabel!("Count")
title!("difFUBAR ω1")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_difFUBAR_omega_1.png")


x7 = df_col_prefix(filtered_joined_res, "mean(ω2)")[!, "difFUBAR_mean(ω2)"]
histogram(x7, label="mean(ω2) (no filter) P(ω1 ≠ ω2) > 0.99 / contrastFEL_1-Pvalue > 0.01", alpha=0.5, density=true)
xlabel!("mean(ω2)")
ylabel!("Count")
title!("difFUBAR ω2")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "2_histogram_filtered_bottom_right_difFUBAR_omega_2.png")

######################################
# Deeper analysis of difFUBAR issues #
# Filter 
# "difFUBAR_P(ω1 ≠ ω2)"].<0.5
# "contrastFEL_1-Pvalue"].>0.9
######################################

filtered_joined_res = joined_res[joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"].<0.50, :]
filtered_joined_res = filtered_joined_res[filtered_joined_res[!, "contrastFEL_1-Pvalue"].>0.99, :]

scatter(filtered_joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"], filtered_joined_res[!, "contrastFEL_1-Pvalue"],
    xlabel="difFUBAR_P(ω1 ≠ ω2)", ylabel="contrastFEL_1-Pvalue",
    title="P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, filtered_joined_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)
savefig(plot_dir * "3_coherence_scatterplot_difFUBAR_uncertain_contrastFEL_certain.png")
# Low confidence from difFUBAR and high confidence from contrastFEL are actual hits that are incorrectly classified by contrastFEL and correctly called wrong by difFUBAR


x1 = df_col_prefix(filtered_joined_res, "alpha")[!, "contrastFEL_alpha"]
histogram(x1, label="Alpha P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true, nbins=50)
xlabel!("Alpha")
ylabel!("Count")
title!("contrastFEL alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_contrastFEL_alpha.png")


x2 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (TEST)"]
histogram(x2, label="Beta(test) P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (test)")
ylabel!("Count")
title!("contrastFEL beta (test)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_contrastFEL_beta_test.png")


x3 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (background)"]
histogram(x3, label="Beta (background) P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (background)")
ylabel!("Count")
title!("contrastFEL beta (background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_contrastFEL_beta_background.png")


# Make analysis of what max alpha, beta (background), beta (test)

x = filtered_joined_res[:, ["contrastFEL_alpha", "contrastFEL_beta (TEST)", "contrastFEL_beta (background)"]]
x4 = maximum.(eachrow(x))

histogram(x4, label="max(alpha, beta(test), beta(background)) P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true, nbins=50)
xlabel!("Max of alpha, beta(test), beta(background)")
ylabel!("Count")
title!("Max of alpha, beta(test), beta(background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_contrastFEL_max_of_alpha_betas.png")



x5 = df_col_prefix(filtered_joined_res, "mean(α)")[!, "difFUBAR_mean(α)"]
histogram(x5, label="Alpha (no filter) P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true)
xlabel!("mean(α)")
ylabel!("Count")
title!("difFUBAR alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_difFUBAR_alpha.png")


x6 = df_col_prefix(filtered_joined_res, "mean(ω1)")[!, "difFUBAR_mean(ω1)"]
histogram(x6, label="mean(ω1) (no filter) P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true)
xlabel!("mean(ω1)")
ylabel!("Count")
title!("difFUBAR ω1")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_difFUBAR_omega_1.png")


x7 = df_col_prefix(filtered_joined_res, "mean(ω2)")[!, "difFUBAR_mean(ω2)"]
histogram(x7, label="mean(ω2) (no filter) P(ω1 ≠ ω2) < 0.5 / contrastFEL_1-Pvalue > 0.99", alpha=0.5, density=true)
xlabel!("mean(ω2)")
ylabel!("Count")
title!("difFUBAR ω2")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "3_histogram_filtered_top_difFUBAR_omega_2.png")


#####################
# check right axis ##
#####################

filtered_joined_res = joined_res[joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"].>0.999, :]

scatter(filtered_joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"], filtered_joined_res[!, "contrastFEL_1-Pvalue"],
    xlabel="difFUBAR_P(ω1 ≠ ω2)", ylabel="contrastFEL_1-Pvalue",
    title="P(ω1 ≠ ω2) > 0.999",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, filtered_joined_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)
savefig(plot_dir * "4_coherence_high_difFUBAR_certainty_all_contrastFEL.png")


x1 = df_col_prefix(filtered_joined_res, "alpha")[!, "contrastFEL_alpha"]
histogram(x1, label="Alpha P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true, nbins=50)
xlabel!("Alpha")
ylabel!("Count")
title!("contrastFEL alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_contrastFEL_alpha.png")


x2 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (TEST)"]
histogram(x2, label="Beta(test) P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (test)")
ylabel!("Count")
title!("contrastFEL beta (test)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_contrastFEL_beta_test.png")


x3 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (background)"]
histogram(x3, label="Beta (background) P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (background)")
ylabel!("Count")
title!("contrastFEL beta (background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_contrastFEL_beta_background.png")


# Make analysis of what max alpha, beta (background), beta (test)

x = filtered_joined_res[:, ["contrastFEL_alpha", "contrastFEL_beta (TEST)", "contrastFEL_beta (background)"]]
x4 = maximum.(eachrow(x))

histogram(x4, label="max(alpha, beta(test), beta(background)) P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true, nbins=50)
xlabel!("Max of alpha, beta(test), beta(background)")
ylabel!("Count")
title!("Max of alpha, beta(test), beta(background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_contrastFEL_max_of_alpha_betas.png")



x5 = df_col_prefix(filtered_joined_res, "mean(α)")[!, "difFUBAR_mean(α)"]
histogram(x5, label="Alpha (no filter) P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true)
xlabel!("mean(α)")
ylabel!("Count")
title!("difFUBAR alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_difFUBAR_alpha.png")


x6 = df_col_prefix(filtered_joined_res, "mean(ω1)")[!, "difFUBAR_mean(ω1)"]
histogram(x6, label="mean(ω1) (no filter) P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true)
xlabel!("mean(ω1)")
ylabel!("Count")
title!("difFUBAR ω1")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_difFUBAR_omega_1.png")


x7 = df_col_prefix(filtered_joined_res, "mean(ω2)")[!, "difFUBAR_mean(ω2)"]
histogram(x7, label="mean(ω2) (no filter) P(ω1 ≠ ω2) > 0.999", alpha=0.5, density=true)
xlabel!("mean(ω2)")
ylabel!("Count")
title!("difFUBAR ω2")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "4_histogram_filtered_right_difFUBAR_omega_2.png")

#####################
# check bottom axis ##
#####################

filtered_joined_res = joined_res[joined_res[!, "contrastFEL_1-Pvalue"].<0.001, :]

scatter(filtered_joined_res[!, "difFUBAR_P(ω1 ≠ ω2)"], filtered_joined_res[!, "contrastFEL_1-Pvalue"],
    xlabel="difFUBAR_P(ω1 ≠ ω2)", ylabel="contrastFEL_1-Pvalue",
    title="1-Pvalue < 0.001",
    legend=false, alpha=0.02, size=(800, 600),
    color=map(x -> x == 1 ? :blue : :red, filtered_joined_res[!, "contrastFEL_actual_difference"])
    #color=[:blue, :red][joined_res[!, "difFUBAR_actual_difference"].+1])
)
savefig(plot_dir * "5_coherence_all_difFUBAR_certainty_low_contrastFEL.png")


x1 = df_col_prefix(filtered_joined_res, "alpha")[!, "contrastFEL_alpha"]
histogram(x1, label="Alpha 1-Pvalue < 0.001", alpha=0.5, density=true, nbins=50)
xlabel!("Alpha")
ylabel!("Count")
title!("contrastFEL alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_contrastFEL_alpha.png")


x2 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (TEST)"]
histogram(x2, label="Beta(test) 1-Pvalue < 0.001", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (test)")
ylabel!("Count")
title!("contrastFEL beta (test)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_contrastFEL_beta_test.png")


x3 = df_col_prefix(filtered_joined_res, "beta")[!, "contrastFEL_beta (background)"]
histogram(x3, label="Beta (background) 1-Pvalue < 0.001", alpha=0.5, density=true, nbins=50)
xlabel!("Beta (background)")
ylabel!("Count")
title!("contrastFEL beta (background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_contrastFEL_beta_background.png")


# Make analysis of what max alpha, beta (background), beta (test)

x = filtered_joined_res[:, ["contrastFEL_alpha", "contrastFEL_beta (TEST)", "contrastFEL_beta (background)"]]
x4 = maximum.(eachrow(x))

histogram(x4, label="max(alpha, beta(test), beta(background)) 1-Pvalue < 0.001", alpha=0.5, density=true, nbins=50)
xlabel!("Max of alpha, beta(test), beta(background)")
ylabel!("Count")
title!("Max of alpha, beta(test), beta(background)")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_contrastFEL_max_of_alpha_betas.png")



x5 = df_col_prefix(filtered_joined_res, "mean(α)")[!, "difFUBAR_mean(α)"]
histogram(x5, label="Alpha (no filter) 1-Pvalue < 0.001", alpha=0.5, density=true)
xlabel!("mean(α)")
ylabel!("Count")
title!("difFUBAR alpha")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_difFUBAR_alpha.png")


x6 = df_col_prefix(filtered_joined_res, "mean(ω1)")[!, "difFUBAR_mean(ω1)"]
histogram(x6, label="mean(ω1) (no filter) 1-Pvalue < 0.001", alpha=0.5, density=true)
xlabel!("mean(ω1)")
ylabel!("Count")
title!("difFUBAR ω1")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_difFUBAR_omega_1.png")


x7 = df_col_prefix(filtered_joined_res, "mean(ω2)")[!, "difFUBAR_mean(ω2)"]
histogram(x7, label="mean(ω2) (no filter) 1-Pvalue < 0.001", alpha=0.5, density=true)
xlabel!("mean(ω2)")
ylabel!("Count")
title!("difFUBAR ω2")
plot!(legend=:topleft)  # Position the legend
savefig(plot_dir * "5_histogram_filtered_bottom_difFUBAR_omega_2.png")


###############################################
# How can contrastFEL be so certain on these?
# Do we need to know this as well ?  #
###############################################

filtered_joined_res[:, ["contrastFEL_actual_directional_effect_difference", "contrastFEL_alpha", "contrastFEL_beta (TEST)", "contrastFEL_beta (background)"]]

names(filtered_joined_res)

# This should be fine no?

