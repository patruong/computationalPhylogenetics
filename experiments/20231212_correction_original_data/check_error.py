import os 
import pandas as pd
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20231212_correction_original_data/output_full_tree/")

baseline = pd.read_csv("baseline/_con_lik_matrix.csv", sep = ",", header = None)
final = pd.read_csv("final/_con_lik_matrix.csv", sep = ",", header = None)

# Check error
(final - baseline).sum().sum() #check total error
final - baseline # error 
(final - baseline).max().max()
(final - baseline).min().min()
final == baseline

# time hard-coded speedup
baseline_time = 535*60 + 40
final_time = 32*60 + 32
baseline_time / final_time
# contrast-FEL runs on pruned nwk tree with 320 sequences.
run_1 = 953*60 + 42
run_2 = 873*60 + 7

run_1 / final_time
run_2 / final_time



