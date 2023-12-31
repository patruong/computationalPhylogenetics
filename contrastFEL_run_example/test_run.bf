# Load your data
msa = ReadDataFile("/home/ptruong/git/computationalPhylogenetics/contrastFEL_run_example/flu_sim_seq_and_bm.fasta")
tree = ReadTree("/home/ptruong/git/computationalPhylogenetics/contrastFEL_run_example/flu_sim.tre")

# Run the contrast-fel analysis
contrasts = ContrastFEL(msa, tree)

# Specify the partitions and models if needed
# partitions = [1, 2, 3] # List of partitions
# models = ["MG94", "GY94"] # List of models

# Fit models to data
# FitModel(contrasts, partitions, models)

# Execute the analysis
ExecuteAFile(contrasts)
