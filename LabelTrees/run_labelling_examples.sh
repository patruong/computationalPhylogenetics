#!/bin/bash

# Define the directory path
output_directory="data/output/"

# Step 1: Remove the entire data/output/ directory and its contents
echo "Removing the entire $output_directory directory..."
rm -rf "$output_directory"

# Step 2: Create a new empty data/output/ directory
echo "Creating a new empty $output_directory directory..."
mkdir -p "$output_directory"

# Step 3: Run scripts containing "label_sequence"
for script in label_sequence_*.sh; do
    if [ -f "$script" ]; then
        echo "Running script: $script"
        bash "$script"
        echo "Script finished: $script"
        echo "------------------------"
    fi
done

# Additional step (optional): Run label-tree.bf
if [ -f "label-tree.bf" ]; then
    echo "Running script: label-tree.bf"
    hyphy label-tree.bf
    echo "label-tree.bf finished."
    echo "------------------------"
fi

