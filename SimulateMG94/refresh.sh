#!/bin/bash

# List of directory paths and filenames
directories=("data/example1" "data/example2" "data/example3" "data/example4")
filenames=("example1" "example2" "example3" "example4")

# Function to remove existing files and directories
remove_existing() {
    for dir in "${directories[@]}"; do
        if [ -d "$dir" ]; then
            echo "Removing directory: $dir"
            rm -rf "$dir"
        elif [ -f "$dir" ]; then
            echo "Removing file: $dir"
            rm -f "$dir"
        fi
    done
}

# Function to create new empty directories and files
create_directories_and_files() {
    for dir in "${directories[@]}"; do
        echo "Creating directory: $dir"
        mkdir -p "$dir"
    done
    
    for ((i=0; i<${#directories[@]}; i++)); do
        dir="${directories[i]}"
        file="${filenames[i]}"
        echo "Creating empty file: $dir/$file"
        touch "$dir/$file"
    done
}

# Main script
echo "Refreshing directory files..."

# Remove existing files and directories
remove_existing

# Create new empty directories and files
create_directories_and_files

echo "Refresh complete."
