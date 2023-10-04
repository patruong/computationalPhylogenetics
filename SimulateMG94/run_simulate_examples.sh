#!/bin/bash

# Get the current directory
script_directory=$(pwd)

# Loop through the scripts in the directory and execute them
for script in "$script_directory"/simulate*.sh; do
    if [ -f "$script" ]; then
        echo "Running script: $script"
        bash "$script"
        echo "Script finished: $script"
        echo "------------------------"
    fi
done
