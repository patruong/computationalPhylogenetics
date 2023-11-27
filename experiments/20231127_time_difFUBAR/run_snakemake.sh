#!/bin/bash

# Define the command you want to run
SNAKEMAKE_COMMAND="snakemake -c39 --rerun-incomplete"

# Start an infinite loop
while true; do
    # Run the Snakemake command
    $SNAKEMAKE_COMMAND

    # Check the exit status of the command
    if [ $? -eq 0 ]; then
        # If the command exited successfully, break out of the loop
        break
    else
        # If the command crashed, wait for a while before running it again
        echo "Snakemake command crashed. Restarting in 5 seconds..."
        sleep 5
    fi
done
