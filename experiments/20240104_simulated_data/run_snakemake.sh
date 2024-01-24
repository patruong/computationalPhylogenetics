#!/bin/bash

while true; do
    snakemake -c10 --rerun-incomplete
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo "Command executed successfully."
        break
    else
        echo "Command failed. Retrying..."
        sleep 5  # Adjust the sleep duration as needed
    fi
done
