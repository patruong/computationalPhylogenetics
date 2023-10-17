#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 15 00:20:57 2023

@author: ptruong
"""

import re
import json
import pandas as pd
import argparse


def replace_single_value_curly_brackets(s):
    r = r"\{(\d*\.?\d*)\}" # the regular expression to match a number inside curly brackets
    t = "[\\1]" # the replacement string to replace the matched number with square brackets
    result = re.sub(r, t, s) # the result of the replacement
    return result # print the result
     
def replace_vector_value_curly_brackets(s):
    #s = "{0.77, 0, 0, 0, 0}" # the original string  
    r = r"\{([\d\., ]*)\}" # the regular expression to match a comma-separated list of numbers inside curly brackets
    t = "[\\1]" # the replacement string to replace the matched list with square brackets
    result = re.sub(r, t, s) # the result of the replacement
    return result # print the result

def fix_Q_key_curly_brackets(s):
    # Define a regular expression to match and replace inner curly brackets
    pattern = r'(\{.*?\})'
    # Replace inner curly brackets with square brackets
    output_text = re.sub(pattern, lambda m: m.group(0).replace('{', '[').replace('}', ']'), s)
    return output_text

def add_comma_between_square_brackets(s):
    return re.sub(r'(?<=\])\s+(?=\[)', ',', s)

def fix_curly_brackets_in_some_keys(s):
    # Define a list of keys to transform
    keys_to_transform = ["EFV", "Q", "alphabet", "bases","stop", 
                         "rates", "meta", "model_list"]
    
    # Use a loop to transform the keys in the sample
    transformed_sample = s
    
    for key in keys_to_transform:
        # Use regex to perform the transformation for the specific key
        pattern = rf'"{key}"\s*:\s*{{([^}}]*)}}'
        replacement = rf'"{key}": [\1]'
        transformed_sample = re.sub(pattern, replacement, transformed_sample)
    
    
    transformed_sample = re.sub(r',(?=\s*})', '', transformed_sample)
    return transformed_sample


def convert_hyphy_simulation_settings_to_JSON(input_file):
    # Open the file for reading
    with open(input_file, 'r') as file:
        # Read the entire file contents into a string
        file_contents = file.read()
        
    res = replace_single_value_curly_brackets(file_contents)
    res = replace_vector_value_curly_brackets(res)
    res = fix_Q_key_curly_brackets(res) 
    res = add_comma_between_square_brackets(res)
    res = fix_curly_brackets_in_some_keys(res)
    
    
    json_file = json.loads(res)
    json_dict= dict(json_file)
    return json_dict

    # Define the output file path
    #output_file_path = "output.json"
    
    # Write the data to the output file
    #with open(output_file_path, "w") as output_file:
    #    json.dump(json_file, output_file, indent=2)
        
def increment_beta_index_by_1(input_dict):
    new_dict = {}
    
    for old_key, value in input_dict.items():
        new_key = str(int(old_key) + 1)  # Increment the old key by 1
        new_dict[new_key] = value
    
    return new_dict

def get_simulation_settings_from_hyphy(input_file):
    data = convert_hyphy_simulation_settings_to_JSON(input_file)       
    
    header = {'0': "alpha"}
    header_betas = increment_beta_index_by_1(data["simulator.site.profile"]["variables"])
    header.update(header_betas)
    
    sim_settings = pd.DataFrame(data["simulator.site.profile"]["rates"], columns=header.values()).reset_index().rename({"index":"site"}, axis=1)
    sim_settings["site"] += 1
    return sim_settings


def main(input_file, output = "hyphy_sim_settings.tsv"):
    df = get_simulation_settings_from_hyphy(input_file)
    df.to_csv(output, sep = "\t", index=False)
    print(f'HyPhy simulator.site.profile written to {output}.')

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process simulation settings from HyPhy")

    # Positional argument for the input file
    parser.add_argument("input_file", help="Input file containing simulation settings")
    
    # Optional argument for specifying the output file (default is "hyphy_sim_settings.tsv")
    parser.add_argument("--output", default="hyphy_sim_settings.tsv", help="Output file name")
    
    args = parser.parse_args()
    
    # Call the main function with the provided arguments
    main(args.input_file, args.output)

