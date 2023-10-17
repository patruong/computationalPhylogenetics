#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 23:39:06 2023

@author: ptruong
"""

import json

# Open the file for reading
with open('settings_sample.file', 'r') as file:
    # Read the entire file contents into a string
    file_contents = file.read()

# Print the file contents
print(file_contents)

data = file_contents

# Clean up and convert the data to a Python dictionary
data = data.replace("{\n", "[\n")
data = data.replace("}\n", "],\n")
data = data.replace(":", ": ")
data = data.replace("{{", "[[")
data = data.replace("}}", "]],")
data = data.replace("{", "[")
data = data.replace("}", "]")
#data = data.replace("
data = data.rstrip(',')  
                    
# Convert the cleaned data to a Python list
data_list = eval(data)

# Convert the list to a JSON string
json_data = json.dumps(data, indent=4)

# You now have the data in valid JSON format as a string in the 'json_data' variable
print(json_data)

# Save the JSON data to a file
with open("output.json", "w") as json_file:
    json_file.write(json_data)

print("JSON data saved to 'output.json'")

dict(json_data)