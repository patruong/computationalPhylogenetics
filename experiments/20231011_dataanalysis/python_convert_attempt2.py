#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 14 23:39:06 2023

@author: ptruong
"""

import json
import re

# Open the file for reading
with open('settings_sample.file', 'r') as file:
    # Read the entire file contents into a string
    file_contents = file.read()

# Print the file contents
print(file_contents)

def parse_text(text):
    # Remove unnecessary whitespace
    text = re.sub(r'\s+', ' ', text)

    # Initialize variables
    result = {}
    current_dict = result
    stack = []

    # Regular expression pattern to match key-value pairs
    pattern = r'(\w+):\s*{([^}]*)}|"([^"]+)":\s*"(.*?)"'

    # Find all key-value pairs and iterate through them
    for match in re.finditer(pattern, text):
        if match.group(1):
            # This is a dictionary element
            key, value = match.group(1), match.group(2)
            current_dict[key] = parse_list(value)
        else:
            # This is a key-value pair
            key, value = match.group(3), match.group(4)
            current_dict[key] = value

    return result

def parse_list(text):
    # Remove curly braces and split by whitespace
    elements = re.sub(r'[{}]', '', text).split()
    return [element.strip() for element in elements]

input_text = file_contents
parsed_data = parse_text(input_text)
print(parsed_data)
