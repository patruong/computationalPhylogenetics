#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 12:09:25 2023

@author: ptruong
"""

import os
import pandas as pd

experiment_dir = "output/omnibus/pos_thresh_0.95"
df = pd.DataFrame()
for file_dir in os.listdir(experiment_dir):
    if "time_results.csv" in os.listdir(experiment_dir + "/" + file_dir):
        partial_df = pd.read_csv(experiment_dir + "/" + file_dir + "/time_results.csv")
        df = pd.concat([df, partial_df])
        
df.to_csv("times.csv", sep = "\t")
