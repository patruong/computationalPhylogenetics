#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 22:54:54 2023

@author: ptruong
"""

smaller_sample_incorrect = """
{
    "model":{
      "EFV":    {
   [0.01929014233974128],[0.01622765628659735],[0.01769443322434847],[0.01623235001696338],[0.01622972861556908],[0.01365311115696728],[0.01488718145152558],[0.0136570602190682],[0.01769333756976861],[0.01488435883925134],[0.01622971850752687],[0.01488866403069913],[0.01623137334257151],[0.01365449476855555],[0.01488869012427688],[0.01365844423085596],[0.01929347190999687]
       },
       "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
       "ID":"simulator.substitution_model",
       "Q":    {
   ["", "simulator.substitution_model.theta_AC*beta*0.23367778810339", "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", "simulator.substitution_model.theta_AT*beta*0.2337453776869083", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],["simulator.substitution_model.theta_AC*beta*0.2777775", "", "simulator.substitution_model.theta_CG*beta*0.2547993342097017", "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2990038356551464", "", "", "", "", "", "", "", "", "", "", "", ""]
       }}}
"""

      
smaller_sample_corrected= """
{
    "model":{
        "EFV":    { "values" : [
            [0.01929014233974128],[0.01622765628659735],[0.01769443322434847],
            [0.01623235001696338],[0.01622972861556908],[0.01365311115696728],
            [0.01488718145152558],[0.0136570602190682],[0.01769333756976861],
            [0.01488435883925134],[0.01622971850752687],[0.01488866403069913],
            [0.01623137334257151],[0.01365449476855555],[0.01488869012427688],
            [0.01365844423085596],[0.01929347190999687]
            ]
        },
       "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
       "ID":"simulator.substitution_model",
       "Q":    { "values" : [
   ["", "simulator.substitution_model.theta_AC*beta*0.23367778810339", 
    "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", 
    "simulator.substitution_model.theta_AT*beta*0.2337453776869083", 
    "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", 
    "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", 
    "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", 
    "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", 
    "", "", "", "", "", "", "", "", "", "", "", "", 
    "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", 
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],
    ["simulator.substitution_model.theta_AC*beta*0.2777775", "", 
     "simulator.substitution_model.theta_CG*beta*0.2547993342097017", 
     "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", 
     "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", 
     "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", 
     "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", 
     "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", 
     "", "", "", "", "", "", "", "", "", "", "", "", 
     "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", 
     "", "", "", "", "", "", "", "", "", "", "", 
     "simulator.substitution_model.theta_AT*beta*0.2990038356551464", "", "", "", 
     "", "", "", "", "", "", "", "", ""]
    ]
       }
    }
}
"""


smaller_sample_corrected_version_2= """
{
    "model":{
        "EFV":    [
            [0.01929014233974128],[0.01622765628659735],[0.01769443322434847],
            [0.01623235001696338],[0.01622972861556908],[0.01365311115696728],
            [0.01488718145152558],[0.0136570602190682],[0.01769333756976861],
            [0.01488435883925134],[0.01622971850752687],[0.01488866403069913],
            [0.01623137334257151],[0.01365449476855555],[0.01488869012427688],
            [0.01365844423085596],[0.01929347190999687]
            ],
       "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
       "ID":"simulator.substitution_model",
       "Q":    [
   ["", "simulator.substitution_model.theta_AC*beta*0.23367778810339", 
    "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", 
    "simulator.substitution_model.theta_AT*beta*0.2337453776869083", 
    "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", 
    "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", 
    "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", 
    "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", 
    "", "", "", "", "", "", "", "", "", "", "", "", 
    "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", 
    "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],
    ["simulator.substitution_model.theta_AC*beta*0.2777775", "", 
     "simulator.substitution_model.theta_CG*beta*0.2547993342097017", 
     "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", 
     "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", 
     "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", 
     "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", 
     "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", 
     "", "", "", "", "", "", "", "", "", "", "", "", 
     "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", 
     "", "", "", "", "", "", "", "", "", "", "", 
     "simulator.substitution_model.theta_AT*beta*0.2990038356551464", "", "", "", 
     "", "", "", "", "", "", "", "", ""]
    ]
    }
}
"""

# Print the corrected JSON
print(smaller_sample)
json_file = json.loads(smaller_sample_corrected_version_2)


type(json_file)




#3####

import re

# Sample input
smaller_sample_1 = """
{
    "model":{
      "EFV":    {
   [0.01929014233974128],[0.01622765628659735],[0.01769443322434847],[0.01623235001696338],[0.01622972861556908],[0.01365311115696728],[0.01488718145152558],[0.0136570602190682],[0.01769333756976861],[0.01488435883925134],[0.01622971850752687],[0.01488866403069913],[0.01623137334257151],[0.01365449476855555],[0.01488869012427688],[0.01365844423085596],[0.01929347190999687]
       },
       "Equilibrium frequency estimator":"Corrected 3x4 frequency estimator",
       "ID":"simulator.substitution_model",
       "Q":    {
   ["", "simulator.substitution_model.theta_AC*beta*0.23367778810339", "simulator.substitution_model.theta_AG*alpha*0.2547993342097017", "simulator.substitution_model.theta_AT*beta*0.2337453776869083", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""],["simulator.substitution_model.theta_AC*beta*0.2777775", "", "simulator.substitution_model.theta_CG*beta*0.2547993342097017", "simulator.substitution_model.theta_CT*alpha*0.2337453776869083", "", "simulator.substitution_model.theta_AC*beta*0.2337076295815298", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2547835568149736", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2337313136034966", "", "", "", "simulator.substitution_model.theta_AC*beta*0.2336921683896183", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AG*beta*0.2336521569835915", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "simulator.substitution_model.theta_AT*beta*0.2990038356551464", "", "", "", "", "", "", "", "", "", "", "", ""]
       }}}
"""

# Define a list of keys to transform
keys_to_transform = ["EFV", "Q"]

# Use a loop to transform the keys in the sample
transformed_sample = smaller_sample_1

for key in keys_to_transform:
    # Use regex to perform the transformation for the specific key
    pattern = rf'"{key}"\s*:\s*{{([^}}]*)}}'
    replacement = rf'"{key}": [\1]'
    transformed_sample = re.sub(pattern, replacement, transformed_sample)

# Print the result
print(transformed_sample)

json.loads(transformed_sample)


import re

text = """
"model_map":{
  "N0":"GROUP_1",
  "N1":"GROUP_2",
  "N10":"GROUP_3",
  "N100":"GROUP_3",
  "N101":"GROUP_3",
  "N102":"GROUP_2",
  "N103":"GROUP_0",
  "N104":"GROUP_0",
  "N105":"GROUP_0",
  "N106":"GROUP_3",
  "N107":"GROUP_3",
  "N108":"GROUP_3",
  "N109":"GROUP_2",
  "N11":"GROUP_3",
  "N110":"GROUP_0",
  "N111":"GROUP_2",
  "N112":"GROUP_0",
  "N113":"GROUP_0",
  "N114":"GROUP_3",
  "N115":"GROUP_2",
  "N116":"GROUP_2",
  "N117":"GROUP_2",
  "N118":"GROUP_1",
  "N119":"GROUP_2"
 }
"""

# Use regex to remove commas before }
result = re.sub(r',(?=\s*})', '', text)

# Print the result
print(result)
