# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


# find subtree from tree.

group = "TEST"
tree_str = "((((N172{TEST},N173{REFERENCE})N301{TEST},(N174{REFERENCE},N175{TEST})N302{TEST})N365{TEST},((N176{TEST},N177{TEST})N303{TEST},(N178{TEST},N179{REFERENCE})N304{TEST})N366{REFERENCE})N397{TEST},(((N180{REFERENCE},N181{REFERENCE})N305{TEST},(N182{TEST},N183{REFERENCE})N306{TEST})N367{TEST},((N184{TEST},N185{REFERENCE})N307{REFERENCE},(N186{TEST},N187{TEST})N308{TEST})N368{TEST})N398{TEST})"
tree_str.split(",")

# if everything inside () is TEST return everything inside parenthesis.