#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:27:56 2019

@author: Caroline
"""

import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import csv

### Once the data is log-transformed, this is where you will read in the gene abundance data file (saved as a csv)
df = pd.read_csv('transformed_data.csv',sep='\t')

#This transposes the data, switching the columns and the rows
df = df.T
#print(df[0:10])

### make a list of all the genes in the study, make sure to rename 'newList to something more sensible
geneList = list(df.iloc[0])

### renames all columns so that they begin with the name of each gene
df.columns = df.iloc[0]
df = df[1:]
#the 1: is reading everyhing after the 0 (or 1?) column


### make two list variables to hold the spearman correlation and the p values of each calculation
list_correlations = []
list_pvalues = []



### loop through the data and calculate spearman correlation and p value for each gene against ENOG4111PNQ (baiE)
for i in range(0, len(geneList)):
    list_correlations.append(scipy.stats.spearmanr(df['COG1902'], df[geneList[i]])[0])
    list_pvalues.append(scipy.stats.spearmanr(df['COG1902'], df[geneList[i]])[1])
### This is where you'll change the name of the gene for each of the 10-11 genes we're interested in.
### 'ENOG411PNQ' will be replaced by each other gene in the list of genes. That should be all you need to change in this loop.

### Write results to a csv file. Column 1 is the correlation, Column 2 is the P value
### Should change the name of the file that saves all of the p values and correlation coefficients for each gene.
with open('COG1902_correlations_and_p_values.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(list_correlations, list_pvalues))
    quit()