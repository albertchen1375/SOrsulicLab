import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from pandas import ExcelWriter


df = pd.read_csv('Cell Line Gene Expression Batch  Corrected Transposed.csv')

with open('protein_coding.txt', 'r') as file:
    lines = file.readlines()

# Initialize an empty dictionary
gene_dict = {}

# Process each line
for line in lines:
    # Strip any extra whitespace characters like newline
    gene_symbol = line.strip()
    # Add the gene symbol as a key in the dictionary with a placeholder value
    if gene_symbol:
        gene_dict[gene_symbol] = 1  # Or assign any value you need


df.set_index('depmap_id', inplace=True)

# Split the DataFrame into rows before and after index 9
df_before_row_6 = df.iloc[:7]  

df_after_row_6 = df.iloc[7:]   

# Create a boolean mask where True represents rows to keep
mask = df_after_row_6.index.to_series().apply(lambda x: x in gene_dict.keys())

# Apply the mask to filter the DataFrame
df_non_protein_coding_after_row_9 = df_after_row_6[mask]

df_combined = pd.concat([df_before_row_6, df_non_protein_coding_after_row_9])



# Reset index to turn 'depmap_id' back into a column if needed
df_combined.reset_index(inplace=True)

df_combined.to_csv('Cell Line Gene Expression_non_protein_coding_removed.csv', header= True)

