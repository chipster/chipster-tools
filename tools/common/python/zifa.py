# Import required modules
from ZIFA import block_ZIFA
import pandas as pd
import numpy as np
import sys

# Read arguments
# Arguments:  
# 1: latent_dimensions, 
# 2: p0_thresh, 
# 3: single_sigma
# 4: path to input tsv
# 5: path to output, rotated_data
# 6: path to output, rotation_matrix
latent_dimensions = int(sys.argv[1])
#python.assign('n_blocks', 0)
p0_thresh = float(sys.argv[2])
single_sigma = int(sys.argv[3])
input_tsv_path = sys.argv[4]
rotated_data_path = sys.argv[5]
rotation_matrix_path = sys.argv[6]
# Read the data from a file
zifa_data = pd.read_csv(input_tsv_path, sep='\t', header=0)
# Convert expression data to correct format, count -> log2
zifa_data = np.log2(zifa_data + 1)
# Filter the non expressive genes away the same way as it is done in ZIFA.
# By doing this we do not lose the names of the samples when running ZIFA.
zifa_data = zifa_data.loc[np.asarray((np.abs(zifa_data.values) < 1e-6).mean(axis = 1) <= p0_thresh),:]
# Read the cells's and genes' names
cells = zifa_data.columns
genes = zifa_data.index
# Transpose the data, because ZIFA wants the data in different form
zifa_data = zifa_data.T
# Run ZIFA
rotated_data, params = block_ZIFA.fitModel(zifa_data.values, latent_dimensions, p0_thresh=p0_thresh, singleSigma=single_sigma)
# Write results into a file
# Create a dataframe, so  we can write cells' names
rotated_data = pd.DataFrame(data = rotated_data, index = cells)
# Open 'zifa_results.tsv' file and write the df dataframe into it with a to_csv function using tabs and writing also the cells' names (indexes)
with open(rotated_data_path, 'w') as f: 
    rotated_data.to_csv(f, sep='\t', index=True, mode='w')
# We want also the rotation matrix, that is used to produce the low dimensional linear combination from the original data
rotation_matrix = pd.DataFrame(data = params['A'], index = genes)
# Write the rotation matrix in the same way as the rotated data aka zifa_results.tsv above
with open(rotation_matrix_path, 'w') as f:
    rotation_matrix.to_csv(f, sep='\t', index=True, mode='w')
