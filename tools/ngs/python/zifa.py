# TOOL zifa.py: "Zero Inflated Factor Analysis" (Dimensionality reduction tool for single cell data.)
# INPUT data: "Observed zero-inflated data" TYPE GENERIC
# OUTPUT OPTIONAL zifa.log
# OUTPUT OPTIONAL output.tsv
# OUTPUT OPTIONAL plot.pdf
# PARAMETER OPTIONAL latentDimension: "How many latent dimensions is wanted" TYPE INTEGER FROM 1 TO 100 DEFAULT 3 (Number of latent dimensions)
# PARAMETER OPTIONAL block: "To fit or not to fit with block algorithm" TYPE [yes:"yes", no:"no"] DEFAULT no (Is block algorithm used to fit ZIFA or not.)
# PARAMETER OPTIONAL blockNumber: "To how many blocks genes are divided into." TYPE INTEGER FROM 0 TO 100000 DEFAULT 0 (How many blocks is used when running in block mode. If zero, the default block size is used that is: number of genes divided by 500.)
# PARAMETER OPTIONAL dataType: "Is the input data just count data or is it the logarithm of the count data." TYPE [count:"count", logarithmic:"logarithmic"] DEFAULT count (Is the input data just count data or is it the logarithm of the count data.)
# PARAMETER OPTIONAL header: "Is the first row and column of input data text, not actual data." TYPE [header:"header", noHeader:"no header"] DEFAULT header (Is there a header in the input data; first row and columns are not actual data)
# PARAMETER OPTIONAL p0Thresh: "Filters out genes that are zero in more than this proportion of samples" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.95 (Filters out genes that are zero in more than this proportion of samples)
# PARAMETER OPTIONAL singleSigma: "SingleSigma, if True, fit only a single variance parameter." TYPE INTEGER FROM 0 TO 1 DEFAULT 0 (If True, fit only a single variance parameter zero-inflated PPCA rather than a different on for every gene.)

# AO 4.7.2017 First version

import sys
import matplotlib
# Force matplotlib to not use any XWindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from numpy import genfromtxt, savetxt, log2, transpose
from math import log
import pandas as pd

from ZIFA import ZIFA
from ZIFA import block_ZIFA

# Copy std out to zifa.log
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("zifa.log", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        #this flush method is needed for python 3 compatibility.
        #this handles the flush command by doing nothing.
        #you might want to specify some extra behavior here.
        pass    

# Initialize Logger
sys.stdout = Logger() 

## Check parameters
if header == "header":
    # Read the first row and column, without the first object that is the corner
    # read a file 'data' that is in .tsv format
#    myData = genfromtxt('data', delimiter = '\t', skip_header = 1)
    myData = pd.read_csv('data', sep = '\t', header = 1)
#    cells = myData.columns.values
    # delete first column, first row is already trimmed away
#    myData = myData[1:,:]
else:
    # TODO: use also pandas
    # No header
    myData = genfromtxt('data', delimiter = '\t', skip_header = 0)
if dataType == "count": 
    myData = log2(myData + 1)
if block == "no":
    # run ZIFA
    Zhat, params = ZIFA.fitModel(transpose(myData), latentDimension, singleSigma) 
else:
    if blockNumber == 0:
        # run block_ZIFA, with default number of blocks
        Zhat, params = block_ZIFA.fitModel(transpose(myData), latentDimension, p0_thresh = p0Thresh, singleSigma = singleSigma)
    else:
        # run block_ZIFA 
        Zhat, params = block_ZIFA.fitModel(transpose(myData), latentDimension, n_blocks = blockNumber, p0_thresh = p0Thresh, singleSigma = singleSigma) 

## Write results to a file
# open a file called workfile, with write permissions 'w'
f = open('output.tsv', 'w') 

# print myData, the input, to std out that is also copied into zifa.log
print(myData) 

# print the output of ZIFA, this shows in zifa.log
print(Zhat) 

## OLD TODO: remove this comment
# save results into tsv file,delimiter is \t and the file name that is associated with file handle f is output.tsv. fmt='%f' means that numbers are in floating point format.

# if there was a header in the input include it also in output
#savetxt(f, Zhat, delimiter='\t', fmt='%f', header = 
# Use pandas module to write "indecies" aka row names into file, lets first write it to a csv, that is the filetype that pandas know how to handle, and then convert it to a tsv that is the default filetype in Chipster.
# Convert myData into data frame, that can be saved as tsv with index

# TODO: add pca display option
if header == "header":
#    df = pd.DataFrame(data = Zhat, index = cells)
    df = pd.DataFrame(data = Zhat)
    df.to_csv(f, sep = '\t')
# Write without index
else:
    df = DataFrame(data = Zhat)
    df.to_csv(f, sep = '\t')

fig = plt.figure(figsize = [5, 5])
plt.scatter(Zhat[:,0], Zhat[:,1], color = 'red', s = 4)    
plt.xlim([-4,4])
plt.ylim([-4,4])
plt.title('ZIFA Estimated Latent Positions')
fig.savefig('plot.pdf')


#EOF
