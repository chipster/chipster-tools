# TOOL seurat_zifa.R: "Zero Inflated Factor Analysis" (Dimensionality reduction tool for single cell data.)
# INPUT OPTIONAL seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL ZIFAplots.pdf
# OUTPUT OPTIONAL ZIFAgenes.txt
# OUTPUT OPTIONAL seurat_obj_2.Robj
# OUTPUT OPTIONAL coordinates.tsv
# PARAMETER OPTIONAL latentDimensions: "How many latent dimensions is wanted" TYPE INTEGER FROM 1 TO 100 DEFAULT 3 (Number of latent dimensions)
# PARAMETER OPTIONAL p0Thresh: "Filters out genes that are zero in more than this proportion of samples" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.95 (Filters out genes that are zero in more than this proportion of samples)
# RUNTIME R-3.3.2

# PARAMETER OPTIONAL blockNumber: "To how many blocks genes are divided into." TYPE INTEGER FROM 0 TO 100000 DEFAULT 0 (How many blocks is used when running in block mode. If zero, the default block size is used that is: number of genes divided by 500.)
# PARAMETER OPTIONAL block: "To fit or not to fit with block algorithm" TYPE [yes:"yes", no:"no"] DEFAULT no (Is block algorithm used to fit ZIFA or not.)
# PARAMETER OPTIONAL singleSigma: "SingleSigma, if True, fit only a single variance parameter." TYPE INTEGER FROM 0 TO 1 DEFAULT 0 (If True, fit only a single variance parameter zero-inflated PPCA rather than a different on for every gene.)


# This script consists from 4 components:
# 1. Import data from Seurat objects to ZIFA as a file
# 2. Run ZIFA in Python
# 3. Move ZIFA's results back to Seurat object
# 4. Visualize the results

# In addition to R packages this script requires also Python because ZIFA is written in Python.
# The following python modules are required:
# ZIFA, from https://github.com/epierson9/ZIFA
# pandas
# numpy
#
# ZIFA also requires: 
# pylab
# scipy
# numpy
# scikits.learn

library(Seurat) # For Seurat objects
library(rPython) # For Python integration
library(plyr) # For generic name conversion 


# 1.
# Move data from Seurat object into a file so ZIFA can use it.

# Load the object that contains the data
load("seurat_obj.Robj")
# Get a list of wanted genes
genes <- seurat_obj@var.genes
# Filter genes in the same way as in Seurat's PCA script
# Take only unique ones
genes <- unique(genes[genes%in%rownames(seurat_obj@scale.data)])
# Calculate variance
genes.var <- apply(seurat_obj@scale.data[genes,],1,var)
# Use only genes that have larger than 0 variance, i.e. there are different values
genes.use <- genes[genes.var>0]
# Take only non NAs
genes.use <- genes.use[!is.na(genes.use)]
# Read data from raw.data, take only genes specified by var.genes, and the same cells as in scale.data
# ZIFA wants raw count data or log2 data, lets use raw count data so we avoid conversion errors
zifa.data <- seurat_obj@raw.data[genes.use, colnames(seurat_obj@scale.data)]
# Create matrix
zifa.matrix  <- matrix(data=zifa.data, nrow=zifa.data@Dim[1], byrow=FALSE, dimnames=zifa.data@Dimnames)
# Save the .tsv for ZIFA
# This file's columns are cells and rows are genes
# File contains the names of the cells and genes
write.table(zifa.matrix, file='seurat_to_zifa.tsv', quote=FALSE, sep='\t', row.names=TRUE)


# 2.
# Here we use rPython to run ZIFA. ZIFA is written in Python therefore we have to use python.
# ZIFA's results are saved into files: rotated_data.tsv and rotation_matrix.tsv
# rotated_data.tsv contains the new coordinates of the cells in the latent dimensions
# rotation_matrix is the matrix that is used to map high-dimensional gene space data to the latent space. 
# The map from gene space to latent space is just a matrix multiplication: rotation_matrix * gene_data = rotated_data
# User defines the values of latentDimensions and p0Thresh

python.assign('latent_dimensions', latentDimensions)
#python.assign('n_blocks', 0)
python.assign('p0_thresh', p0Thresh)
python.assign('single_sigma', 0)
# Import required modules
python.exec('from ZIFA import ZIFA, block_ZIFA')
python.exec('import pandas as pd')
python.exec('import numpy as np')
# Read the data from a file
# We have to escape the apostrophes in the command, because it is wrapped inside of quotation itself
python.exec('zifa_data = pd.read_csv(\'seurat_to_zifa.tsv\', sep=\'\\t\', header=0)')
# Convert expression data to correct format, count -> log2
python.exec('zifa_data = np.log2(zifa_data + 1)')
# Filter the non expressive genes away the same way as it is done in ZIFA.
# By doing this we do not lose the names of the samples when running ZIFA.
python.exec('zifa_data = zifa_data.loc[np.asarray((np.abs(zifa_data.values) < 1e-6).mean(axis = 1) <= p0_thresh),:]')
# Read the cells's and genes' names
python.exec('cells = zifa_data.columns')
python.exec('genes = zifa_data.index')
# Transpose the data, because ZIFA wants the data in different form
python.exec('zifa_data = zifa_data.T')
# Run ZIFA
python.exec('rotated_data, params = block_ZIFA.fitModel(zifa_data.values, latent_dimensions, p0_thresh=p0_thresh, singleSigma=single_sigma)')
# Write results into a file
# Create a dataframe, so  we can write cells' names
python.exec('rotated_data = pd.DataFrame(data = rotated_data, index = cells)')
# This monster opens a file 'zifa_results.tsv' and writes the df dataframe into it with a to_csv function using tabs and writing also the cells' names (indexes)
# \n is linebreak and \t is tab
python.exec('with open(\'rotated_data.tsv\', \'w\') as f: \n\t rotated_data.to_csv(f, sep=\'\\t\', index=True, mode=\'w\')')
# We want also the rotation matrix, that is used to produce the low dimensional linear combination from the original data
python.exec('rotation_matrix = pd.DataFrame(data = params[\'A\'], index = genes)')
# Write the rotation matrix in the same way as the rotated data aka zifa_results.tsv above
python.exec('with open(\'rotation_matrix.tsv\', \'w\') as f: \n\t rotation_matrix.to_csv(f, sep=\'\\t\', index=True, mode=\'w\')')


# 3.
# Read ZIFA's outputs in R and insert the data into Seurat object

zifa.rotated.data <- read.table(file='rotated_data.tsv', sep='\t', header=TRUE, row.names=1)
zifa.rotation.matrix <- read.table(file='rotation_matrix.tsv', sep='\t', header=TRUE, row.names=1)
# Rename coordinates, so they are in the same format that PCA produces in Seurat e.g. "PC1", "PC2", ...
latent.dimensions <- ncol(zifa.rotated.data)
for (dimension in 1:latent.dimensions) {
	new.name <- paste('PC', dimension, sep="")
	colnames(zifa.rotated.data)[dimension] <- new.name
	colnames(zifa.rotation.matrix)[dimension] <- new.name
}
# Calculate the standard deviation of the cells in each dimension
# Standard deviation is not so important for ZIFA as it is for PCA, but 
# lets calculate it, so we can visualize it if we want to
# We have to also rename it, so it is compatible with Seurat
zifa.sdev <- apply(zifa.rotated.data, 2, sd)
zifa.sdev <- data.frame(zifa.sdev)
zifa.sdev <- rename(zifa.sdev, c('zifa.sdev'='sdev'))
# Assign rotated data into Seurat objects pca.rot
# NOTE Seurat naming conventions are confusing, .rot is usually the rotation matrix, not the rotated data
seurat_obj@pca.rot <- zifa.rotated.data
seurat_obj@pca.x <- zifa.rotation.matrix
seurat_obj@pca.obj <- list(zifa.sdev)


# 4. 
# Visualize and print results the same way as in the PCA tool

# Print results into a ZIFAgenes.txt
sink('ZIFAgenes.txt')
PrintPCA(seurat_obj, pcs.print=1:2, genes.print=5, use.full=FALSE)
sink()
# Save visualizations into a pdf
pdf(file="ZIFAplots.pdf")
VizPCA(seurat_obj, pcs.use=1:2, num.genes=30, use.full=FALSE)
PCAPlot(seurat_obj, 1, 2)
PCHeatmap(seurat_obj, pc.use=1, cells.use=20, do.balanced=TRUE)
PCElbowPlot(seurat_obj)
# Close the pdf file
dev.off()
# Save the Roj for the next tool
save(seurat_obj, file="seurat_obj_2.Robj")

#EOF