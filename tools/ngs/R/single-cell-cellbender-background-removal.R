# TOOL single-cell-cellbender-background-removal.R: "Seurat v4 -Remove background contamination with CellBender" (This tool estimates non-empty cells from raw 10x feature-barcode matrices and removes systematic background contamination from the estimated cells. The raw feature-barcode matrix input must be given in hdf5 file format. The CellBender filtered feature-barcode matrix output of this tool can be used as input to the "Setup and QC" tool.)
# INPUT raw_fb_matrix.h5: "Raw 10x feature-barcode matrix in hdf5 format" TYPE GENERIC ()
# OUTPUT cellbender_fb_matrix_report.html
# OUTPUT cellbender_fb_matrix_filtered.h5
# OUTPUT cellbender_fb_matrix.log
# PARAMETER OPTIONAL cell_num: "Expected number of cells" TYPE STRING DEFAULT "auto" (This is the number of droplets that reasonably surely contain cells. If "auto", CellBender automatically estimates this based on the input dataset. Otherwise, the value should be an integer chosen based on the "UMI curve" plot. )
# PARAMETER OPTIONAL droplet_num: "Total number of included droplets" TYPE STRING DEFAULT "auto" (This is a number that goes a few thousand barcodes into the "empty droplet plateau" and includes some droplets that are reasonably surely empty. If "auto", CellBender automatically estimates this based on the input dataset. Otherwise, the value should be an integer chosen based on the "UMI curve" plot. Keep in mind that the algorithm takes longer to run with larger values.)
# PARAMETER OPTIONAL epoch_num: "Number of epochs" TYPE INTEGER DEFAULT 150 (This is the number of complete passes through the entire training dataset when training the neural network in CellBender. Typically, 150 is a good choice. As a rule of thumb, the value should not exceed 300.)
# PARAMETER OPTIONAL lr: "Learning rate" TYPE DECIMAL DEFAULT 0.0001 (This parameter controls the amount by which the neural network's parameters are updated in the opposite direction of the gradient of the loss function when training the neural network. The default is typically a good choice, but it can be reduced if there are large downward dips in the ELBO score. The value should be between 0 and 1.)
# PARAMETER OPTIONAL fpr: "Nominal false positive rate" TYPE DECIMAL DEFAULT 0.01 (This controls the trade-off between removing noise and retaining signal. Larger values correspond to removing more noise at the expense of more signal. The value should be between 0 and 1. )
# RUNTIME R-4.2.3-cellbender
# SLOTS 2

# The output cellbender_fb_matrix_filtered.h5 of this CellBender tool is formatted 
# exactly like a 10x Cell Ranger v3 hdf5 file so that it is compatible with 
# Seurat v4 dataloader Read10X_h5() in the "Seurat v4 -Setup and QC" -tool

library(rhdf5)

# Cellbender needs to find jupyter from python bin to create html report
path <- Sys.getenv('PATH')
path <- paste(path, ':/opt/chipster/tools-bin/python-3.7.17/bin', sep='')
Sys.setenv(PATH = path)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
cellbender <- paste(chipster.tools.path, "python-3.7.17/bin/cellbender", sep = "/")

# Document version numbers
cellbender.version.command <- paste(cellbender, "-v")
version <- system(cellbender.version.command, intern = TRUE)
documentVersion("CellBender", version)

if (cell_num != "auto" && !grepl("^\\d+$", cell_num)) {
  stop("CHIPSTER-NOTE: The expected number of cells should either be auto or an integer.")
}
if (droplet_num != "auto" && !grepl("^\\d+$", droplet_num)) {
  stop("CHIPSTER-NOTE: The total number of droplets should either be auto or an integer.")
}
if (lr < 0 | lr > 1) {
  stop("CHIPSTER-NOTE: The learning rate should be a value between 0 and 1.")
}
if (fpr < 0 | fpr > 1) {
  stop("CHIPSTER-NOTE: The nominal false positive rate should be a value between 0 and 1.")
}
if (droplet_num < cell_num) {
  stop("CHIPSTER-NOTE: The total number of included droplets should be greater than the expected number of cells.")
}

command <- c(
cellbender, "remove-background",
"--input", "raw_fb_matrix.h5",
"--output", "cellbender_fb_matrix.h5",
"--epochs", as.character(epoch_num),
"--learning-rate", as.character(lr),
"--fpr", as.character(fpr),
"--cpu-threads", chipster.threads.max)

command <- paste(command, collapse=" ")

if (cell_num != "auto") {
  command <- paste(command, paste(c("--expected-cells", cell_num), collapse=" "))
}
if (droplet_num != "auto") {
  command <- paste(command, paste(c("--total-droplets-included", droplet_num), collapse=" "))
}

documentCommand(command)

runExternal(command)

# Extra information in CellBender's cellbender_fb_matrix_filtered.h5 must be deleted so 
# that it is formatted exactly like a Cell Ranger v3 h5 file 

h5f <- H5Fopen("cellbender_fb_matrix_filtered.h5")

if (!is.null(h5f$droplet_latents)) {
  h5delete(h5f, "droplet_latents")
}
if (!is.null(h5f$global_latents)) {
  h5delete(h5f, "global_latents")
}
if (!is.null(h5f$metadata)) {
  h5delete(h5f, "metadata")
}

if (is.null(h5f$matrix)) {
  stop("CHIPSTER-NOTE: The CellBender output cellbender_fb_matrix_filtered.h5 does not contain matrix group.
  Output is not compatible with Seurat v4 dataloader Read10X_h5()")
}

if (!file.exists("cellbender_fb_matrix_report.html")) {
  print("This is the .log file output:")
  system("cat cellbender_fb_matrix.log")
  stop("CHIPSTER-NOTE: CellBender should output an html file. However, CellBender does not output the html file which might indicate
  that something went wrong with the CellBender method. Look at the screen output of this job and look for the .log file output that 
  is shown there to see possible warnings and error messages. Also, try rerunning the CellBender method again using the same parameter
  values. If this error occurs again after rerunning with the same parameter values, please ensure that your parameter values are 
  reasonable.")
} 

H5Fclose(h5f) 

# EOF





