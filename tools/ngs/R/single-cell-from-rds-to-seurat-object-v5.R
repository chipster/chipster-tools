# TOOL single-cell-from-rds-to-seurat-object-v5.R: "Seurat v5 -From RDS file to Seurat object" (This tool converts an RDS file into a Seurat object. The RDS file given as input must contain Seurat object information created with Seurat v5 to be converted into a Seurat object.)
# INPUT rds_file: "RDS file" TYPE GENERIC 
# OUTPUT seurat_obj.Robj
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""

options(Seurat.object.assay.version = "v5")
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Read input names
inputnames <- read_input_definitions()

# Check that RDS file is given as input
if (grepl(".rds", inputnames$rds_file, ignore.case=TRUE) == FALSE) {
    stop("CHIPSTER-NOTE: The input file must be an RDS file.")
} 

# Convert rds file into Seurat object
seurat_obj <- readRDS("rds_file") 

# Check that converted RDS file contains a Seurat object
if (class(seurat_obj) != "Seurat") {
  stop("CHIPSTER-NOTE: The RDS file does not contain Seurat object information.")
}

# Check that Seurat object was created with Seurat v4
if (as.integer(substr(seurat_obj@version,1,1)) != 5) {
  stop("CHIPSTER-NOTE: The RDS file does not contain Seurat object information created with Seurat v5. The given RDS file was created with Seurat version ", as.integer(substr(seurat_obj@version,1,1)), ".")
}

# Save the Seurat object for the next tool
save(seurat_obj, file = "seurat_obj.Robj")

# EOF

