# TOOL single-cell-seurat-rename-sample-v5.R: "Seurat v5 BETA -(Re)name sample and samplegroup" (With this tool you can name or rename your samples and samplegroups in your Seurat object..)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL seurat_obj_renamed.Robj
# PARAMETER project.name: "Project name for plotting" TYPE STRING DEFAULT My_project_name (You can give your project a name. The name will appear on the plots. Do not use underscore _ in the names!)
# PARAMETER sample_name: "Sample name" TYPE STRING DEFAULT control1 (Type the sample name or identifier here. For example control1, cancer3a. Do not use underscore _ in the names! Fill this field if you are combining samples later.)
# PARAMETER sample.group: "Sample group" TYPE STRING DEFAULT CTRL (Type the sample name or identifier here. For example CTRL, STIM, TREAT. Do not use underscore _ in the names! Fill this field if you are combining samples later.)
# RUNTIME R-4.3.2-single-cell
# TOOLS_BIN ""

# 2024-11-28 ML


library(Seurat)
options(Seurat.object.assay.version = "v5")

# Load the R-Seurat-object (called seurat_obj)
load("seurat_obj.Robj")

# Replace empty spaces in sample and project name with underscore _ :
project.name <- gsub(" ", "_", project.name)
sample_name <- gsub(" ", "_", sample_name)
sample.group <- gsub(" ", "_", sample.group)

# (Re)name:
seurat_obj@meta.data$stim <- sample_name
seurat_obj$type <- sample.group
seurat_obj@project.name <- project.name

# Save the Robj for the next tool
save(seurat_obj, file = "seurat_obj_renamed.Robj")

## EOF
