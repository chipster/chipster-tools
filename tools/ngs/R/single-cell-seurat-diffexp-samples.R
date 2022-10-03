# TOOL single-cell-seurat-diffexp-samples.R: "Seurat v4 -Find conserved cluster markers and DE genes in multiple samples" (This tool lists the cell type markers that are conserved across the conditions, and the differentially expressed genes between the conditions for a user defined cluster. In case of more than 2 samples, each sample is compared to all the other samples. This tool can be used for Seurat objects with 2 or more samples in them.) 
# INPUT OPTIONAL combined_seurat_obj.Robj: "Combined Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL de-list.tsv
# OUTPUT OPTIONAL de-list_{...}.tsv
# OUTPUT OPTIONAL conserved_markers.tsv
# PARAMETER OPTIONAL normalisation.method: "Normalisation method used previously" TYPE [LogNormalize:"Global scaling normalization", SCT:"SCTransform"] DEFAULT LogNormalize (Which normalisation method was used in preprocessing, Global scaling normalization \(default, NormalizeData function used\) or SCTransform.)
# PARAMETER cluster: "Name of the cluster" TYPE STRING DEFAULT 3 (Name of the cluster of which you want to identify the differentially expressed of. By default, the clusters are named with numbers starting from 0.)
# PARAMETER OPTIONAL only.positive: "Return only positive marker genes" TYPE [FALSE, TRUE] DEFAULT TRUE (Tool only returns positive markers as default. Change the parameter here if you want to also include the negative markers.)
# PARAMETER OPTIONAL logFC.conserved: "Fold change threshold for conserved markers in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis. )
# PARAMETER OPTIONAL logFC.de: "Fold change threshold for differentially expressed genes in log2 scale" TYPE DECIMAL FROM 0 TO 5 DEFAULT 0.25 (Genes with an average fold change smaller than this are not included in the analysis.)
# PARAMETER OPTIONAL pval.cutoff.conserved: "p-value cutoff for conserved markers" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the p-value of the conserved cluster marker genes: by default, p-values bigger than 0.05 in any sample are filtered out. This is looking at the max_pval column of the results table.)
# PARAMETER OPTIONAL pval.cutoff.de: "Adjusted p-value cutoff for differentially expressed genes" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (Cutoff for the adjusted p-value of the DE genes: by default, adjusted p-values bigger than 0.05 are filtered out.)
# PARAMETER OPTIONAL mincellsconserved: "Minimum number of cells in one of the groups for conserved markers" TYPE INTEGER DEFAULT 3 (How many cells at least there needs to be in a each sample in the cluster in question.)
# RUNTIME R-4.2.0-single-cell



# 2018-16-05 ML 
# 11.07.2019 ML Seurat v3
# 23.09.2019 EK Add only.pos = TRUE
# 30.10.2019 ML Add filtering parameters
# 02.12.2019 EK Change FC filter to use logfc.threshold prefiltering, add adjusted p-value filtering for DE genes
# 2021-10-04 ML Update to Seurat v4
# 2022-07-21 ML Tune for SCTransform data
# 2022-09-22 ML Fix conserved markers filtering for multiple sample case, add sanity check for cluster number
# 2022-09-20 ML Add min.cells.group parameter  


library(Seurat)

# Load the R-Seurat-objects
load("combined_seurat_obj.Robj")

if (exists("seurat_obj")) {
  data.combined <- seurat_obj
}

# When SCTransform was used to normalise the data, do a prep step:
if (normalisation.method == "SCT"){
  data.combined <- PrepSCTFindMarkers(data.combined)
}

# Check that this cluster is available in data:
if is.na(match(cluster, levels(Idents(data.combined))) ){
  stop("CHIPSTER-NOTE: Cluster given as input can not be found in the Seurat object!")
}

# Identify conserved cell type markers
# (uses package "metap" instead of metaDE since Seurat version 2.3.0)
DefaultAssay(data.combined) <- "RNA" # this is very crucial.
cluster.markers <- FindConservedMarkers(data.combined, ident.1 = cluster, grouping.var = "stim", only.pos = only.positive,  
    verbose = FALSE, logfc.threshold = logFC.conserved, min.cells.group = mincellsconserved)


# Filter conserved marker genes based on adj p-val:
# dat2 <- subset(cluster.markers, (CTRL_p_val_adj < pval.cutoff.conserved & STIM_p_val_adj < pval.cutoff.conserved))
# dat2 <- subset(cluster.markers, (cluster.markers[,5] < pval.cutoff.conserved & cluster.markers[,10] < pval.cutoff.conserved))
dat2 <- subset(cluster.markers, cluster.markers[,"max_pval"] < pval.cutoff.conserved)

# Write to table
write.table(dat2, file = "conserved_markers.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


# Differentially expressed genes across conditions for the cluster (defined by the user, for example cluster 3 -> "3")
data.combined$celltype.stim <- paste(Idents(data.combined), data.combined$stim, sep = "_")
data.combined$celltype <- Idents(data.combined)
Idents(data.combined) <- "celltype.stim"

lvls <- levels(as.factor(data.combined$stim))

if (length(lvls) < 2) { 
  	stop("CHIPSTER-NOTE: There are fewer than 2 samples in the data.")

# If there are only two samples:
} else if (length(lvls) == 2) { 
  ident1 <- paste(cluster, "_", lvls[1], sep = "")
  ident2 <- paste(cluster, "_", lvls[2], sep = "")
  if (normalisation.method == "SCT"){
    cluster_response <- FindMarkers(data.combined, assay = "SCT", ident.1 = ident1, ident.2 = ident2, verbose = FALSE, logfc.threshold = logFC.de)
  } else { 
    cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, logfc.threshold = logFC.de)
  }
  # Filter DE genes based on adj p-val:
  de2 <- subset(cluster_response, (p_val_adj < pval.cutoff.de))

  # Write to table
  write.table(de2, file = "de-list.tsv", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

# If there are more than 2 samples in the data:
} else { 

  for (i in 1:length(lvls) ) { 
    # i:th sample vs all the others (= -i)
    ident1 <- paste(cluster, "_", lvls[i], sep = "")
    ident2 <- paste(cluster, "_", lvls[-i], sep = "")
    cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, logfc.threshold = logFC.de)
    if (normalisation.method == "SCT"){
      cluster_response <- FindMarkers(data.combined, assay = "SCT", ident.1 = ident1, ident.2 = ident2, verbose = FALSE, logfc.threshold = logFC.de)
    } else { 
    cluster_response <- FindMarkers(data.combined, ident.1 = ident1, ident.2 = ident2, verbose = FALSE, logfc.threshold = logFC.de)
    }
    # Filter DE genes based on adj p-val:
    de2 <- subset(cluster_response, (p_val_adj < pval.cutoff.de) )

    # Comparison name for the output file:
    comparison.name <- paste(lvls[i], "vsAllOthers", sep="")
    name.for.output.file <- paste("de-list_", comparison.name, ".tsv", sep="")

    # Write to table
    write.table(de2, file = name.for.output.file, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  } 
} 
# Save the Robj for the next tool
# save(combined_seurat_obj, file="seurat_obj_combined.Robj")

## EOF



