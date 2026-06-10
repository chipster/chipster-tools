# TOOL single-cell-seurat-annotate-cells-sctype-v5.R: "Seurat v5 - Annotate cells with ScType" (You can use this tool to annotate clusters using ScType)
# INPUT seurat_obj.Robj: "Seurat object" TYPE GENERIC
# OUTPUT OPTIONAL DimPlot.pdf
# PARAMETER tissuetype: "Tissue type" TYPE [Lung: "Lung", "Immune system": "Immune system"] DEFAULT "Immune system"  (Which CellDex reference to use for annotations.)
# RUNTIME R-4.5.1-seurat5
# TOOLS_BIN ""



#Load needed packages
library("Seurat")
library("dplyr")
library("HGNChelper")
library("openxlsx")
library("ggraph")
library("igraph")
library("tidyverse")
library("data.tree")

# PARAMETER tissuetype: "Tissue type of seurat object" TYPE [Immune System: "Immune system"] DEFAULT Lung (Choose the tissue type of your study)

#tissuetype <- as.character(tissuetype)

load("seurat_obj.Robj")

if (exists("data.combined")) {
  seurat_obj <- data.combined
}

#The following functions are from https://github.com/IanevskiAleksandr/sc-type

gene_sets_prepare <- function(path_to_db_file, cell_type){
  
  cell_markers = openxlsx::read.xlsx(path_to_db_file)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}



#Actual score calculation
sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
  
  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }
  
  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T); 
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)
  
  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }
  
  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){ 
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
  
  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
  
  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }
  
  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
  
  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){ 
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  })) 
  
  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
  
  es.max
}


# Run ScType function:


run_sctype <- function(seurat_object, known_tissue_type = NULL, assay = "RNA", scaled = TRUE, custom_marker_file = NULL, plot = FALSE, name = "sctype_classification") {
  db_=sctype_source()
  # Check for missing arguments
  if (is.null(seurat_object)) {
    stop("Argument 'seurat_object' is missing")
  }
  if (!inherits(seurat_object, "Seurat")) {
    stop("Argument 'seurat_object' must be a Seurat object")
  }
  # Set default custom marker file
  if (is.null(custom_marker_file)) {
    custom_marker_file = db_
  }
  # Auto-detect tissue type if not provided
  if (is.null(known_tissue_type)) {
    print("Guessing tissue type: \n");
    tissue_type = auto_detect_tissue_type(path_to_db_file = custom_marker_file, 
                                          seuratObject = seurat_object, 
                                          scaled = scaled, assay = assay)
    rownames(tissue_type)=NULL
    tissue_type=tissue_type$tissue[1]
  } else {
    tissue_type = known_tissue_type
  }
  
  # Prepare gene sets
  gs_list = gene_sets_prepare(custom_marker_file, tissue_type)
  
  data_type <- if (scaled) "scale.data" else "counts"  
  package_type <- data_type %in% names(attributes(seurat_object[[assay]]))
  
  # Calculate scType scores
  if(package_type){
    
    print("Using Seurat v4 object")
    es.max = sctype_score(scRNAseqData = slot(seurat_object[[assay]], data_type),
                          scaled = TRUE,gs = gs_list$gs_positive, 
                          gs2 = gs_list$gs_negative)   
    
  } else{
    
    print("Using Seurat v5 object")
    
    if (data_type == "scale.data") {
      scRNAseqData <- seurat_object[[assay]]$scale.data
    } else {
      scRNAseqData <- seurat_object[[assay]]$counts
    }
    
    es.max = sctype_score(scRNAseqData = as.matrix(scRNAseqData),
                          scaled = TRUE,gs = gs_list$gs_positive, 
                          gs2 = gs_list$gs_negative)       
  }
  
  # Extract top cell types for each cluster
  cL_resutls = do.call("rbind", lapply(unique(seurat_object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_object@meta.data[seurat_object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  seurat_object_res=seurat_object
  seurat_object_res@meta.data[name] = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat_object_res@meta.data[seurat_object_res@meta.data$seurat_clusters == j,name] = as.character(cl_type$type[1])
  }
  if(plot){
    plot_ = DimPlot(seurat_object_res, reduction = "umap", group.by = name)   
    print(plot_)
  }
  text_=paste("New metadata added: ",name)
  print(text_)
  return(seurat_object_res)
} 


matrix <- GetAssayData(seurat_obj, layer = "counts")

sctype_source <- function(){
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  return(db_)
}



#Get database and tissuetype 
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- tissuetype # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

db_

tissue

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)



# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_obj[["RNA"]])));

#No need for this check
#print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seurat_obj[["RNA"]]$scale.data) else as.matrix(seurat_obj[["RNA"]]@scale.data)

es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)


cL_resutls <- do.call("rbind", lapply(unique(seurat_obj@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data$seurat_clusters==cl)), 10)
}))

sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  


sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"





#Get the celltype annotation to your seurat object
seurat_obj@meta.data$sctype_classification = ""

for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_obj@meta.data$sctype_classification[seurat_obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification')        


#Run ScType on your seurat object
seurat_obj <- run_sctype(seurat_obj, known_tissue_type = tissue, custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",name="sctype_classification",plot=TRUE)





# prepare edges
cL_resutls <- cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 <- sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss <- c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes <- rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db <- openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr <- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")


pdf(file = "DimPlot.pdf")

DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification') 
DimPlot(seurat_obj, reduction = "umap", label = F, repel = F, group.by = 'sctype_classification') 
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
print(gggr)
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)+gggr+DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification') 
DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_classification') + gggr



dev.off()

# EOF