# TOOL metabarcoding-distance-ordination.R: "Distance matrices and ordinations" (Produces a Euclidean or Bray-Curtis distance matrix using a phyloseq object. Can also be used to compute Aitchison distances \(if selecting Euclidean distances for a CLR-transformed data set\). The matrix is used to perform either an 1\) unconstrained ordination \(non-metric multidimensional scaling, nMDS\) or 2\) a constrained ordination \(distance-based redundancy analysis, db-RDA\). To perform db-RDA, it is necessary to specify at least one phenodata variable describing the study design \(with the ordination axes constrained to linear combinations of variables\). It is possible to list up to five variables. Distances and a data frame containing sample variables are saved as an Rda file. Requires a phyloseq object in Rda format as the input. Note: groupings currently support only discrete \(i.e. non-numerical\) labels.)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT ps_ordi.txt: ps_ordi.txt
# OUTPUT ps_ordiplot.pdf: ps_ordiplot.pdf
# OUTPUT ps_dist.Rda: ps_dist.Rda
# PARAMETER disttype: "Type of distance measure" TYPE [euclidean: "Euclidean", bray: "Bray-Curtis"] DEFAULT euclidean (Choice between Euclidean and Bray-Curtis distances)
# PARAMETER orditype: "Type of ordination" TYPE [nmds: "nMDS", dbrda: "db-RDA"] DEFAULT nmds (Choice between using non-metric multidimensional scaling \(nMDS\) or distance-based redundancy analysis \(db-RDA\))
# PARAMETER samplevar: "Phenodata variable with sequencing sample IDs" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable with unique IDs for each community profile.)
# PARAMETER OPTIONAL samplenames: "Show sample IDs in ordination?" TYPE [yes: "Yes", no: "No"] DEFAULT yes (Should sample labels be plotted next to data points in the ordination?)
# PARAMETER OPTIONAL group_colour: "Phenodata variable for grouping ordination points by colour" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable used for grouping ordination points by colour.)
# PARAMETER OPTIONAL group_shape: "Phenodata variable for grouping ordination points by shape" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable used for grouping ordination points by shape.)
# PARAMETER OPTIONAL cap_variable1: "Phenodata variable 1 for db-RDA formula specification" TYPE METACOLUMN_SEL DEFAULT EMPTY (1st phenodata variable used in the \"formula\" argument when performing db-RDA \(minimum requirement is 1 variable\))
# PARAMETER OPTIONAL cap_variable2: "Phenodata variable 2 for db-RDA formula specification" TYPE METACOLUMN_SEL DEFAULT EMPTY (2nd phenodata variable used in the \"formula\" argument when performing db-RDA)
# PARAMETER OPTIONAL cap_variable3: "Phenodata variable 3 for db-RDA formula specification" TYPE METACOLUMN_SEL DEFAULT EMPTY (3rd phenodata variable used in the \"formula\" argument when performing db-RDA)
# PARAMETER OPTIONAL cap_variable4: "Phenodata variable 4 for db-RDA formula specification" TYPE METACOLUMN_SEL DEFAULT EMPTY (4th phenodata variable used in the \"formula\" argument when performing db-RDA)
# PARAMETER OPTIONAL cap_variable5: "Phenodata variable 5 for db-RDA formula specification" TYPE METACOLUMN_SEL DEFAULT EMPTY (5th phenodata variable used in the \"formula\" argument when performing db-RDA)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(ggplot2)
library(ggrepel)

# Load phyloseq object
load("ps.Rda")

# Stop if phenodata column 1 is not specified when using db-RDA
if (orditype == "dbrda" && 
	cap_variable1 == "EMPTY"){
	stop("db-RDA phenodata variable 1 not selected - at least one variable is required")
}

# nMDS or db-RDA 
# Also create two further objects for downstream statistics: 
# 1) a data frame of the ps object sample data
# 2) Euclidean or Bray-Curtis distances

if (disttype == "euclidean" && orditype == "nmds"){
	set.seed(1)
	ps_ordi <- ordinate(physeq = ps, 
		method = "NMDS", distance = "euclidean") # ordination
	ps_df <- data.frame(sample_data(ps)) # data frame
	set.seed(1)
	ps_dist <- phyloseq::distance(ps, method = 'euclidean') # distances
}
if (disttype == "euclidean" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 == "EMPTY" &&
		cap_variable3 == "EMPTY" &&
		cap_variable4 == "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "euclidean",
				formula = ~ get(cap_variable1))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'euclidean')
}
if (disttype == "euclidean" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 == "EMPTY" &&
		cap_variable4 == "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "euclidean",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'euclidean')
}
if (disttype == "euclidean" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 != "EMPTY" &&
		cap_variable4 == "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "euclidean",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2) +
						get(cap_variable3))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'euclidean')
}
if (disttype == "euclidean" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 != "EMPTY" &&
		cap_variable4 != "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "euclidean",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2) +
						get(cap_variable3) +
						get(cap_variable4))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'euclidean')
}
if (disttype == "euclidean" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 != "EMPTY" &&
		cap_variable4 != "EMPTY" &&
		cap_variable5 != "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "euclidean",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2) +
						get(cap_variable3) +
						get(cap_variable4) +
						get(cap_variable5))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'euclidean')
}
if (disttype == "bray" && orditype == "nmds"){
	set.seed(1)
	ps_ordi <- ordinate(physeq = ps, 
		method = "NMDS", distance = "bray")
	ps_df <- data.frame(sample_data(ps))
	set.seed(1)
	ps_dist <- phyloseq::distance(ps, method = 'bray')
}
if (disttype == "bray" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 == "EMPTY" &&
		cap_variable3 == "EMPTY" &&
		cap_variable4 == "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "bray",
				formula = ~ get(cap_variable1))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'bray')
}
if (disttype == "bray" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 == "EMPTY" &&
		cap_variable4 == "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "bray",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'bray')
}
if (disttype == "bray" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 != "EMPTY" &&
		cap_variable4 == "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "bray",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2) +
						get(cap_variable3))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'bray')
}
if (disttype == "bray" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 != "EMPTY" &&
		cap_variable4 != "EMPTY" &&
		cap_variable5 == "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "bray",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2) +
						get(cap_variable3) +
						get(cap_variable4))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'bray')
}
if (disttype == "bray" && orditype == "dbrda" &&
		cap_variable1 != "EMPTY" &&
		cap_variable2 != "EMPTY" &&
		cap_variable3 != "EMPTY" &&
		cap_variable4 != "EMPTY" &&
		cap_variable5 != "EMPTY"){
			set.seed(1)
			ps_ordi <- ordinate(physeq = ps, 
				method = "CAP", distance = "bray",
				formula = ~ get(cap_variable1) + 
						get(cap_variable2) +
						get(cap_variable3) +
						get(cap_variable4) +
						get(cap_variable5))
			ps_df <- data.frame(sample_data(ps))
			set.seed(1)
			ps_dist <- phyloseq::distance(ps, method = 'bray')
}

# Ordination plots

# Define a colour palette
# (Using a custom set containing up to 22 colours)
colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
		"#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
		"#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
		"#771122", "#AA4455", "#DD7788", "#CBD588")

if (group_colour != "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "no" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				color = group_colour) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("nMDS")
}
if (group_colour == "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "no" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("nMDS")
}
if (group_colour != "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "no" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi,
				color = group_colour, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("nMDS")
}
if (group_colour == "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "no" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi)+ 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("nMDS")
}
if (group_colour != "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "yes" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				color = group_colour) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("nMDS")
}
if (group_colour == "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "yes" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("nMDS")
}
if (group_colour != "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "yes" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi,
				color = group_colour, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("nMDS")
}
if (group_colour == "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "yes" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi)+ 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("nMDS")
}
if (group_colour != "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "no" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				color = group_colour) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("db-RDA")
}
if (group_colour == "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "no" &&
		orditype == "nmds"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("nMDS")
}
if (group_colour != "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "no" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi,
				color = group_colour, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("db-RDA")
}
if (group_colour == "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "no" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi)+ 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			ggtitle("db-RDA")
}
if (group_colour != "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "yes" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				color = group_colour) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("db-RDA")
}
if (group_colour == "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "yes" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("db-RDA")
}
if (group_colour != "EMPTY" &&
		group_shape != "EMPTY" &&
		samplenames == "yes" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi,
				color = group_colour, 
				shape = group_shape) + 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("db-RDA")
}
if (group_colour == "EMPTY" &&
		group_shape == "EMPTY" &&
		samplenames == "yes" &&
		orditype == "dbrda"){
		ordiplot <- plot_ordination(ps, ps_ordi)+ 
			theme(aspect.ratio = 1) +
			geom_point(size = 6.5, alpha = 0.75) +
			scale_fill_manual(values = colours) +
			theme(axis.title = element_blank()) +
			theme(axis.ticks = element_blank()) +
			theme(axis.text = element_blank()) +
			theme(legend.title = element_blank()) +
			theme(legend.position = "bottom") +
			geom_text_repel(aes(label = get(samplevar)),
                color = "gray45",
                min.segment.length = 0, 
                seed = 42, 
                box.padding = 0.6,
                size = 5) +
			ggtitle("db-RDA")
}

# Open a report PDF
pdf("ps_ordiplot.pdf")

ordiplot

# Close the report PDF
dev.off()

# Print out ordination summary
sink("ps_ordi.txt")
	cat("\n\n\n")
	cat("### Ordination summary ###\n")
	cat("\n\n\n")
	print(ps_ordi)
	cat("\n\n\n")
sink()

# Export Rda file
save(list = c("ps_df", "ps_dist"), 
	file = "ps_dist.Rda")
