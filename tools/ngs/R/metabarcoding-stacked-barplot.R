# TOOL metabarcoding-stacked-barplot.R: "OTU relative abundance bar plots" (This tool is used to visualize OTU relative abundance data within a phyloseq object using ggplot2. The tool produces stacked bar plots that are built using a user-specified level of biological organization \(with a maximum of displaying 22 taxa\). Options may be specified for faceting \(up to two phenodata variables\) and using a minimal OTU relative abundance threshold. Requires a phyloseq object in Rda format as the input. To use this tool, OTU abundances must be provided as relative abundances \(with per-column OTU abundances adding up to 1\).)
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT ps_barplot.pdf
# PARAMETER xvar: "Phenodata variable with sequencing sample IDs" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata variable with unique IDs for each community profile.)
# PARAMETER OPTIONAL thresh: "Relative abundance cut-off threshold \(%\) for excluding OTUs" TYPE INTEGER FROM 0 TO 99 DEFAULT 0
# PARAMETER OPTIONAL type: "Level of biological organization for tabulating taxon composition" TYPE [phylum: "Phylum", class: "Class", order: "Order", family: "Family", genus: "Genus"] DEFAULT phylum (Level of biological organization for plot construction \(default is phylum\))
# PARAMETER OPTIONAL facet1: "Phenodata variable 1 for plot faceting" TYPE METACOLUMN_SEL DEFAULT EMPTY (First phenodata variable used for plot faceting.)
# PARAMETER OPTIONAL facet2: "Phenodata variable 2 for plot faceting" TYPE METACOLUMN_SEL DEFAULT EMPTY (Second phenodata variable used for plot faceting.)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(phyloseq)
library(ggplot2)

# Load phyloseq object
load("ps.Rda")

# Check if per-column OTU counts add up to 1; if not, stop
if (any(sample_sums(ps) != 1)) {
	stop("CHIPSTER-NOTE: Per-column counts do not add up to 1, please convert your data to relative abundances.")
}

# Agglomerate data at the desired taxonomic level
# Functions also include some (likely superfluous) data cleaning
if (type == "phylum"){
	ps <- tax_glom(ps, taxrank = rank_names(ps)[2],
		NArm = TRUE, bad_empty = c(NA,""," ","\t"))
	} else if (type == "class"){
		ps <- tax_glom(ps, taxrank = rank_names(ps)[3],
			NArm = TRUE, bad_empty = c(NA,""," ","\t"))
	} else if (type == "order"){
		ps <- tax_glom(ps, taxrank = rank_names(ps)[4],
			NArm = TRUE, bad_empty = c(NA,""," ","\t"))
	} else if (type == "family"){
		ps <- tax_glom(ps, taxrank = rank_names(ps)[5],
			NArm = TRUE, bad_empty = c(NA,""," ","\t"))
	} else if (type == "genus") {
		ps <- tax_glom(ps, taxrank = rank_names(ps)[6],
			NArm = TRUE, bad_empty = c(NA,""," ","\t"))
}

# Melt to long format (also converts ps object to data frame)
ps_df <- psmelt(ps)

# Exclude OTUs under user-specified relative abundance threshold
if (thresh != "0"){
	thresh <- thresh / 100 # Divide by 100 to get correct format (e.g. 3% -> 0.03)
	ps_df <- subset(ps_df, Abundance > thresh)
}

# Define a colour palette
# (Using a custom set containing up to 22 colours)
colours <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", 
		"#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", 
		"#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
		"#771122", "#AA4455", "#DD7788", "#CBD588")

# Bar plots

if (facet1 == "EMPTY" &&
	facet2 != "EMPTY"){
	stop("Facet 2 specified without facet 1; if faceting by a single variable, use facet 1")
}
if (type == "phylum" && 
	facet1 != "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Phylum)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_wrap(~ get(facet1), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "class" && 
	facet1 != "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Class)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_wrap(~ get(facet1), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "order" && 
	facet1 != "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Order)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_wrap(~ get(facet1), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "family" && 
	facet1 != "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Family)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_wrap(~ get(facet1), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "genus" && 
	facet1 != "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Genus)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_wrap(~ get(facet1), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "phylum" && 
	facet1 != "EMPTY" &&
	facet2 != "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Phylum)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_grid(get(facet1) ~ get(facet2), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "class" && 
	facet1 != "EMPTY" &&
	facet2 != "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Class)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_grid(get(facet1) ~ get(facet2), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "order" && 
	facet1 != "EMPTY" &&
	facet2 != "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Order)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_grid(get(facet1) ~ get(facet2), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "family" && 
	facet1 != "EMPTY" &&
	facet2 != "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Family)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_grid(get(facet1) ~ get(facet2), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 8))
}
if (type == "genus" && 
	facet1 != "EMPTY" &&
	facet2 != "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Genus)) + 
	geom_bar(stat = "identity", position = "stack") +
	facet_grid(get(facet1) ~ get(facet2), scales = "free") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "phylum" && 
	facet1 == "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Phylum)) + 
	geom_bar(stat = "identity", position = "stack") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "class" && 
	facet1 == "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Class)) + 
	geom_bar(stat = "identity", position = "stack") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "order" && 
	facet1 == "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Order)) + 
	geom_bar(stat = "identity", position = "stack") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "family" && 
	facet1 == "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Family)) + 
	geom_bar(stat = "identity", position = "stack") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}
if (type == "genus" && 
	facet1 == "EMPTY" &&
	facet2 == "EMPTY"){
psplot <- ggplot(ps_df, aes(x = get(xvar), 
		y = Abundance, fill = Genus)) + 
	geom_bar(stat = "identity", position = "stack") +
	scale_fill_manual(values = colours) +
	scale_y_continuous(name = "Relative abundance (%)\n") +
	theme(axis.text.y = element_text(size = 12, hjust = 0.5)) +
	theme(axis.title.y = element_text(size = 14)) +
	theme(axis.title.x = element_blank()) +
	theme(axis.text.x = element_text(size = 12, 
		angle = 30, vjust = 0.5, hjust = 0.5)) +
	theme(axis.ticks.x = element_blank()) +
	theme(strip.text.x = element_text(size = 11)) +
	theme(legend.title = element_blank()) +
	theme(legend.position = "right") +
	theme(legend.text = element_text(size = 10))
}

# Open a report PDF
pdf("ps_barplot.pdf", width = 12, height = 7)

psplot

# Close the report PDF
dev.off()
