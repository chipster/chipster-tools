# TOOL metabarcoding-rarecurve-alpha.R: "Sequence numbers, rarefaction curve and alpha diversity estimates" (Lists per-sample sequence numbers, plots rarefaction curves and tabulates alpha diversity estimates \(observed no. of OTUs, Chao1 and Shannon indices, and Pielou's evenness\) and calculates means of them if the group parameter is selected. It calculates also Wilcoxon rank sum test between two groups. If you have more than two groups, please specify which groups you want to compare. Requires a phyloseq object in Rda format as the input. Note that the diversity estimates are only reliable when using raw \(untrimmed\) OTU data, as many diversity metrics are dependent on singletons and doubletons in the data set under analysis. )
# INPUT ps.Rda: "Phyloseq object in Rda format" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata" TYPE GENERIC
# OUTPUT ps_rarecurve.pdf
# OUTPUT ps_alphadiv.txt
# PARAMETER OPTIONAL group_column: "Phenodata variable for showing grouping" TYPE METACOLUMN_SEL DEFAULT empty (Phenodata variable describing groping which is added to alpha diversity table for improved readability and for calculating means within groups.)
# PARAMETER OPTIONAL group1: "Group 1 for Wilcoxon rank sum test  (if >2 groups overall)" TYPE STRING DEFAULT empty (First sample group name (one of the sample groups under the phenodata variable used for grouping\) for Wilcoxon rank sum test  )
# PARAMETER OPTIONAL group2: "Group 2 for Wilcoxon rank sum test  (if >2 groups overall)" TYPE STRING DEFAULT empty (Second sample group name (one of the sample groups under the phenodata variable used for grouping\) for Wilcoxon rank sum test  )
# RUNTIME R-3.6.1-phyloseq

# JH 2020
# ES 9.7.2021 alpha diversity estimates for rarefied data
# ES 9.8.2021 calculate means, standard deviation, standard error, wilcox rank sum test
# PARAMETER OPTIONAL type: "Are the data raw (untrimmed) or rarefied to even depth?" TYPE [raw, rarefied] DEFAULT raw (Type of data to be analyzed \(raw vs rarefied; by default, assumes raw data.\))

# Load libraries
library(microbiome)
library(phyloseq)
library(vegan)

# Load data
load("ps.Rda")

# Calculate per-sample sequence no.s
seqno <- sample_sums(ps)

# Alpha diversity estimates and Pielou's evenness
if (group_column == "empty") {
	set.seed(1)
	richness <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
	set.seed(1)
	pielou <- evenness(ps, 'pielou')
	richness <- cbind(richness, pielou)
}
if (group_column != "empty") {
	set.seed(1)
	richness <- estimate_richness(ps, measures = c("Observed", "Chao1", "Shannon"))
	set.seed(1)
	pielou <- evenness(ps, 'pielou')
	richness <- cbind(richness, pielou)	
	colno <- pmatch(group_column, colnames(ps@sam_data)) # Get desired column no. from ps metadata
	richness <- cbind(richness, ps@sam_data[,colno]) # cbind the column to the diversity table
}

# Open a report PDF
pdf("ps_rarecurve.pdf")

# Plot rarefaction curve
set.seed(1)
rarecurve(t(otu_table(ps)), step = 100, 
          cex.lab = 1.5, cex.axis = 1.5, label = FALSE, ylab = "OTUs", xlab = "No. of sequences")

# Close the report PDF
dev.off()
# test variable for goup1 and group2 column
test1 <- FALSE
test2 <- FALSE
# Print out sequence numbers and alpha diversity table
#if (type == "raw"){
sink("ps_alphadiv.txt")
	cat("\n\n\n")
	cat("### Per-sample sequence no.s ###\n")
	cat("\n\n\n")
	print(seqno)
	cat("\n\n\n")
	cat("### Alpha diversity estimates (observed OTUs, Chao1, Shannon's index, Pielou's evenness) ###\n")
	cat("\n\n\n")
	print(richness)
	cat("\n\n\n")

	# calculate means standard error and Wilcoxon rank sum test by ES 
	if (group_column != "empty") {
		cat("### Means, standard deviations and standard errors from alpha diversity estimates by every group specified with the phenodata parameter")
		cat("\n\n")
		
		# Loads phenodata and list all different groups to groupnames
		phenodata<-read.table("phenodata.tsv", header=T, sep="\t")
		groupnames <- unique(phenodata[,group_column])

		
		x<-1
		# calculate means, sd and std.error for every group
		for (name in groupnames){
			cat("Group name:",name, sep = " ")
			cat("\n")
			# take those rows where group == name / last column
			set1 <- subset(richness, richness[,ncol(richness)] == name)
			# Observed OTUs
			
			cat("Observed OTUs\t\tMean:",round(mean(set1[,"Observed"]),3),sep= " ")
			cat("\t\tStandard deviation:",round(sd(set1[,"Observed"]),3),sep=" ")
			cat("\t\tStandard Error:",round(sd(set1[,"Observed"])/sqrt(length(set1[,"Observed"])),3), sep = " ")
			cat("\n")
			# Chao1 
			chao1_m <-mean(set1[,"Chao1"])
			chao1_std <- sd(set1[,"Chao1"])/sqrt(length(set1[,"Chao1"]))
			cat("Chao1\t\t\tMean:",round(chao1_m,3), sep = " ")
			cat("\t\tStandard deviation:",round(sd(set1[,"Chao1"]),3),sep=" ")
			cat("\t\tStandard Error:",round(chao1_std,3), sep = " ")
			cat("\n")
			#sqrt(sum(!is.na(x)))
			# Shannon index
			shannon_m <- mean(set1[,"Shannon"])
			cat("Shannon's index\t\tMean:",round(shannon_m,3), sep = " ")
			cat("\t\tStandard deviation:",round(sd(set1[,"Shannon"]),3),sep=" ")
			cat("\t\tStandard Error:",round(sd(set1[,"Shannon"])/sqrt(length(set1[,"Shannon"])),3), sep = " ")
			cat("\n")
			# Pielou evennes
			pielou_m <- mean(set1[,"pielou"])
			cat("Pielou evenness\t\tMean:",round(pielou_m,3), sep =" ")
			cat("\t\tStandard deviation:",round(sd(set1[,"pielou"]),3),sep=" ")
			cat("\t\tStandard Error:",round(sd(set1[,"pielou"])/sqrt(length(set1[,"pielou"])),3), sep = " ")
			cat("\n\n")
			# if more than 2 groups use parameters group1 and group2
			if (length(groupnames)>2){
				if(group1==name){
					test1 <- TRUE
					chao1_a <- set1[,"Chao1"]
					shannon_a <- set1[,"Shannon"]
					pielou_a <- set1[,"pielou"]
					observed_a <- set1[,"Observed"]
				}
				if (group2 == name){
					test2 <- TRUE
					chao1_b <- set1[,"Chao1"]
					shannon_b <- set1[,"Shannon"]
					pielou_b <- set1[,"pielou"]
					observed_b <- set1[,"Observed"]
				}
				
			}else
				if (x==1){
					observed_a <- set1[,"Observed"]
					chao1_a <- set1[,"Chao1"]
					shannon_a <- set1[,"Shannon"]
					pielou_a <- set1[,"pielou"]
				}else{
				observed_b <- set1[,"Observed"]
				chao1_b <- set1[,"Chao1"]
				shannon_b <- set1[,"Shannon"]
				pielou_b <- set1[,"pielou"]
			}
			x = x+1
			}
		if (length(groupnames)<2) {
			cat("### You have only one group. You need at least 2 groups for Wilcoxon rank sum test ###")
		}else if (length(groupnames)==2){
			names = levels(groupnames)
			cat("### Wilcoxon rank sum test between groups:",names[1],names[2],"###", sep = " ")
			cat("\n\n")

			cat("# For observed OTUs: \n")
			print(wilcox.test(observed_a,observed_b))
		
			cat("# For Chao1: \n")
			print(wilcox.test(chao1_a,chao1_b))
			
			cat("# For Shannon index: \n")
			print(wilcox.test(shannon_a,shannon_b))
		
			cat("# For Pielou's evennes: \n")
			print(wilcox.test(pielou_a,pielou_b))
		}else{
			if (test1==FALSE || test2 == FALSE){
				cat("###You have more than 2 groups. Couldn't find groups",group1,",",group2," from Phenodata file.### \n", sep=" ")
				cat("Give the right groupnames as parameters.")
			}else{
				cat("### Wilcoxon rank sum test between groups:",group1,",",group2, sep = " ")
				cat("\n\n")
				cat("# For observed OTUs: \n")
				print(wilcox.test(observed_a,observed_b))
				cat("# For Chao1: \n")
				print(wilcox.test(chao1_a,chao1_b))
				cat("# For Shannon index: \n")
				print(wilcox.test(shannon_a,shannon_b))
				cat("# For Pielou's evennes: \n")
				print(wilcox.test(pielou_a,pielou_b))
			}
		}
		}else{
			cat("### For calculating means and Wilcoxon rank sum test, you need to select the groupcolumn in the parameter tab.###")
		}
	cat("\n\n")
sink()


#}
# commented out by ES
#if (type == "rarefied"){
#sink("ps_alphadiv.txt")
#	cat("\n\n\n")
#	cat("### Per-sample sequence no.s ###\n")
#	cat("\n\n\n")
#	print(seqno)
#	cat("\n\n\n")
#	cat("### Alpha diversity estimates (observed OTUs, Chao1, Shannon's index, Pielou's evenness) ###\n")
#	cat("\n\n\n")
#	cat("N/A (unavailable due to use of rarefied data)\n")
#	cat("\n\n\n")
#sink()
#}



