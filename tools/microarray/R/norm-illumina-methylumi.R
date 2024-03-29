# TOOL norm-illumina-methylumi.R: "Illumina - methylumi pipeline" (Illumina methylation assay normalization using FinalReport files and lumi methodology. You need to import the FinalReport file DIRECTLY, not using the Import tool. The FinalReport file has to contain sample methylation profiles, group profile will not work with this tool.)
# INPUT FinalReport_sample_methylation_profile.txt: FinalReport_sample_methylation_profile.txt TYPE GENERIC
# OUTPUT normalized.tsv: normalized.tsv
# OUTPUT unmethylated.tsv: unmethylated.tsv
# OUTPUT methylated.tsv: methylated.tsv
# OUTPUT META phenodata.tsv: phenodata.tsv
# OUTPUT OPTIONAL QC-plot.pdf: QC-plot.pdf
# PARAMETER chiptype: Chiptype TYPE [FDb.InfiniumMethylation.hg19:FDb.InfiniumMethylation.hg19] DEFAULT FDb.InfiniumMethylation.hg19 (Select the annotation package)
# PARAMETER target: "Illumina ID" TYPE STRING DEFAULT TargetID (Name of the TargetID column)
# PARAMETER signala: "Signal_A/Grn pattern" TYPE STRING DEFAULT Signal_A (Pattern identifying and common to all unmethylated data columns)
# PARAMETER signalb: "Signal_B/Red pattern" TYPE STRING DEFAULT Signal_B (Pattern identifying and common to all methylated data columns)
# PARAMETER OPTIONAL normalization: Normalization TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile ()
# PARAMETER OPTIONAL color.balance.adjustment: "Color balance adjustment" TYPE [none: none, quantile: quantile, ssn: ssn] DEFAULT quantile (Adjustment of color balance)
# PARAMETER OPTIONAL background.correction: "Background correction" TYPE [none: none, bgAdjust2C: bgAdjust2C, forcePositive: forcePositive] DEFAULT none (Should background adjustment be applied)
# PARAMETER OPTIONAL QCplots: QCplots TYPE [yes: yes, no: no] DEFAULT yes (Do you want quality control plots)
# PARAMETER OPTIONAL image.width: "Image width" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Width of the QC image)
# PARAMETER OPTIONAL image.height: "Image height" TYPE INTEGER FROM 200 TO 3200 DEFAULT 600 (Height of the QC image)

# JTT: 02.02.2011: Illumina methylation array data preprocessing and normalization for FinalReport file
# JTT: 25.09.2012: Modified
# JTT: 28.10.2012: Modified
# MK: 20.11.2013: Added parameters to enalbe analysis of data having atypical columnames
# TH: 10.11.2016: Disable chiptype HumanMethylation450: HumanMethylation450 for now
# ML: 27.2.2017: Fix, use FDb.InfiniumMethylation.hg19 package. Disable HumanMethylation27: HumanMethylation27 package.

# setwd("C://Users//Jarno Tuimala//Desktop//methylumi data")
# color.balance.adjustment<-c("quantile")
# background.correction<-c("none")
# normalization<-c("quantile")
# chiptype<-c("HumanMethylation450", "FDb.InfiniumMethylation.hg19")
# QCplots<-"yes"
# image.width<-c(600)
# image.height<-c(600)

# Renaming variables
w <- image.width
h <- image.height

system(paste("perl -p -i -e \'s/", target, "/TargetID/g\' FinalReport_sample_methylation_profile.txt", sep = ""))
system(paste("perl -p -i -e \'s/", signala, "/Signal_A/g\' FinalReport_sample_methylation_profile.txt", sep = ""))
system(paste("perl -p -i -e \'s/", signalb, "/Signal_B/g\' FinalReport_sample_methylation_profile.txt", sep = ""))

# Loading libraries
library(lumi)
library(methylumi)
library(annotate)

# Converting to the correct chiptype
# if(chiptype=="HumanMethylation27") {
# 	chiptype<-c("IlluminaHumanMethylation27k")
# 	chiptype<-paste(chiptype, ".db", sep="")
# }
# if(chiptype=="HumanMethylation450") {
# 	chiptype<-c("IlluminaHumanMethylation450k")
# 	chiptype<-paste(chiptype, ".db", sep="")
# }
if (chiptype == "FDb.InfiniumMethylation.hg19") {
   chiptype <- "FDb.InfiniumMethylation.hg19"
}

# Loading data files.
# Note that if the pattern given in target-variable, following wrong error is returned:
# in .getFileSeparator(filename) :
#   Cannot determine separator used in the file, please manually set the "sep" parameter!
# Calls: methylumiR -> .getFileSeparator

dat <- methylumiR("FinalReport_sample_methylation_profile.txt", lib = chiptype, sep = "\t")
methyLumiM <- as(dat, "MethyLumiM")
methyLumiM <- addAnnotationInfo(methyLumiM, lib = chiptype)

# Color balance adjustment
dat2 <- lumiMethyC(methyLumiM, method = color.balance.adjustment)

# Background adjustment
if (color.balance.adjustment != "none") {
   dat3 <- lumiMethyB(dat2, method = background.correction, separateColor = FALSE)
} else {
   dat3 <- lumiMethyB(dat2, method = background.correction, separateColor = TRUE)
}

# Normalization
dat4 <- lumiMethyN(dat3, method = normalization)
# dat4<-normalizeMethyLumiSet(dat)

# QC plots
if (QCplots == "yes") {
   pdf(file = "QC-plot.pdf", width = w / 72, height = h / 72)
   par(mfrow = c(2, 2))
   plotColorBias1D(dat, main = "Unpreprocessed")
   plotColorBias1D(dat4, main = "Preprocessed")
   boxplotColorBias(dat, main = "Unpreprocessed")
   boxplotColorBias(dat4, main = "Preprocessed")
   dev.off()
}

# Convert sample names to Chipster style
dat5 <- exprs(dat4)
sample.names <- colnames(dat5)
sample.names <- paste("chip.", sample.names, sep = "")
names(dat5) <- sample.names
colnames(dat5) <- sample.names

dat6 <- methylated(dat4)
sample.names <- colnames(dat6)
sample.names <- paste("chip.", sample.names, sep = "")
names(dat6) <- sample.names
colnames(dat6) <- sample.names

dat7 <- unmethylated(dat4)
sample.names <- colnames(dat7)
sample.names <- paste("chip.", sample.names, sep = "")
names(dat7) <- sample.names
colnames(dat7) <- sample.names

# Annotations
# Not up to date
# library(chiptype, character.only=T)
# symbols <- unlist (lookUp(rownames(dat5), chiptype, what="SYMBOL"))
# genenames <- unlist (lookUp(rownames(dat5), chiptype, what="GENENAME"))
# symbols <- gsub("#", "", symbols)
# genenames <- gsub("#", "", genenames)
# symbols <- gsub("'", "", symbols)
# genenames <- gsub("'", "", genenames)

# Write out a phenodata
group <- c(rep("", ncol(dat5)))
training <- c(rep("", ncol(dat5)))
write.table(data.frame(sample = sample.names, chiptype = chiptype, group = group), file = "phenodata.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

# Write out expression data
write.table(data.frame(dat5), file = "normalized.tsv", col.names = T, quote = F, sep = "\t", row.names = T)
write.table(data.frame(dat6), file = "methylated.tsv", col.names = T, quote = F, sep = "\t", row.names = T)
write.table(data.frame(dat7), file = "unmethylated.tsv", col.names = T, quote = F, sep = "\t", row.names = T)

# write.table(data.frame(symbol=symbols, description=genenames, dat5), file="normalized.tsv", col.names=T, quote=F, sep="\t", row.names=T)
# write.table(data.frame(symbol=symbols, description=genenames, dat6), file="methylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
# write.table(data.frame(symbol=symbols, description=genenames, dat7), file="unmethylated.tsv", col.names=T, quote=F, sep="\t", row.names=T)
