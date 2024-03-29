# TOOL norm-affy-RPA.R: Affymetrix (Affymetrix preprocessing using CEL-files. Probe sets are automatically flagged using P A M flags. Variance stabilization can be applied only with MAS5 or Plier preprocessing methods. Custom chiptype can't be used with Plier preprocessing.)
# INPUT microarray{...}.cel: microarray{...}.cel TYPE AFFY
# OUTPUT normalized.tsv: normalized.tsv
# OUTPUT META phenodata.tsv: phenodata.tsv
# PARAMETER normalization.method: normalization.method TYPE [mas5: mas5, plier: plier, rma: rma, rpa: rpa, gcrma: gcrma, li-wong: li-wong] DEFAULT rma (Preprocessing method)
# PARAMETER stabilize.variance: stabilize.variance TYPE [yes: yes, no: no] DEFAULT no (Variance stabilazing normalization)
# PARAMETER custom.chiptype: custom.chiptype TYPE [empty: empty, hgu133ahsentrezg(hgu133a): hgu133ahsentrezg(hgu133a), hgu133a2hsentrezg(hgu133av2): hgu133a2hsentrezg(hgu133av2), hgu133phsentrezg(hgu133plus): hgu133phsentrezg(hgu133plus), hgu133plus2hsentrezg(hgu133plus2): hgu133plus2hsentrezg(hgu133plus2), hgu133bhsentrezg(hgu133b): hgu133bhsentrezg(hgu133b), hgu95av2hsentrezg(hgu95av2): hgu95av2hsentrezg(hgu95av2), moe430ammentrezg(moe430a): moe430ammentrezg(moe430a), moe430bmmentrezg(moe430b): moe430bmmentrezg(moe430b), mouse430a2mmentrezg(mouse430a2): mouse430a2mmentrezg(mouse430a2), mouse4302mmentrezg(mouse4302): mouse4302mmentrezg(mouse4302), mm74av1mmentrezg(mgu74a): mm74av1mmentrezg(mgu74a), mgu74av2mmentrezg(mgu74av2): mgu74av2mmentrezg(mgu74av2), mgu74bv2mmentrezg(mgu74bv2): mgu74bv2mmentrezg(mgu74bv2), mgu74cv2mmentrezg(mgu74cv2): mgu74cv2mmentrezg(mgu74cv2), rae230arnentrezg(rae230a): rae230arnentrezg(rae230a), rae230brnentrezg(rae230b): rae230brnentrezg(rae230b), rat2302rnentrezg(rat2302): rat2302rnentrezg(rat2302), rgu34arnentrezg(rgu34a): rgu34arnentrezg(rgu34a), rgu34brnentrezg(rgu34b): rgu34brnentrezg(rgu34b), rgu34crnentrezg(rgu34c): rgu34crnentrezg(rgu34c)] DEFAULT empty (custom chiptype)


# Affymetrix normalization
# JTT 8.6.2006
# Changes to column naming on 29.6.2006
# Changes to phenodata table writing on 29.1.2007
# Modified to work with R 2.9.0 on 12th May 2009
# Modifications by MG
# Changes to custom.chiptype PARAMETER to account for changes in naming of custom CDF packages (version 12), 23.9.2009
# Changes to cope with dropped custom package support for certain array types, 12.11.2009
# RPA method added 4.1.2011 (Leo Lahti)
# 25.6.2015 ML Moved the script to the "attick"

# Renaming variables
norm <- normalization.method
stabvar <- stabilize.variance


# Initializes analyses
library(affy)
library(gcrma)
library(RPA)

# Reads in data
dat <- ReadAffy()

# Modifies the data objects to take custom chiptype into account
if (custom.chiptype != "empty") {
   chiptype <- custom.chiptype
   library(Biostrings)
   # chiptype<-substr(x=chiptype, start=1, stop=(matchPattern("(", chiptype)@start-1))
   chiptype <- gsub("\\(.*\\)", "", chiptype, perl = T)
   dat@annotation <- chiptype
   dat@cdfName <- chiptype
} else {
   chiptype <- dat@annotation
}

# Since custom annotation packages are no longer supported for hgu133plu
# and for mm74av1 the packages for hgu133plu2 and for mm74av2 are going to
# be used instead
if (chiptype == "hgu133phsentrezg") {
   chiptype <- "hgu133plus2hsentrezg("
}
if (chiptype == "mm74av1mmentrezg") {
   chiptype <- "mm74av2mmentrezg"
}

# Normalizations (MAS5, RMA, GCRMA or Li-Wong (dChip))

# MAS5 normalization
if (norm == "mas5") {
   dat2 <- exprs(mas5(dat))
   calls <- exprs(mas5calls(dat))
   if (stabvar == "yes") {
      if (ncol(dat2) < 2) {
         stop("You need to have at least two chip to be able to use VSN!")
      }
      library(vsn)
      dat2 <- exprs(vsn(dat2))
   } else {
      dat2 <- log2(dat2)
   }
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}

# PLIER normalization
if (norm == "plier" & custom.chiptype == "empty") {
   library(plier)
   dat2 <- exprs(justPlier(eset = dat, replicate = 1:length(dat), get.affinities = FALSE, normalize = FALSE, norm.type = c("together"), augmentation = 0.1, defaultaffinity = 1.0, defaultconcentration = 1.0, attenuation = 0.005, seaconvergence = 0.000001, seaiteration = 3000, gmcutoff = 0.15, probepenalty = 0.001, concpenalty = 0.000001, usemm = TRUE, usemodel = FALSE, fitaffinity = T, plierconvergence = 0.000001, plieriteration = 3000, dropmax = 3.0, lambdalimit = 0.01, optimization = 0))
   calls <- exprs(mas5calls(dat))
   if (stabvar == "yes") {
      if (ncol(dat2) < 2) {
         stop("You need to have at least two chip to be able to use VSN!")
      }
      library(vsn)
      dat2 <- exprs(vsn(dat2))
   } else {
      dat2 <- log2(dat2)
   }
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}
if (norm == "plier" & custom.chiptype != "empty") {
   stop("Custom chipstypes can't be used with Plier! Use some other preprocessing method.")
}

# RMA normalization
if (norm == "rma" & custom.chiptype == "empty") {
   dat2 <- exprs(justRMA())
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}
if (norm == "rma" & custom.chiptype != "empty") {
   dat2 <- exprs(rma(dat))
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}

# RPA normalization
if (norm == "rpa" & custom.chiptype == "empty") {
   dat2 <- exprs(rpa(dat))
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}
if (norm == "rpa" & custom.chiptype != "empty") {
   dat2 <- exprs(rpa(dat))
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}

# GCRMA normalization
if (norm == "gcrma" & custom.chiptype == "empty") {
   dat2 <- exprs(justGCRMA(type = c("fullmodel"), fast = T, optimize.by = c("speed")))
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}
if (norm == "gcrma" & custom.chiptype != "empty") {
   dat2 <- exprs(gcrma(dat))
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}


# Li-Wong (dChip) normalization
if (norm == "li-wong") {
   dat2 <- exprs(expresso(dat, normalize.method = "invariantset", bg.correct = FALSE, pmcorrect.method = "pmonly", summary.method = "liwong"))
   calls <- exprs(mas5calls(dat))
   dat2 <- as.data.frame(round(dat2, digits = 2))
   calls <- as.data.frame(calls)
   names(dat2) <- paste("chip.", names(dat2), sep = "")
   names(calls) <- paste("flag.", names(calls), sep = "")
   dat2 <- data.frame(dat2, calls)
}

# Writes out a phenodata table
sample <- rownames(pData(dat))
group <- c(rep("", nrow(pData(dat))))
training <- c(rep("", nrow(pData(dat))))
time <- c(rep("", nrow(pData(dat))))
random <- c(rep("", nrow(pData(dat))))
chiptype <- paste(chiptype, ".db", sep = "")
write.table(data.frame(sample = sample, chiptype = chiptype, group = group), file = "phenodata.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

# Writing out data
a <- try(library(chiptype, character.only = T))
if (chiptype != "empty" & class(a) != "try-error") {
   # Including gene names to data
   lib2 <- sub(".db", "", chiptype)
   symbol <- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep = "")))))[rownames(dat2), ])
   genename <- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep = "")))))[rownames(dat2), ])
   # Fxes an issue introduced in BioC2.4 where the "#" character is introduced in some gene names
   genename <- gsub("#", "", genename)
   symbol <- gsub("'", "", symbol)
   genename <- gsub("'", "", genename)
   # Writes the results into a file
   write.table(data.frame(symbol, description = genename, dat2), file = "normalized.tsv", col.names = T, quote = F, sep = "\t", row.names = T)
}

if (chiptype == "empty" | class(a) == "try-error") {
   write.table(data.frame(dat2), file = "normalized.tsv", col.names = T, quote = F, sep = "\t", row.names = T)
}
