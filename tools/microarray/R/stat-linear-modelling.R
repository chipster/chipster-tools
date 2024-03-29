# TOOL stat-linear-modelling.R: "Linear modelling" (Analyzes the data using linear modelling as implemented in the R package limma. You can have a maximum of three main effects and their interactions in the model. In addition, you can specify technical replication and pairing of the samples. Main effects can be fitted as such or as categorical variables, i.e., factors. Fold changes and p-values are reported for all effects and interactions.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT limma.tsv: limma.tsv
# OUTPUT limma-design.tsv: limma-design.tsv
# PARAMETER main.effect1: "Main effect 1" TYPE METACOLUMN_SEL DEFAULT group (Main effect 1)
# PARAMETER OPTIONAL main.effect2: "Main effect 2" TYPE METACOLUMN_SEL DEFAULT EMPTY (Main effect 2)
# PARAMETER OPTIONAL main.effect3: "Main effect 3" TYPE METACOLUMN_SEL DEFAULT EMPTY (Main effect 3)
# PARAMETER OPTIONAL technical.replication: "Technical replication" TYPE METACOLUMN_SEL DEFAULT EMPTY (Technical replication)
# PARAMETER OPTIONAL pairing: "Pairing" TYPE METACOLUMN_SEL DEFAULT EMPTY (Paired samples)
# PARAMETER OPTIONAL treat.main.effect1.as.factor: "Treat main effect 1 as factor" TYPE [no: no, yes: yes] DEFAULT no (Should main.effect1 be treated as a factor)
# PARAMETER OPTIONAL treat.main.effect2.as.factor: "Treat main effect 2 as factor" TYPE [no: no, yes: yes] DEFAULT no (Should main.effect2 be treated as a factor)
# PARAMETER OPTIONAL treat.main.effect3.as.factor: "Treat main effect 3 as factor" TYPE [no: no, yes: yes] DEFAULT no (Should main.effect3 be treated as a factor)
# PARAMETER OPTIONAL adjust.p.values: "Adjust p-values" TYPE [yes: yes, no: no] DEFAULT yes (Should the p-values be adjusted for multiple comparisons)
# PARAMETER OPTIONAL p.value.adjustment.method: "p-value adjustment method" TYPE [none: none, bonferroni: Bonferroni, holm: Holm, hochberg: Hochberg, BH: BH, BY: BY] DEFAULT BH (Multiple testing correction method)
# PARAMETER OPTIONAL interactions: "Interactions" TYPE [main: "main effects", two-way: "main effects and interactions"] DEFAULT main (What to include in the model)

# PARAMETER significance [main, interactions] DEFAULT main (Which p-values to return)

# JTT, 22.10.2007: Linear Modelling using limma
# MG, 28.03.2012, modified to handle NUID:s
# MK, 22.02.2014, modified so that does not report anymore coef or p-vals for pairing information
# OH, 12.02.2015, getting columns from phenodata using which rather than grep in order to get exact matches
# ML, 07.11.2016, simplify the outputs
# EK, 03.12.2020, specify optional parameters as required by Web App.

# main.effect1<-"group"
# main.effect2<-"gender"
# main.effect3<-"EMPTY"
# technical.replication<-"EMPTY"
# pairing<-"EMPTY"
# treat.main.effect1.as.factor<-"yes"
# treat.main.effect2.as.factor<-"yes"
# treat.main.effect3.as.factor<-"no"
# adjust.p.values<-"yes"
# p.value.adjustment.method<-"BH"
# interactions<-"two-way"
# significance<-"interactions"


# Loads the libraries
library(limma)

# Loads the normalized data
file <- c("normalized.tsv")
dat <- read.table(file, header = T, sep = "\t", row.names = 1)

# Separates expression values and flags
dat2 <- dat[, grep("chip", names(dat))]

# Loads phenodata
phenodata <- read.table("phenodata.tsv", header = T, sep = "\t")

# Sanity checks
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY") {
   stop("You need to specify at least one main effect! Please modify the setting accordingly, and rerun.")
}
if ((main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & interactions == "two-way") | (main.effect2 == "EMPTY" & interactions == "three-way") | (main.effect3 == "EMPTY" & interactions == "three-way")) {
   print("Only one main effect specified with interactions! No interactions specified for the model.")
   interactions <- c("main")
}
# if(interactions=="main" & significance=="interactions") {
#   print("Interactions can't be tested, since only main effect specified for the model!")
#   interactions<-c("main")
# }

# Extracting the variables from phenodata
main1 <- phenodata[, which(main.effect1 == colnames(phenodata))]
main2 <- phenodata[, which(main.effect2 == colnames(phenodata))]
main3 <- phenodata[, which(main.effect3 == colnames(phenodata))]
techrep <- phenodata[, which(technical.replication == colnames(phenodata))]
pair <- phenodata[, which(pairing == colnames(phenodata))]

# Converting vectors to factor, if needed
if (treat.main.effect1.as.factor == "yes") {
   main1 <- factor(main1)
}
if (treat.main.effect2.as.factor == "yes") {
   main2 <- factor(main2)
}
if (treat.main.effect3.as.factor == "yes") {
   main3 <- factor(main3)
}

# Specifying the models

# The basic models

# One main effect
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~main1)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~main2)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Two main effects
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main2)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main2 + main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Three main effects
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Two main effects and interactions
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main2 * main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Three main effects
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Three main effects and interactions
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2 * main3)
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}


# The basic models and biological replication

# Only technical replication
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# One main effect and technical replication
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~main1)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~main2)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Two main effects and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main2)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main2 + main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Three main effects and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Two main effects and interactions and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main2 * main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Three main effects and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Three main effects and interactions and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing == "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2 * main3)
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}


# Paired models

# One main effect
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main2 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Two main effects
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main2 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main2 + main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Three main effects
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Two main effects and interactions
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main2 * main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Three main effects
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# Three main effects and interactions
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication == "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2 * main3 + factor(pair))
   fit <- lmFit(dat2, design)
   fit <- eBayes(fit)
}

# The basic models and biological replication

# Only technical replication
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# One main effect and technical replication
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main2 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Two main effects and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main2 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main2 + main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Three main effects and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Two main effects and interactions and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 == "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 != "EMPTY" & main.effect2 == "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}
if (main.effect1 == "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main2 * main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Three main effects and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY") {
   design <- model.matrix(~ main1 + main2 + main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}

# Three main effects and interactions and technical replication
if (main.effect1 != "EMPTY" & main.effect2 != "EMPTY" & main.effect3 != "EMPTY" & technical.replication != "EMPTY" & pairing != "EMPTY" & interactions == "two-way") {
   design <- model.matrix(~ main1 * main2 * main3 + factor(pair))
   corfit <- duplicateCorrelation(dat2, ndups = 1, block = techrep)
   fit <- lmFit(dat2, design, block = techrep, cor = corfit$consensus)
   fit <- eBayes(fit)
}


# Extracting data
## New method ##
m <- matrix(nrow = nrow(dat2), ncol = ncol(design))
mm <- matrix(nrow = nrow(dat2), ncol = ncol(design))
for (i in 1:ncol(design)) {
   pp <- toptable(fit, coef = i, number = nrow(dat2), adjust.method = p.value.adjustment.method, sort.by = "none")
   if (adjust.p.values == "yes") {
      pp2 <- pp$adj.P.Val
   }
   if (adjust.p.values == "no") {
      pp2 <- pp$P.Value
   }
   pp3 <- pp$logFC
   m[, (i)] <- pp2
   mm[, (i)] <- pp3
   rownames(pp) <- rownames(dat2)
}


# Writing the data to disk
m <- data.frame(m)
mm <- data.frame(mm)

# Fold change
fc <- mm
colnames(fc) <- paste("FC.", colnames(design), sep = "")
fc2 <- mm
fc2 <- data.frame(round(fc2[, -1], digits = 2))
colnames(fc2) <- paste("chip.", colnames(design)[-1], sep = "")
rownames(fc2) <- rownames(dat)

# P-values
pvalues <- m
colnames(pvalues) <- paste("p.adjusted.", colnames(design), sep = "")
pvalues2 <- data.frame(pvalues)
pvalues2 <- data.frame(round(pvalues2[, -1], digits = 6))
colnames(pvalues2) <- paste("chip.p.adjusted.", colnames(design)[-1], sep = "")
rownames(pvalues2) <- rownames(dat)
pvaluesrounded <- round(pvalues, digits = 6)
fcrounded <- round(fc, digits = 2)

# Create data frames
dat3 <- cbind(dat, pvaluesrounded, fcrounded)

# Remove patient interactions
data_pairs <- grep("chip.factor(pair)", colnames(dat3))
if (length(data_pairs) > 0) {
   dat3 <- dat3[, -data_pairs]
}
data_pairs <- grep("chip.factor(pair)", colnames(fc2))
if (length(data_pairs) > 0) {
   fc2 <- fc2[, -data_pairs]
}
data_pairs <- grep("chip.factor(pair)", colnames(pvalues2))
if (length(data_pairs) > 0) {
   pvalues2 <- pvalues2[, -data_pairs]
}

# Write data to disk
write.table(dat3, file = "limma.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
write.table(design, file = "limma-design.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
# write.table(fc2, file="foldchange.tsv", sep="\t", row.names=T, col.names=T, quote=F)
# write.table(pvalues2, file="pvalues.tsv", sep="\t", row.names=T, col.names=T, quote=F)
