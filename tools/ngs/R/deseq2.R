# TOOL deseq2.R: "Differential expression using DESeq2" (Differential expression analysis using the DESeq2 Bioconductor package. You can create the input count table and phenodata file using the tool \"Utilities - Define NGS experiment\". If you have more than two experimental groups, note that the output figures sum up information from all pairwise comparisons.)
# INPUT data.tsv: "Count table" TYPE GENERIC
# INPUT META phenodata.tsv: "Phenodata file" TYPE GENERIC
# OUTPUT OPTIONAL de-list-deseq2.tsv
# OUTPUT OPTIONAL summary.txt
# OUTPUT OPTIONAL deseq2_report.pdf
# OUTPUT OPTIONAL de-list-deseq2.bed
# PARAMETER column: "Column describing groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test.)
# PARAMETER OPTIONAL ad_factor: "Column describing additional experimental factor" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing an additional experimental factor. If given, p-values in the output table are from a likelihood ratio test of a model including the experimental groups and experimental factor, vs a model which only includes the experimental factor.)
# PARAMETER OPTIONAL p.value.cutoff: "Cutoff for the adjusted P-value" TYPE DECIMAL FROM 0 TO 1 DEFAULT 0.05 (The cutoff for Benjamini-Hochberg adjusted p-value. Note that the developers of DESeq2 use 0.1 as a default cut-off.)
# PARAMETER OPTIONAL round_pvals: "Round p-values to this many digits" TYPE INTEGER DEFAULT 4 (Round p-values and adjusted p-values to this many digits.)
# PARAMETER OPTIONAL round_others: "Round other values to this many digits" TYPE INTEGER DEFAULT 2 (Round other computed values to this many digits.)
# PARAMETER OPTIONAL bed: "Create BED file" TYPE [yes,no] DEFAULT no (Create a BED file.)
# PARAMETER OPTIONAL size.factor.estimation.type: "Method for estimating size factors" TYPE [ratio:"ratio", poscounts:"poscounts"] DEFAULT ratio (Which size estimation method to use. Default option "ratio" uses the standard median ratio method introduced in DESeq. "poscounts" offers alternative estimator, which can be used even when all genes contain a sample with a zero.)
# RUNTIME R-4.2.0-phyloseq



## OLD:  RUNTIME R-3.6.1-phyloseq
# For testing:
# OUTPUT OPTIONAL dds.Robj
# OUTPUT OPTIONAL output_table.tsv
# OUTPUT OPTIONAL res.Robj

# MK 15.04.2014, added the possibility to use DESeq2 in dea-deseq.R
# AMS 17.06.2014, split the DESeq2 part to a separate tool
# EK 1.7.2014, clarified the script before moving it to production, and fixed a bug that disabled DESeq2's automatic independent filtering
# EK 9.2.2015, updated to R3.1.2, changed the MA plot, added summary
# AMS 7.4.2015, Join pdf outputs to one
# ML 25.6.2015, Fixed some problems with S4vectors
# EK 10.2.2016, Added alpha to results()
# ML 4.5.2016, Changed the direction of the comparison when more than 2 groups
# ML+SS 18.10.2016, Fixed plotting & rounding and formatting when more than 2 groups
# ML+SS 08.05.2018, comparison possible also when >9 groups
# ML 10.5.2022, Move to R-3.6.1 and DESeq 1.26.0, add the size factor estimation type parameter to cope with cases when all genes contain a sample with a zero.
# ML 02.02.2023, Add the LFC shrink (this was previously inside DESeq function)
# ML 26.02.2024, Add rounding parameters

# column <-"group"
# ad_factor<-"EMPTY"
# p.value.cutoff<-0.05


# Load the library
source(file.path(chipster.common.lib.path, "bed-utils.R"))
library(DESeq2)

# Load the counts data and extract expression value columns
dat <- read.table("data.tsv", header = T, sep = "\t", row.names = 1)
dat2 <- dat[, grep("chip", names(dat))]

# Get the experimental group information from phenodata
phenodata <- read.table("phenodata.tsv", header = T, sep = "\t")
groups <- as.character(phenodata[, pmatch(column, colnames(phenodata))])
group_levels <- levels(as.factor(groups))

# Read the additional experimental factor from phenodata and construct a design matrix from this information
exp_factor <- NULL
if (ad_factor != "EMPTY") {
    exp_factor <- as.character(phenodata[, pmatch(ad_factor, colnames(phenodata))])
    design <- data.frame(condition = as.factor(groups), exp_factor = exp_factor)
    rownames(design) <- colnames(dat2)
}

# Create a DESeqDataSet object
if (ad_factor == "EMPTY") {
    dds <- DESeqDataSetFromMatrix(countData = dat2, colData = data.frame(condition = groups), design = ~condition)
} else if (ad_factor != "EMPTY") {
    dds <- DESeqDataSetFromMatrix(countData = dat2, colData = design, design = ~ exp_factor + condition)
}

# Vector / variable that holds comparison names
results_name <- NULL

dds <- DESeq(dds, sfType = size.factor.estimation.type)

# save(dds, file="dds.Robj")


# Calculate statistic for differential expression, merge with original data table, keep significant DEGs, remove NAs and sort by FDR.
# If there are more than 2 groups, get pairwise results for each comparison.

# Two groups:
if (length(unique(groups)) == 2) {
    res <- results(dds, alpha = p.value.cutoff)

    # LFC shrink, needs to be after results function. (This step was previously inside DESeq function)
    resLFC <- lfcShrink(dds, res = res, type = "ashr")
    res <- resLFC

    sig <- cbind(dat, res)
    sig <- as.data.frame(sig)
    # sig <- sig[! (is.na(sig$padj)), ]
    sig <- sig[sig$padj <= p.value.cutoff, ]
    sig <- sig[!(is.na(sig$padj)), ]
    sig <- sig[order(sig$padj), ]
    # Open pdf file for output
    pdf(file = "deseq2_report.pdf")

    ## MA plot before lfcShrink:
    # plotMA(dds, alpha=p.value.cutoff, main=c("DESeq2 MA-plot, FDR =", p.value.cutoff), ylim=c(-2,2))

    # MA plot after lfcShrink:
    plotMA(resLFC, alpha = p.value.cutoff, main = c("DESeq2 MA-plot with LFC shrink (used in result table), FDR =", p.value.cutoff), ylim = c(-2, 2))

    # Summary into a txt file
    sink("summary.txt")
    summary(res, alpha = p.value.cutoff)
    sink()

    # Output significant DEGs
    if (dim(sig)[1] > 0) {
        ndat <- ncol(dat)
        nmax <- ncol(sig)
        write.table(cbind(sig[, 1:ndat], round(sig[, (ndat + 1):(nmax - 2)], digits = round_pvals), format(sig[, (nmax - 1):nmax], digits = round_others, scientific = T)), file = "de-list-deseq2.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
    }

    # write.table(sig, file = "output_table.tsv", sep = "\t", row.names = T, col.names = T, quote = F)


    # More than 2 groups:
} else if (length(unique(groups)) > 2) {
    test_results <- dds

    # Open pdf file for output
    pdf(file = "deseq2_report.pdf")

    ## MA plot before lfcShrink:
    # plotMA(dds,alpha=p.value.cutoff,main=c("DESeq2 MA-plot, FDR =", p.value.cutoff),ylim=c(-2,2))

    res <- NULL
    # going through all the pairwise comparisons (i vs j):
    conditions <- colData(test_results)$condition
    # i goes from the first group to the one before last:
    for (i in 1:(nlevels(conditions) - 1)) {
        # j goes through all groups starting from i:
        for  (j in (i + 1):nlevels(conditions)) {
            i_label <- levels(conditions)[i]
            j_label <- levels(conditions)[j]


            # Old version, 1/2023:
            # pairwise_results <- as.data.frame(results(test_results, contrast=c("condition",j_label,i_label))) # note: j,i => j-i
            res2 <- results(test_results, contrast = c("condition", j_label, i_label)) # note: j,i => j-i

            # LFC shrink, needs to be after results function. (This step was previously inside DESeq function)
            resLFC <- lfcShrink(dds, res = res2, type = "ashr")
            res2 <- resLFC

            # Summary into a txt file
            sink("summary.txt", append = TRUE)
            print(paste("Comparison: ", i_label, "_vs_", j_label, sep = ""))
            summary(res2, alpha = p.value.cutoff)
            sink()

            pairwise_results <- as.data.frame(res2)

            # building the table:
            if (is.null(res)) res <- pairwise_results else res <- cbind(res, pairwise_results)

            results_name <- c(results_name, paste(i_label, "_vs_", j_label, sep = ""))

            # MA plot after LFC shrink:
            plotMA(resLFC, alpha = p.value.cutoff, main = c("DESeq2 MA-plot with LFC shrink, FDR =", p.value.cutoff, " comparison = ", results_name[length(results_name)]), ylim = c(-2, 2)) # , cex = 3
        }
    }

    # "stat" column no longer there, 1/2023. Old version: colnames(res) <- paste(colnames(res), rep(results_name,each=6), sep=".")
    number_of_cols_in_res <- length(pairwise_results)
    colnames(res) <- paste(colnames(res), rep(results_name, each = number_of_cols_in_res), sep = ".")
    min_padj <- apply(res[, grep("padj", colnames(res))], 1, min)
    sig <- cbind(dat, res, min_padj = min_padj)
    sig <- sig[(sig$min_padj <= p.value.cutoff), ]
    sig <- sig[!(is.na(sig$min_padj)), ]
    sig <- sig[order(sig$min_padj), ]
    sig <- sig[, -grep("min_padj", colnames(sig))]


    # Output significant DEGs
    if (dim(sig)[1] > 0) {
        ndat <- ncol(dat)
        npairdat <- ncol(pairwise_results) # number of columns in one comparison
        ncomp <- choose(nlevels(conditions), 2) # number of comparisons
        rounded_sig <- sig
        for (i in 1:ncomp) {
            skip <- ndat + (i - 1) * npairdat # how many columns to skip
            round2 <- skip + 1:(npairdat - 2) # everything except the last two columns (p-value & padj), beginning after skipped cols
            round4 <- skip + (npairdat - 1):npairdat # the last two columns (p-value & padj), beginning after skipped cols
            rounded_sig[, round2] <- round(rounded_sig[, round2], digits = round_others) # round to 2 digits everything except p-value & padj
            rounded_sig[, round4] <- format(rounded_sig[, round4], digits = round_pvals, scientific = T) # round to 4 digits p-value & padj
        }
        write.table(rounded_sig, file = "de-list-deseq2.tsv", sep = "\t", row.names = T, col.names = T, quote = F)
    }
}


# Create a template output table for plotting.
# If having N genes and 3 comparisons, this conversion results in a data matrix that has Nx3 rows (cells?)
output_table <- NULL
colnames(res) <- gsub("\\..*$", "", colnames(res))
for (i in grep("baseMean$", colnames(res))) {
    # How many columns per comparison = col_size
    col_size <- grep("padj", colnames(res))[1] - grep("baseMean", colnames(res))[1]
    # Clue together the counts (=dat) and the i:th comparison from the results table (=res)
    output_table <- rbind(output_table, cbind(dat, res[, (i:(i + col_size))]))
}
rownames(output_table) <- make.names(rep(rownames(res), length(grep("baseMean$", colnames(res)))), unique = T)
output_table <- as.data.frame(output_table)

# testing:
#  write.table(output_table, file="output_table.tsv", sep="\t", row.names=T, col.names=T, quote=F)
#  save(output_table, file = "res.Robj")


# If genomic coordinates are present, output a sorted BED file for region matching tools
if (bed == "yes") {
    if ("chr" %in% colnames(dat)) {
        if (dim(sig)[1] > 0) {
            bed <- output_table[, c("chr", "start", "end")]
            bed <- as.data.frame(bed)
            if (is.null(results_name)) {
                gene_names <- rownames(res)
            } else {
                gene_names <- paste(rep(results_name, each = nrow(res)), rownames(res), sep = "")
            }
            bed <- cbind(bed, name = gene_names) # name
            bed <- cbind(bed, score = output_table[, "log2FoldChange"]) # score
            bed <- bed[(output_table$padj <= p.value.cutoff & (!(is.na(output_table$padj)))), ]
            bed <- sort.bed(bed)
            write.table(bed, file = "de-list-deseq2.bed", sep = "\t", row.names = F, col.names = F, quote = F)
        }
    }
}


# Combined MA-plot, when there are more than 2 groups.
if (length(unique(groups)) > 2) {
    # Define function for making MA-plot.
    plotDE <- function(res) {
        plot(res$baseMean, res$log2FoldChange,
            log = "x", ylim = c(-2, 2), pch = 20, cex = .25, col = ifelse(res$padj < p.value.cutoff, "blue", "black"),
            main = "MA plot with LFC shrink, all comparisons", xlab = "mean counts", ylab = "log2(fold change)"
        )
    }
    # Make MA-plot
    plotDE(unique(output_table))
    legend(x = "topleft", legend = c("significant based on adj. p-val", "not significant"), col = c("blue", "black"), cex = 1, pch = 19)
    abline(h = 0, col = "darkgreen", lwd = 2)
    # abline(h = c(-1, 0, 1), col = c("dodgerblue", "darkgreen", "dodgerblue"), lwd = 2) #
}

# Make dispersion plot
plotDispEsts(dds, main = "Dispersion plot", cex = 0.2)
legend(x = "topright", legend = "fitted dispersion", col = "red", cex = 1, pch = "-")


## Make histogram of p-values with overlaid significance cutoff. When more than two groups, min.pvalue is taken over all comparisons for genes
hist(output_table$pval, breaks = 100, col = "blue", border = "slateblue", freq = FALSE, main = "P-value distribution", xlab = "p-value", ylab = "proportion (%)")
hist(output_table$padj, breaks = 100, col = "red", border = "slateblue", add = TRUE, freq = FALSE)
abline(v = p.value.cutoff, lwd = 2, lty = 2, col = "black")
legend(x = "topright", legend = c("p-values", "adjusted p-values", "significance cutoff"), col = c("blue", "red", "black"), cex = 1, pch = 15)


# Close pdf
dev.off()

# EOF
