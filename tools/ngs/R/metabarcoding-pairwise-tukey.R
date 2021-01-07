# TOOL metabarcoding-pairwise-tukey.R: "Post-hoc pairwise Tukey's HSD for OTU abundance data" (Performs post-hoc pairwise comparisons of group dispersions using a Tukey's Honestly Significant Difference \(HSD\) test. Should only be used after a significant PERMDISP result. Performs up to three tests \(for three different phenodata variables\), depending on preceding steps taken when using the PERMDISP tool. Requires an Rda file \(ps_disp.Rda\) produced by the PERMDISP tool as the input.)
# INPUT ps.Rda: "Data set in Rda format" TYPE GENERIC
# OUTPUT pairwise_tukey_table.txt: pairwise_tukey_table.txt
# OUTPUT pairwise_tukey_plot.pdf: pairwise_tukey_plot.pdf
# PARAMETER howmany: "No. of PERMDISP analyses in preceding step (PERMDISP tool)?" TYPE INTEGER FROM 1 TO 3 (Number of PERMDISP analyses performed during the preceding analysis step, i.e. when using the PERMDISP tool.)
# RUNTIME R-3.6.1-phyloseq

# JH 2020

# Load packages
library(vegan)

# Load data
load("ps.Rda")

# Tukey HSD
if (howmany == "1" &&
	permutest_beta1$tab[1,6] < 0.05	){
		tukey_beta1 <- TukeyHSD(beta1)
}
if (howmany == "1" &&
	permutest_beta1$tab[1,6] > 0.05	){
		stop("No significant PERMDISP result (analysis aborted)")
}
if (howmany == "2" &&
	permutest_beta1$tab[1,6] < 0.05 &&
	permutest_beta2$tab[1,6] < 0.05){
		tukey_beta1 <- TukeyHSD(beta1)
		tukey_beta2 <- TukeyHSD(beta2)
}
if (howmany == "2" &&
	permutest_beta1$tab[1,6] > 0.05 &&
	permutest_beta2$tab[1,6] < 0.05){
		tukey_beta1 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta2 <- TukeyHSD(beta2)
}
if (howmany == "2" &&
	permutest_beta1$tab[1,6] < 0.05 &&
	permutest_beta2$tab[1,6] > 0.05){
		tukey_beta1 <- TukeyHSD(beta1)
		tukey_beta2 <- "PERMDISP p > 0.05; post-hoc test not performed"
}
if (howmany == "2" &&
	permutest_beta1$tab[1,6] > 0.05 &&
	permutest_beta2$tab[1,6] > 0.05){
		stop("No significant PERMDISP results (analysis aborted)")
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] < 0.05 &&
	permutest_beta2$tab[1,6] < 0.05 &&
	permutest_beta3$tab[1,6] < 0.05){
		tukey_beta1 <- TukeyHSD(beta1)
		tukey_beta2 <- TukeyHSD(beta2)
		tukey_beta3 <- TukeyHSD(beta3)
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] > 0.05 &&
	permutest_beta2$tab[1,6] < 0.05 &&
	permutest_beta3$tab[1,6] < 0.05){
		tukey_beta1 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta2 <- TukeyHSD(beta2)
		tukey_beta3 <- TukeyHSD(beta3)
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] < 0.05 &&
	permutest_beta2$tab[1,6] > 0.05 &&
	permutest_beta3$tab[1,6] < 0.05){
		tukey_beta1 <- TukeyHSD(beta1)
		tukey_beta2 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta3 <- TukeyHSD(beta3)
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] < 0.05 &&
	permutest_beta2$tab[1,6] < 0.05 &&
	permutest_beta3$tab[1,6] > 0.05){
		tukey_beta1 <- TukeyHSD(beta1)
		tukey_beta2 <- TukeyHSD(beta2)
		tukey_beta3 <- "PERMDISP p > 0.05; post-hoc test not performed"
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] > 0.05 &&
	permutest_beta2$tab[1,6] > 0.05 &&
	permutest_beta3$tab[1,6] < 0.05){
		tukey_beta1 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta2 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta3 <- TukeyHSD(beta3)
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] > 0.05 &&
	permutest_beta2$tab[1,6] < 0.05 &&
	permutest_beta3$tab[1,6] > 0.05){
		tukey_beta1 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta2 <- TukeyHSD(beta2)
		tukey_beta3 <- "PERMDISP p > 0.05; post-hoc test not performed"
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] < 0.05 &&
	permutest_beta2$tab[1,6] > 0.05 &&
	permutest_beta3$tab[1,6] > 0.05){
		tukey_beta1 <- TukeyHSD(beta1)
		tukey_beta2 <- "PERMDISP p > 0.05; post-hoc test not performed"
		tukey_beta3 <- "PERMDISP p > 0.05; post-hoc test not performed"
}
if (howmany == "3" &&
	permutest_beta1$tab[1,6] > 0.05 &&
	permutest_beta2$tab[1,6] > 0.05 &&
	permutest_beta3$tab[1,6] > 0.05){
		stop("No significant PERMDISP results (analysis aborted)")
}

# Plots showing differences in mean levels of group (95% family-wise confidence level)

# Open a report PDF
pdf("pairwise_tukey_plot.pdf")

if (howmany == "1" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
}
if (howmany == "2" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
		plot(tukey_beta2, las = 1)
}
if (howmany == "2" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 == "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
}
if (howmany == "2" &&
	tukey_beta1 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta2, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
		plot(tukey_beta2, las = 1)
		plot(tukey_beta3, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta2, las = 1)
		plot(tukey_beta3, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
		plot(tukey_beta3, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 == "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
		plot(tukey_beta2, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 == "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta1, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 != "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 == "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta2, las = 1)
}
if (howmany == "3" &&
	tukey_beta1 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta2 == "PERMDISP p > 0.05; post-hoc test not performed" &&
	tukey_beta3 != "PERMDISP p > 0.05; post-hoc test not performed"){
		plot(tukey_beta3, las = 1)
}

# Close the report PDF
dev.off()

# Print results tables
if (howmany == "1"){
sink("pairwise_tukey_table.txt")
	cat("\n\n\n")
	cat("### Pairwise Tukey HSD comparisons ###\n")
	cat("\n\n\n")
	print(tukey_beta1)
	cat("\n\n\n")
sink()
}
if (howmany == "2"){
sink("pairwise_tukey_table.txt")
	cat("\n\n\n")
	cat("### Pairwise Tukey HSD comparisons ###\n")
	cat("\n\n\n")
	cat("### Test 1 ###\n")
	cat("\n\n\n")
	print(tukey_beta1)
	cat("\n\n\n")
	cat("### Test 2 ###\n")
	cat("\n\n\n")
	print(tukey_beta2)
	cat("\n\n\n")
sink()
}
if (howmany == "3"){
sink("pairwise_tukey_table.txt")
	cat("\n\n\n")
	cat("### Pairwise Tukey HSD comparisons ###\n")
	cat("\n\n\n")
	cat("### Test 1 ###\n")
	cat("\n\n\n")
	print(tukey_beta1)
	cat("\n\n\n")
	cat("### Test 2 ###\n")
	cat("\n\n\n")
	print(tukey_beta2)
	cat("\n\n\n")
	cat("### Test 3 ###\n")
	cat("\n\n\n")
	print(tukey_beta3)
	cat("\n\n\n")
sink()
}
