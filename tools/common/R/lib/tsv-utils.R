# First sort by chromosome, then by start coordinate.
# Chromosomes that are numeric are compared numerically. Chromosomes that
# are non-numeric are compared lexically (in their normalised form).
# Numeric names are always considered smaller than non-numeric.
#
sort.tsv <- function(input, output, chr) {
    system(paste("java -cp  '", chipster.java.libs.path, "/*' fi.csc.chipster.tools.ngs.SortTsv ", input, " ", output, " ", chr, sep = ""))
}
