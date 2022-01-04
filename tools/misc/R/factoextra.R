# TOOL factoextra.R: "factoextra example" (factoextra example)
# INPUT input.tsv: "TSV file" TYPE GENERIC
# OUTPUT output.tsv
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.1-statistics

source(file.path(chipster.common.path, "tool-utils.R"))

library("factoextra")

library("ggrepel")

runExternal(paste("mv", "input.tsv", "output.tsv"))

library("ggplot2")