# TOOL sort_gtf.R: "Sort GTF" (Sort GTF file)
# INPUT input.gtf: "GTF file" TYPE GTF
# OUTPUT sorted.gtf

# source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "gtf-utils.R"))

sort.gtf("input.gtf", "sorted.gtf")
