# TOOL sort_gtf.R: "Sort GTF" (Sort GTF file)
# INPUT input.gtf: "GTF file" TYPE GTF
# OUTPUT sorted.gtf  

# source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "gtf-utils.R"))

sort.gtf("input.gtf", "sorted.gtf")
