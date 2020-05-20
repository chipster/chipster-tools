# TOOL gtftobed.R: "Convert GTF to BED" (Converts GTF to BED format.)
# INPUT file.gtf: "GTF file" TYPE GENERIC
# OUTPUT OPTIONAL gtftobed.bed

source(file.path(chipster.common.path,"tool-utils.R"))
source(file.path(chipster.common.path,"zip-utils.R"))

# binary
binary <- c(file.path(chipster.tools.path,"gtf2bed","gtf2bed.pl"))

# check out if the file is compressed and if so unzip it
unzipIfGZipFile("file.gtf")

# command
command <- paste(binary,"file.gtf > gtftobed.bed")
documentCommand(command)

# run
system(command)

# read input names
inputnames <- read_input_definitions()
base <- remove_extension(inputnames$file.gtf)

# Make a matrix of output names
outputnames <- matrix(NA,nrow = 2,ncol = 2)
outputnames[1,] <- c("gtftobed.bed",paste(base,".bed",sep = ""))

# Write output definitions file
write_output_definitions(outputnames)