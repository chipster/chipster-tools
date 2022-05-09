# TOOL picard-addreadgroup.R: "Add Read Group to BAM" (Adds Read Group line \(\@RG\) to BAM header. This line contains additional information that is required by some applications like Mutect2.)
# INPUT input.bam : "Reads" TYPE BAM
# OUTPUT output.bam
# PARAMETER rgid: "Read group identifier" TYPE STRING DEFAULT "1" (Read group identifier.)
# PARAMETER rgsm: "Sample name for read group" TYPE STRING DEFAULT "Sample" (The name of the sample sequenced in this read group.)
# PARAMETER rglb: "Library identifier for read group" TYPE STRING DEFAULT "1" (DNA preparation library identifier. The Mark Duplicates tool uses this field to determine which read groups might contain molecular duplicates, in case the same DNA library was sequenced on multiple lanes.)
# PARAMETER rgpl: "Platform for read group" TYPE [ none: "Not defined",CAPILLARY, DNBSEQ, HELICOS, ILLUMINA, IONTORRENT, LS454, ONT, PACBIO, SOLID] DEFAULT ILLUMINA (Platform\/technology used to produce the read.)
# PARAMETER rgpu: "Read Group platform unit" TYPE STRING DEFAULT "1" (Read Group platform unit \(e.g., flowcell-barcode.lane for Illumina or slide for SOLiD\).)

source(file.path(chipster.common.path, "tool-utils.R"))

# Picard binary
picard.binary <- file.path(chipster.tools.path, "picard-tools", "picard.jar")

# Picard basic options
picard.options <- paste("I=input.bam O=output.bam")

#Read group definitions
picard.options <- paste(picard.options," ","RGID=",rgid,sep="")
picard.options <- paste(picard.options," ","RGSM=",rgsm,sep="")
picard.options <- paste(picard.options," ","RGLB=",rglb,sep="")
picard.options <- paste(picard.options," ","RGPL=",rgpl,sep="")
picard.options <- paste(picard.options," ","RGPU=",rgpu,sep="")

# Picard command
picard.command <- paste("java -jar", picard.binary, "AddOrReplaceReadGroups")

# Picard Version
system(paste(picard.command,"--version 2> version.tmp"))
version <- system("cat version.tmp",intern = TRUE)
documentVersion("Picard", version)

# Run command
picard.command <- paste(picard.command, picard.options)
documentCommand(picard.command)
runExternal(picard.command)

# Handle output names
inputnames <- read_input_definitions()
name <- inputnames$input.bam

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("output.bam", paste(name, sep=""))

# Write output definitions file
write_output_definitions(outputnames)
