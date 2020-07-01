# TOOL mothur-sffinfo.R: "Extract sequences from a SFF file" (Extract a FASTA sequence file and a QUAL quality file from a SFF file. Sequences can optionally be trimmed for quality. This tool is based on the Mothur tool sffinfo.)
# INPUT reads.sff: "SFF file" TYPE GENERIC
# OUTPUT OPTIONAL reads.fasta
# OUTPUT OPTIONAL reads.qual
# OUTPUT OPTIONAL reads.raw.fasta
# OUTPUT OPTIONAL reads.raw.qual
# PARAMETER OPTIONAL trim: "Trim reads for quality" TYPE [yes, no] DEFAULT no (Trim sequences and quality scores to the clipQualLeft and clipQualRight values.)

# AMS 19.06.2013

source(file.path(chipster.common.path,"tool-utils.R"))

# binary
binary <- c(file.path(chipster.tools.path,"mothur","mothur"))
version <- system(paste(binary,"--version"),intern = TRUE)
documentVersion("Mothur",version)

# Options
sffinfo.options <- ""
sffinfo.options <- paste(sffinfo.options,"sffinfo(sff=reads.sff")
if (trim == "no") {
  sffinfo.options <- paste(sffinfo.options,", trim=F",sep = ",")
}
sffinfo.options <- paste(sffinfo.options,")",sep = "")

# Write batch file
documentCommand(sffinfo.options)
write(sffinfo.options,"sffinfo.mth",append = FALSE)

# command
command <- paste(binary,"sffinfo.mth")

# run
system(command)

