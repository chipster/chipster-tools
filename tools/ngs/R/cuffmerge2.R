# TOOL cuffmerge2.R: "Merge transcript assemblies with Cuffmerge" (Given several transcript GTF files obtained by Cufflinks, Cuffmerge merges them into one. The merged GTF file can be used in differential expression analysis with Cuffdiff.)
# INPUT annotation{...}.gtf: "GTF files" TYPE GTF
# OUTPUT OPTIONAL merged.gtf
# OUTPUT OPTIONAL cuffmerge.log
# PARAMETER OPTIONAL exclude.duplicates: "Exclude transcripts with duplicate ids" TYPE [yes, no] DEFAULT no (Exclude all entries with duplicated transcript ids. You should only do this if you get an error message about duplicated ids.)

# AMS 21.11.2012
# AMS 11.01.2013 Removed unnecessary outputs
# EK 21.1.2013
# AMS 11.11.2013 Added thread support

# binary
cuffmerge.binary <- c(file.path(chipster.tools.path, "cufflinks2", "cuffmerge"))

# Set PATH so cuffmerge can find gtf_to_sam
cufflinkspath <- c(file.path(chipster.tools.path, "cufflinks2"))
setpathcommand <- paste("export PATH=$PATH:", cufflinkspath, ";", sep="")

# Remove duplicates
if (exclude.duplicates == "yes"){
	# awk command to list duplicates
	listduplicates.awk <- paste("| awk -F \'transcript_id\' \'{print $2}\' |awk -F\'\\\"\' \'{print $2}\' | sort | uniq -c | grep -v -w \"1\" |awk \'{print $2}\'")
	system(paste("echo \"Excluded duplicate ids:\" >> exclude.tmp"))
	# Loop through the inputs
	input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
	for (i in 1:nrow(input.names)) {
		# Generate names for the list files and new gtf files
		list.name <- paste(input.names[i,1], ".list", sep="")
		gtf.name <- paste(input.names[i,1], ".nd.gtf", sep="")
		# Make a list of duplicated transcript ids
		system(paste("grep -w transcript", input.names[i,1], listduplicates.awk, "> ", list.name))
		# Make a gtf without the duplicated entries
		system(paste("grep -v -w -f", list.name, input.names[i,1], ">", gtf.name))
		# Log file
		system(paste("printf \"", input.names[i,2], ":\t\`wc -l < ", list.name, "\`\n\" >> exclude.tmp", sep=""))
	}
	system("mv exclude.tmp cuffmerge.log")
	# Use the new gtf files
	system("ls *.nd.gtf > assemblies.txt")
} else{
	# Use the original gtf files
	system("ls *.gtf > assemblies.txt")
}

# Cuffmerge command
command <- paste(setpathcommand, cuffmerge.binary, "-p", chipster.threads.max, "assemblies.txt", "2>> log.tmp")

# run
#stop(paste('CHIPSTER-NOTE: ', command))
system(command)

# Rename files
if (file.exists("merged_asm/merged.gtf") && file.info("merged_asm/merged.gtf")$size > 0) {
	system("mv merged_asm/merged.gtf merged.gtf")
}else{
	system("cat log.tmp >> cuffmerge.log")
}
