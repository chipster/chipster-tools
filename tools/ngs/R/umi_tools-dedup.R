# TOOL umi_tools-dedup.R: "UMI-tools deduplicate" (UMI-tools dedup)
# INPUT input.bam: "BAM file" TYPE GENERIC
# OUTPUT deduplicated.bam
# OUTPUT stats_edit_distance.tsv
# OUTPUT stats_per_umi.tsv
# OUTPUT stats_per_umi_per_position.tsv
# PARAMETER grouping.method: "Grouping method" TYPE [unique, directional] DEFAULT unique (What method to use to identify group of reads with the same or similar UMIs. Please see the manual page for details.)
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.1

source(file.path(chipster.common.path, "tool-utils.R"))

# Index BAM
samtools.binary <- c(file.path(chipster.tools.path, "samtools-1.2", "samtools"))
system(paste(samtools.binary, "index input.bam > input.bam.bai"))

# active Python virtual environment "venv"
venv_root <- "/opt/chipster/tools/umi-tools/venv"
venv_path <- paste(Sys.getenv("PATH"), paste(venv_root, "bin", sep="/"), sep = ":")
Sys.setenv(PATH = venv_path, VIRTUAL_ENV = venv_root)
 
version <- system("umi_tools --version | cut -d : -f 2-",intern = TRUE)
documentVersion("UMI-tools",version)
umi.command <- paste("umi_tools dedup -I input.bam --output-stats=stats -S deduplicated.bam --method",grouping.method)
documentCommand(umi.command)
runExternal(umi.command)

# Output names
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$input.bam)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=4, ncol=2)
outputnames[1,] <- c("deduplicated.bam", paste(basename, "_dedup.bam", sep=""))
outputnames[2,] <- c("stats_edit_distance.tsv", paste(basename, "_edit_distance.tsv", sep=""))
outputnames[3,] <- c("stats_per_umi.tsv", paste(basename, "_per_umi.tsv", sep=""))
outputnames[4,] <- c("# OUTPUT stats_per_umi_per_position.tsv", paste(basename, "_per_umi_per_position.tsv", sep=""))
# Write output definitions file
write_output_definitions(outputnames)