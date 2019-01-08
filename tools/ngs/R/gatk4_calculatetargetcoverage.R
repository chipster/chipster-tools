# TOOL gatk4_calculatetargetcoverage.R: "Collect proportional coverage using target intervals and read data" (GATK4)
# INPUT reads.bam: "Reads BAM file" TYPE BAM
# INPUT targets.tsv: "Padded targets TSV" TYPE GENERIC
# OUTPUT toss_coverage.tsv
# PARAMETER gatk4.transform: "Transformation" TYPE [RAW, PCOV] DEFAULT RAW (Transformation to perform to the individual count output values.)
# PARAMETER gatk4.groupby: "Group by" TYPE [COHORT, SAMPLE, READGROUP] DEFAULT COHORT (Group counts by Cohort \(all samples in one\), Sample or ReadGroup.)
# PARAMETER gatk4.targetinfo: "Target info" TYPE [COORDS, NAME, FULL] DEFAULT COORDS (Target identification information to be showed in outputs.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))
unzipIfGZipFile("targets.tsv")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk-launch"))
samtools.binary <- c(file.path(chipster.tools.path, "samtools", "samtools"))

# Index BAM
system(paste(samtools.binary, "index reads.bam > reads.bam.bai"))

command <- paste(gatk.binary, "CalculateTargetCoverage", "-I reads.bam", "-T targets.tsv", "--transform", gatk4.transform, "--groupBy", gatk4.groupby, "--targetInformationColumns", gatk4.targetinfo, "--disableReadFilter NotDuplicateReadFilter", "-O toss_coverage.t.tsv", "")

runExternal(command)

# Cleand headers to make compatible with Chipster tsv viewer
system("grep -v '#' toss_coverage.t.tsv > toss_coverage.tsv")