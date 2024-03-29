# TOOL gatk4-filterbyorientation.R: "GATK4 -Filter Mutect2 calls by orientation bias" (Additionally filter Mutect2 somatic variant calls for sequence context-dependent artifacts, e.g. OxoG or FFPE deamination. Tool is based on GATK4 FilterByOrientationBias.)
# INPUT input.vcf: "Mutect calls" TYPE GENERIC
# INPUT metrics_file: "Detail metrics file" TYPE GENERIC
# OUTPUT OPTIONAL filtered.vcf
# OUTPUT OPTIONAL filtered.vcf.summary
# OUTPUT OPTIONAL gatk_log.txt
# PARAMETER OPTIONAL gatk.artifactmodes: "Artifacts on the forward strand" TYPE STRING (PreAdapter Detail artifacts of interest on the forward strand. Format CtoA for a single artifact. If you want to give multiple artifacts, use comma to separate them:  CtoA,TtoG. Artifacts must be one base to one base \(e.g. 'CCtoCA' is illegal\). G>T is OxoG.)

source(file.path(chipster.common.lib.path, "gatk-utils.R"))
source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.lib.path, "zip-utils.R"))

unzipIfGZipFile("metrics_file")

formatGatkVcf("input.vcf")

# binaries
gatk.binary <- c(file.path(chipster.tools.path, "GATK4", "gatk"))

# Run GetPileupSummaries
command <- paste(gatk.binary, "FilterByOrientationBias", "-O filtered.vcf", "-V input.vcf.gz", "-P metrics_file")

if (nchar(gatk.artifactmodes) > 0) {
    # Split into individual entries
    artifacts <- strsplit(gatk.artifactmodes, ",")[[1]]
    artifacts <- gsub("to", "/", artifacts)
    for (i in 1:length(artifacts)) {
        artifactmode <- paste("'", trimws(artifacts[i]), "'", sep = "", collapse = "")
        command <- paste(command, "--artifact-modes", artifactmode)
    }
}

runExternal(command)

# Return error message if no result
if (fileNotOk("filtered.vcf")) {
    system("mv stderr.log gatk_log.txt")
}

# read input names
inputnames <- read_input_definitions()

basename <- strip_name(inputnames$input.vcf)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("filtered.vcf", paste(basename, "_filtered_unbiased.vcf", sep = ""))
outputnames[2, ] <- c("filtered.vcf.summary", paste(basename, "_filtered_unbiased_summary.tsv", sep = ""))

# Write output definitions file
write_output_definitions(outputnames)
