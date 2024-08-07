# TOOL umi_tools-extract.R: "Extract UMIs from QuantSeq reads" (Given a Lexogen QuantSeq FASTQ file, this tool removes the TATA spacer and extracts the 6-base or 12-base UMI and adds it to the read name. This tool is based on the UMI-tools extract.)
# INPUT reads.fq.gz: "Reads" TYPE GENERIC
# OUTPUT output.fq.gz
# OUTPUT OPTIONAL log.txt
# PARAMETER umi_length: "UMI length" TYPE [6, 12] DEFAULT 6 (Length of UMI to extract.)
# PARAMETER OPTIONAL log: "Produce log file" TYPE [yes,no] DEFAULT no (Produce log file.)
# RUNTIME R-4.1.1

source(file.path(chipster.common.lib.path, "tool-utils.R"))

# active Python virtual environment "venv"
venv_root <- "/opt/chipster/tools/umi-tools/venv"
venv_path <- paste(Sys.getenv("PATH"), paste(venv_root, "bin", sep = "/"), sep = ":")
Sys.setenv(PATH = venv_path, VIRTUAL_ENV = venv_root)

version <- system("umi_tools --version | cut -d : -f 2-", intern = TRUE)
documentVersion("UMI-tools", version)
if (umi_length == "6"){
  umi.command <- paste("umi_tools extract --extract-method=regex --bc-pattern \"(?P<umi_1>.{6})(?P<discard_1>TATA).*\" -L log.t.txt -I reads.fq.gz -S output.fq.gz")
else{
  umi.command <- paste("umi_tools extract --extract-method=regex --bc-pattern \"(?P<umi_1>.{12})(?P<discard_1>TATA).*\" -L log.t.txt -I reads.fq.gz -S output.fq.gz")
}
documentCommand(umi.command)
runExternal(umi.command)

# Clean log file
if (log == "yes") {
    system("grep INFO log.t.txt |grep -v Parsed > log.txt")
}

# Output names
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$reads.fq.gz)
print(basename)
# Make a matrix of output names
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("output.fq.gz", paste(basename, "_extracted.fq.gz", sep = ""))
outputnames[2, ] <- c("log.txt", paste(basename, "_extracted.log", sep = ""))
# Write output definitions file
write_output_definitions(outputnames)
