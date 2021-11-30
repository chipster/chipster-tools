# TOOL umi-tools.R: "Umi Tools example" (Umi Tools example)
# INPUT input.tsv: "GTF file" TYPE GENERIC
# OUTPUT output.tsv
# IMAGE comp-20.04-r-deps
# RUNTIME R-4.1.1

source(file.path(chipster.common.path, "tool-utils.R"))

runExternal("echo $0")

runExternal("bash -c 'source /opt/chipster/tools/umi-tools/venv/bin/activate; umi_tools --help'")

runExternal(paste("mv", "input.tsv", "output.tsv"))
