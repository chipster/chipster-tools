# TOOL umi-tools.R: "Umi Tools example" (Umi Tools example)
# INPUT input.tsv: "TSV file" TYPE GENERIC
# OUTPUT output.tsv
# RUNTIME R-4.1.1

source(file.path(chipster.common.path, "tool-utils.R"))

# active Python virtual environment "venv"
venv_root <- "/opt/chipster/tools/umi-tools/venv"
venv_path <- paste(Sys.getenv("PATH"), paste(venv_root, "bin", sep="/"), sep = ":")
Sys.setenv(PATH = venv_path, VIRTUAL_ENV = venv_root)

runExternal("umi_tools --help")

runExternal(paste("mv", "input.tsv", "output.tsv"))
