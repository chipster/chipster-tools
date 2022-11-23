# TOOL copy-fasta.py: "Copy genome fasta" ()
# OUTPUT output.fa
# PARAMETER organism: "Genome" TYPE ["FILES genomes/fasta .fa"] DEFAULT "SYMLINK_TARGET genomes/fasta/default .fa" (Genome or transcriptome that you would like copy.)
# RUNTIME python3

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils

import os

output_fa = "output.fa"

fasta_path = chipster_tools_path + "/genomes/fasta/" + organism + ".fa"

os.symlink(fasta_path, output_fa)

output_names = {
	output_fa: organism + ".fa",
}

tool_utils.write_output_definitions(output_names)
