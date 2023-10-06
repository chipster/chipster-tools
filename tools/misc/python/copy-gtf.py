# TOOL copy-gtf.py: "Copy genome GTF" ()
# OUTPUT output.gtf
# PARAMETER organism: "Genome" TYPE ["FILES genomes/gtf .gtf"] DEFAULT "SYMLINK_TARGET genomes/gtf/default .gtf" (Genome or transcriptome that you would like copy.)
# RUNTIME python3

import tool_utils

import os

output_gtf = "output.gtf"

gtf_path = chipster_tools_path + "/genomes/gtf/" + organism + ".gtf"

os.symlink(gtf_path, output_gtf)

output_names = {
    output_gtf: organism + ".gtf",
}

tool_utils.write_output_definitions(output_names)
