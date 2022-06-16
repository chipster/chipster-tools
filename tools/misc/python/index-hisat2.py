# TOOL index-hisat2.py: "Create Hisat2 index" ()
# INPUT input.fa TYPE GENERIC
# INPUT input.gtf TYPE GENERIC
# OUTPUT output{...}
# RUNTIME python3
# SLOTS 5

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils
import subprocess
from typing import Iterable
import os

def main():

    input_fa = "input.fa"
    input_gtf = "input.gtf"

    extract_splice_sites = chipster_tools_path + "/hisat2/extract_splice_sites.py"
    extract_exons = chipster_tools_path + "/hisat2/extract_exons.py"
    hisat2_build = chipster_tools_path + "/hisat2/hisat2-build"
    hisat2_inspect = chipster_tools_path + "/hisat2/hisat2-inspect"

    session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
    gtf_basename = tool_utils.remove_postfix(session_input_gtf, '.gtf')

    # Create a splice site file
    run_bash(extract_splice_sites + " " + input_gtf  + " > " + gtf_basename + ".ss")

    # Create an exon file
    run_bash(extract_exons + " " + input_gtf + " > " + gtf_basename + ".exon")

    # Run hisat2-build with splice site and exon annotations, this consumes lot of ram: ~200G for human
    # rat would require 0.5 TB RAM with --ss and --exon https://www.biostars.org/p/344449/
    if gtf_basename.startswith("Rattus_norvegicus"):
        print("not using splice sites and exon for rat")
        run_process([hisat2_build, "-p", chipster_threads_max, input_fa, "output"])
    else:
        run_process([hisat2_build, "-p", chipster_threads_max, "--ss", gtf_basename + ".ss", "--exon", gtf_basename + ".exon", input_fa, "output"])

    print("inspect index")
    inspect_output = "inspect_output.txt"
    run_bash(hisat2_inspect + " -n output > " + inspect_output)

    fasta_chr = 0
    index_chr = 0

    print("calculate chromosomes in fasta")
    with open(input_fa) as file:
        for line in file:
            if line.startswith(">"):
                fasta_chr += 1

    print("calculate chromosomes in index")
    with open(inspect_output) as file:
        for line in file:
            index_chr += 1

    if fasta_chr != index_chr:
        raise RuntimeError("hisat2 indexing of genome " + gtf_basename + " failed. Chromosomes in fasta: " + fasta_chr + ", chromosomes in index: " + index_chr)
    
    run_process(["ls", "-lah"])

    # rename index files to fixed output names and define real names for the client

    output_names = {}

    for file in os.listdir("."):
        # .ebwt or .ebwtl
        if file.startswith("output"):
            file_extension = file.replace("output", "")
            output_name = "output" + str(len(output_names))
            os.rename(file, output_name)
            output_names[output_name] = gtf_basename + file_extension

    tool_utils.write_output_definitions(output_names)
    
def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

main()
