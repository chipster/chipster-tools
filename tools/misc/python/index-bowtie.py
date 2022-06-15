# TOOL index-bowtie.py: "Create Bowtie index" ()
# INPUT input.fa TYPE GENERIC
# OUTPUT output{...}
# RUNTIME python3

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils
import subprocess
from typing import Iterable
import os

def main():

    input_fa = "input.fa"

    bowtie_build = chipster_tools_path + "/bowtie/bowtie-build"
    bowtie_inspect = chipster_tools_path + "/bowtie/bowtie-inspect"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, '.fa')

    run_process([bowtie_build, input_fa, fasta_basename])

    print("inspect index")
    inspect_output = "inspect_output.txt"
    run_bash(bowtie_inspect + " -n " + fasta_basename + " > " + inspect_output)

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
        raise RuntimeError("bowtie indexing of genome " + fasta_basename + " failed. Chromosomes in fasta: " + fasta_chr + ", chromosomes in index: " + index_chr)
    
    run_process(["ls", "-lah"])

    # rename index files to fixed output names and define real names for the client

    output_names = {}

    for file in os.listdir("."):
        # .ebwt or .ebwtl
        if ".ebwt" in file:
            output_name = file.replace(fasta_basename, "output") + str(len(output_names))
            os.rename(file, output_name)
            output_names[output_name] = file

    tool_utils.write_output_definitions(output_names)
    
def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

main()
