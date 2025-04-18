# TOOL index-bowtie2-tar.py: "Create Bowtie2 index" ()
# INPUT input.fa TYPE FASTA
# OUTPUT output.tar
# RUNTIME python3
# SLOTS 2

import tool_utils
import subprocess
from typing import Iterable
import os


def main():
    input_fa = "input.fa"
    bowtie2_path = chipster_tools_path + "/bowtie2"
    bowtie2_build = chipster_tools_path + "/bowtie2/bowtie2-build"
    bowtie2_inspect = chipster_tools_path + "/bowtie2/bowtie2-inspect"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, ".fa")

    run_process(
        [bowtie2_build, "--threads", chipster_threads_max, input_fa, fasta_basename]
    )

    print("inspect index")
    inspect_output = "inspect_output.txt"
    run_bash(bowtie2_inspect + " -n " + fasta_basename + " > " + inspect_output)

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
        raise RuntimeError(
            "bowtie2 indexing of genome "
            + fasta_basename
            + " failed. Chromosomes in fasta: "
            + str(fasta_chr)
            + ", chromosomes in index: "
            + str(index_chr)
        )

    run_process(["ls", "-lah"])

    # create tar package

    index_files = []

    # bowtie2 index files
    for file in os.listdir("."):
        if file.endswith(".bt2"):
            index_files.append(file)
            
    run_process(["tar", "-cf", "output.tar"] + index_files)

    tool_utils.write_output_definitions({
        "output.tar": fasta_basename + ".tar"
    })
    
    # save version information
    version = subprocess.check_output([bowtie2_build, "--version"])
    version_number = str(version).split("\\n")[0].split(" ")[2]
    
    version_utils.document_version("Bowtie2", version_number)


def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])


def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError(
            "process failed with return code: "
            + str(process.returncode)
            + ", command: "
            + str(cmd)
        )


main()
