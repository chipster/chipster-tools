# TOOL index-bowtie-tar.py: "Create Bowtie index" ()
# INPUT input.fa TYPE FASTA
# OUTPUT output.tar
# RUNTIME python3

import tool_utils
import subprocess
from typing import Iterable
import os


def main():
    input_fa = "input.fa"

    bowtie_build = chipster_tools_path + "/bowtie/bowtie-build"
    bowtie_inspect = chipster_tools_path + "/bowtie/bowtie-inspect"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, ".fa")

    # our bowtie is too old for --threads
    run_process(
        [bowtie_build, "--threads", chipster_threads_max, input_fa, fasta_basename]
    )

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
        raise RuntimeError(
            "bowtie indexing of genome "
            + fasta_basename
            + " failed. Chromosomes in fasta: "
            + fasta_chr
            + ", chromosomes in index: "
            + index_chr
        )

    run_process(["ls", "-lah"])

    # create tar package

    index_files = []

    for file in os.listdir("."):
        # .ebwt or .ebwtl
        if ".ebwt" in file:
            index_files.append(file)
            
    run_process(["tar", "-cf", "output.tar"] + index_files)

    tool_utils.write_output_definitions({
        "output.tar": fasta_basename + ".tar"
    })
    
    # save version information
    version = subprocess.check_output([bowtie_build, "--version"])
    version_number = str(version).split("\\n")[0].split(" ")[2]
    
    version_utils.document_version("Bowtie", version_number)


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
