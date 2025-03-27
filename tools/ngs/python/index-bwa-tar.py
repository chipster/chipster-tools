# TOOL index-bwa-tar.py: "Create BWA index" ()
# INPUT input.fa TYPE FASTA
# OUTPUT output.tar
# RUNTIME python3

import tool_utils
import subprocess
from typing import Iterable
import os


def main():
    input_fa = "input.fa"

    bwa = chipster_tools_path + "/bwa/bwa"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, ".fa")

    run_process([bwa, "index", "-p", fasta_basename, input_fa])

    index_file_extensions = [
        ".amb",
        ".ann",
        ".bwt",
        ".pac",
        ".sa",
    ]

    for fe in index_file_extensions:
        index_file = fasta_basename + fe
        if not os.path.exists(index_file):
            raise RuntimeError(
                "bwa indexing of genome "
                + fasta_basename
                + " failed. Index file "
                + index_file
                + " not found"
            )

    run_process(["ls", "-lah"])

    # create tar package

    index_files = []

    for fe in index_file_extensions:
        index_file = fasta_basename + fe
        index_files.append(index_file)
            
    run_process(["tar", "-cf", "output.tar"] + index_files)

    tool_utils.write_output_definitions({
        "output.tar": fasta_basename + ".tar"
    })
    
    # save version information
    # bwa prints version to stderr when arguments are not given, but also uses exit code 1
    version = subprocess.check_output([bwa + " 2>&1 || true"], shell=True)
    version_number = str(version).split("\\n")[2].split(" ")[1]
    
    version_utils.document_version("BWA", version_number)


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
