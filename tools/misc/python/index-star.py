# TOOL index-star.py: "Create STAR index without tar package" ()
# INPUT input.fa TYPE GENERIC
# OUTPUT SA
# OUTPUT SAindex
# OUTPUT Genome
# OUTPUT chrLength.txt
# OUTPUT chrName.txt
# OUTPUT chrNameLength.txt
# OUTPUT chrStart.txt
# OUTPUT genomeParameters.txt
# RUNTIME python3
# SLOTS 5

import tool_utils
import subprocess
from typing import Iterable
import os


def main():
    input_fa = "input.fa"

    star = chipster_tools_path + "/STAR/STAR"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fa_basename = tool_utils.remove_postfix(session_input_fa, ".fa")

    run_process(
        [
            star,
            "--runThreadN",
            chipster_threads_max,
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            ".",
            "--genomeFastaFiles",
            input_fa,
        ]
    )

    print("inspect index")

    if not os.path.exists("SA") or not os.path.exists("SAindex"):
        raise RuntimeError("STAR indexing of genome " + fa_basename + " failed")

    run_process(["ls", "-lah"])


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
