# TOOL index-dexseq.py: "Create DEXSeq index without tar package" ()
# INPUT input.gtf TYPE GENERIC
# OUTPUT output.gtf
# RUNTIME python3

import tool_utils
import subprocess
from typing import Iterable
import os


def main():
    input_gtf = "input.gtf"
    output_gtf = "output.gtf"

    # TODO update dexseq_preprare_annotation.py to get rid of old Python 2
    python = chipster_tools_path + "/Python-2.7.12/bin/python"
    dexseq_prepare_annotation = (
        chipster_tools_path + "/dexseq-exoncounts/dexseq_prepare_annotation.py"
    )

    session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
    gtf_basename = tool_utils.remove_postfix(session_input_gtf, ".gtf")

    run_process([python, dexseq_prepare_annotation, input_gtf, output_gtf])

    output_names = {
        output_gtf: gtf_basename + ".DEXSeq.gtf",
    }

    tool_utils.write_output_definitions(output_names)


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
