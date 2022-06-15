# TOOL index-gtf.py: "Create GTF indexes" ()
# INPUT input.gtf TYPE GENERIC
# OUTPUT output.gtf
# OUTPUT output.bed
# RUNTIME python3

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils
import subprocess
from typing import Iterable
import os

def main():

    input_gtf = "input.gtf"
    output_gtf = "output.gtf"
    output_bed = "output.bed"

    bwa = chipster_tools_path + "/bwa/bwa"

    python = chipster_tools_path + "/Python-2.7.12/bin/python"
    dexseq_prepare_annotation = chipster_tools_path + "/dexseq-exoncounts/dexseq_prepare_annotation.py"
    gtf2bed = chipster_tools_path + "/gtf2bed/gtf2bed.pl"

    session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
    gtf_basename = tool_utils.remove_postfix(session_input_gtf, '.gtf')

    run_process([python, dexseq_prepare_annotation, input_gtf, output_gtf])

    run_bash(gtf2bed + " " + input_gtf + ".gtf > " + output_bed)
        
    output_names = {
        output_gtf: gtf_basename + ".DEXSeq.gtf",
        output_bed: gtf_basename + ".bed",
    }

    tool_utils.write_output_definitions(output_names)
    
def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

main()
