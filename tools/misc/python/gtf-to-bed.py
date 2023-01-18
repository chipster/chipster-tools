# TOOL gtf-to-bed.py: "Convert GTF to BED" ()
# INPUT input.gtf TYPE GENERIC
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
    output_bed = "output.bed"

    gtf2bed = chipster_tools_path + "/gtf2bed/gtf2bed.pl"

    session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
    gtf_basename = tool_utils.remove_postfix(session_input_gtf, '.gtf')

    run_bash(gtf2bed + " " + input_gtf + " > " + output_bed)
        
    output_names = {
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
