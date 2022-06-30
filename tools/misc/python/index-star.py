# TOOL index-star.py: "Create STAR index" ()
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

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils
import subprocess
from typing import Iterable
import os

def main():

    input_fa = "input.fa"
    input_gtf = "input.gtf"

    star = chipster_tools_path + "/STAR/STAR"
    
    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fa_basename = tool_utils.remove_postfix(session_input_fa, '.fa')

    run_process([star, "--runThreadN", chipster_threads_max, "--runMode", "genomeGenerate", "--genomeDir", ".", "--genomeFastaFiles", input_fa])

    print("inspect index")

    if not os.path.exists("SA") or not os.path.exists("SAindex"):
        raise RuntimeError("STAR indexing of genome " + fa_basename + " failed")
        
    run_process(["ls", "-lah"])

    # rename index files to fixed output names and define real names for the client

    outputs = [
        "SA",
        "SAindex",
        "Genome",
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "genomeParameters.txt",
    ]

    output_names = {}

    for output in outputs:
        output_names[output] = fa_basename + "/" + output

    tool_utils.write_output_definitions(output_names)
    
def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

main()
