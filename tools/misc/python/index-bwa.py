# TOOL index-bwa.py: "Create BWA index" ()
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

    bwa = chipster_tools_path + "/bwa/bwa"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, '.fa')

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
            raise RuntimeError("bwa indexing of genome " + fasta_basename + " failed. Index file " + index_file + " not found")
    
    run_process(["ls", "-lah"])

    # rename index files to fixed output names and define real names for the client

    output_names = {}

    for fe in index_file_extensions:
        index_file = fasta_basename + fe
        output_name = index_file.replace(fasta_basename, "output") + str(len(output_names))
        os.rename(index_file, output_name)
        output_names[output_name] = index_file

    tool_utils.write_output_definitions(output_names)
    
def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

main()
