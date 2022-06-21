# TOOL index-bowtie2.py: "Create Bowtie2 index" ()
# INPUT input.fa TYPE GENERIC
# INPUT input.gtf TYPE GENERIC
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
    input_gtf = "input.gtf"
    bowtie2_path = chipster_tools_path + "/bowtie2"
    bowtie2_build = chipster_tools_path + "/bowtie2/bowtie2-build"
    bowtie2_inspect = chipster_tools_path + "/bowtie2/bowtie2-inspect"
    tophat = chipster_tools_path + "/tophat2/tophat"

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, '.fa')

    session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
    gtf_basename = tool_utils.remove_postfix(session_input_gtf, '.gtf')

    # gtf_basename contains the genome version and the Ensembl version 
    # (...GRCh38.81), but in the latest releases, the fasta_basename contains
    # only the genome version (...GRCh38). This makes sense, because fasta files 
    # different Ensembl releases should be identical after the patches are removed. 
    # The Bowtie index is based only on fasta, so it in theory it should be 
    # named after the fasta. However, the  Tophat index is based also on gtf and the
    # Tophat script has to find both indexes, so let's use gtf_basename also for the 
    # Bowtie index.
    run_process([bowtie2_build, input_fa, gtf_basename])

    print("inspect index")
    inspect_output = "inspect_output.txt"
    run_bash(bowtie2_inspect + " -n " + gtf_basename + " > " + inspect_output)

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
        raise RuntimeError("bowtie2 indexing of genome " + gtf_basename + " failed. Chromosomes in fasta: " + str(fasta_chr) + ", chromosomes in index: " + str(index_chr))
    
    run_process(["ls", "-lah"])

    # symlink the fasta for the tophat, otherwise it will rebuild it
    run_process(["ln", "-s", input_fa, gtf_basename + ".fa"])
    # tophat seems to take the output name from the gtf file
    run_process(["ln", "-s", input_gtf, gtf_basename + ".gtf"])    
    run_bash("export PATH=\"$PATH:" + bowtie2_path + "\"; " + tophat + " --GTF " + gtf_basename + ".gtf --transcriptome-index " + gtf_basename + " " + gtf_basename)

    run_process(["ls", "-lah"])

    # Why does tophat create a dir named $gtf_basename?
    run_process(["ls", "-lah", gtf_basename])

    # rename index files to fixed output names and define real names for the client

    output_names = {}

    # bowtie2 index files
    for file in os.listdir("."):
        if file.endswith(".bt2"):
            output_name = file.replace(gtf_basename, "output") + str(len(output_names))
            os.rename(file, output_name)
            output_names[output_name] = "bowtie2/" + file

    # tophat index files
    for file in os.listdir(gtf_basename):
        output_name = file.replace(gtf_basename, "output") + str(len(output_names))
        os.rename(gtf_basename + "/" + file, output_name)
        output_names[output_name] = "tophat2/" + file

    tool_utils.write_output_definitions(output_names)
    
def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

main()

