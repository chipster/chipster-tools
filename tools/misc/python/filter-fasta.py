# TOOL filter-fasta.py: "Filter fasta" (Keep only primary chromosomes in fasta. Chromosomes that have karyotype data in Ensembl are considered as the primary chromosomes.)
# INPUT input.fa TYPE GENERIC
# OUTPUT output.fa
# RUNTIME python3

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils
import subprocess
from typing import Iterable
import os
import requests

def main():

    input_fa = "input.fa"
    output_fa = "output.fa"

    samtools = chipster_tools_path + "/samtools-1.2/samtools"

    session_input_fa = tool_utils.read_input_definitions()[input_fa].replace(".dna.toplevel", "")

    karyotype_chr = get_karyotype_chromosomes(session_input_fa)

    print("karyotype chromosomes", karyotype_chr)

    print("index original genome")
    run_process([samtools, "faidx", input_fa])

    fasta_chr = []

    with open(input_fa + ".fai") as file:
        for line in file:
            cols = line.split("\t")
            fasta_chr.append(cols[0])

    fasta_chr_set = set(fasta_chr)

    print("fasta chromosomes", fasta_chr_set)

    if karyotype_chr == None or len(karyotype_chr) == 0:
        print("no karyotype chromosomes, keeping all")

        #run_bash("bgzip --decompress --stdout " + input_fa + " > " + output_fa)
        os.rename(input_fa, output_fa)
    
    else:

        for chromosome in karyotype_chr:
            if chromosome in fasta_chr_set:
                print("copy chromosome", chromosome)
                run_bash(samtools + " faidx " + input_fa + " " + chromosome + " >> " + output_fa)
            else:
                # if the chromosome isn't found, samtools returns an emtpy sequence
		        # we don't want to add new chromosomes to old fastas, so check that the 
		        # chromosome exists before running samtools 
                print("no chromosome " + chromosome + " in fasta, skipping")

    # write better file names for client
    

    output_names = {
        output_fa: session_input_fa,
    }

    tool_utils.write_output_definitions(output_names)

def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])

def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError("process failed with return code: " + str(process.returncode) + ", command: " + str(cmd))

def get_karyotype_chromosomes(name: str):

    # a string before the first period 
    # (drosophila_melanogaster)
    species_name = name.split(".", 1)[0].lower()

    server = "https://rest.ensembl.org"
    ext = "/info/assembly/" + species_name + "?"
 
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 
    if not r.ok:
        r.raise_for_status()

    decoded = r.json()

    karyotype = decoded["karyotype"]

    return set(karyotype)

main()
