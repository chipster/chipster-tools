# TOOL index-hisat2-tar.py: "Create Hisat2 index" ()
# INPUT input.fa TYPE FASTA
# INPUT OPTIONAL input.gtf TYPE GTF
# OUTPUT output.tar
# RUNTIME python3
# SLOTS 9

import tool_utils
import subprocess
from typing import Iterable
import os
import os.path

def main():
    input_fa = "input.fa"
    input_gtf = "input.gtf"

    extract_splice_sites = chipster_tools_path + "/hisat2/extract_splice_sites.py"
    extract_exons = chipster_tools_path + "/hisat2/extract_exons.py"
    hisat2_build = chipster_tools_path + "/hisat2/hisat2-build"
    hisat2_inspect = chipster_tools_path + "/hisat2/hisat2-inspect"
    
    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    genome_basename = tool_utils.remove_postfix(session_input_fa, ".fa")
    
    if not os.path.isfile("input.gtf"):
        print("not using splice site and exon annotations, because gtf file was not given")
        run_process([hisat2_build, "-p", chipster_threads_max, input_fa, "output"])

    else:
        # Run hisat2-build with splice site and exon annotations, this consumes lot of ram: ~200G for human
        # rat would require 0.5 TB RAM with --ss and --exon https://www.biostars.org/p/344449/
        
        session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
        genome_basename = tool_utils.remove_postfix(session_input_gtf, ".gtf")
        
        # Create a splice site file
        run_bash(extract_splice_sites + " " + input_gtf + " > " + genome_basename + ".ss")

        # Create an exon file
        run_bash(extract_exons + " " + input_gtf + " > " + genome_basename + ".exon")

        run_process(
            [
                hisat2_build,
                "-p",
                chipster_threads_max,
                "--ss",
                genome_basename + ".ss",
                "--exon",
                genome_basename + ".exon",
                input_fa,
                "output",
            ]
        )

    print("inspect index")
    inspect_output = "inspect_output.txt"
    run_bash(hisat2_inspect + " -n output > " + inspect_output)

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
            "hisat2 indexing of genome "
            + genome_basename
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
        if file.startswith("output"):
            file_extension = file.replace("output", "")
            new_name = genome_basename + file_extension
            os.rename(file, new_name)
            index_files.append(new_name)
            
    run_process(["tar", "-cf", "output.tar"] + index_files)

    tool_utils.write_output_definitions({
        "output.tar": genome_basename + ".tar"
    })
    
    # save version information
    version = subprocess.check_output([hisat2_build, "--version"])
    version_number = str(version).split("\\n")[0].split(" ")[2]
    
    version_utils.document_version("HISAT2", version_number)

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
