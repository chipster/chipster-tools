# TOOL index-tophat2-tar.py: "Create TopHat2 index" ()
# INPUT input.fa TYPE FASTA
# INPUT input.gtf TYPE GTF
# OUTPUT output.tar
# RUNTIME python
# SLOTS 2

# use old runtime "python", because tophat2 is written in Python 2

import tool_utils_python2 as tool_utils
import subprocess
import os

def main():
    input_fa = "input.fa"
    input_gtf = "input.gtf"
    bowtie2_path = chipster_tools_path + "/bowtie2"
    bowtie2_build = chipster_tools_path + "/bowtie2/bowtie2-build"
    bowtie2_inspect = chipster_tools_path + "/bowtie2/bowtie2-inspect"
    tophat2 = chipster_tools_path + "/tophat2/tophat"
    
    # add bowtie2 to PATH so that tophat can find it
    env = os.environ.copy()
    # relative path, but works if we don't change working dir
    env["PATH"] = env['PATH'] + ":" + bowtie2_path

    session_input_fa = tool_utils.read_input_definitions()[input_fa]
    session_input_gtf = tool_utils.read_input_definitions()[input_gtf]
    fasta_basename = tool_utils.remove_postfix(session_input_fa, ".fa")    
    gtf_basename = tool_utils.remove_postfix(session_input_gtf, ".gtf")

    run_process(
        [bowtie2_build, "--threads", chipster_threads_max, input_fa, gtf_basename]
    )
    
    # tophat2 will generate the fasta from index if it doesn't find it
    run_process(["ln", "-s", input_fa, gtf_basename + ".fa"])
        
    run_process(
        [tophat2, "-G", input_gtf, "--transcriptome-index", gtf_basename, gtf_basename], env
    )
    
    print("inspect bowtie2 index")
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
        raise RuntimeError(
            "bowtie2 indexing of genome "
            + fasta_basename
            + " failed. Chromosomes in fasta: "
            + str(fasta_chr)
            + ", chromosomes in index: "
            + str(index_chr)
        )

    run_process(["ls", "-lah"])
    
    run_process(["ls", "-lah", gtf_basename])

    # create tar package

    index_files = []

    # tophat2 index files
    for file in os.listdir(gtf_basename):
        new_name = file.replace("input", gtf_basename)
        os.rename(gtf_basename + "/" + file, gtf_basename + "/" + new_name)
        # without directory, becaue tar will run in that directory
        index_files.append(new_name)
            
    run_process(["tar", "-cf", "../output.tar"] + index_files, cwd=gtf_basename)

    tool_utils.write_output_definitions({
        "output.tar": fasta_basename + ".tar"
    })
    
    # save version information
    version = subprocess.check_output([bowtie2_build, "--version"])
    # backward slash "\" has to be escaped with another backward slash in python3, but not here in python2
    version_number = str(version).split("\n")[0].split(" ")[2]
    
    version_utils.document_version("Bowtie2", version_number)
    
    version = subprocess.check_output([tophat2, "--version"])
    version_number = str(version).split("\n")[0].split(" ")[1]
    
    version_utils.document_version("TopHat2", version_number)


def run_bash(cmd):
    run_process(["bash", "-c", cmd])


def run_process(cmd, env=None, cwd="."):
    # there is no run() in python2
    #process = subprocess.run(cmd)
    
    if not env:
        env = os.environ.copy()
        
    returncode = subprocess.call(cmd, env=env, cwd=cwd)
        
    if returncode != 0:
        raise RuntimeError(
            "process failed with return code: "
            + str(process.returncode)
            + ", command: "
            + str(cmd)
        )


main()
