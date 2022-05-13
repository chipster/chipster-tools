# TOOL genome-download.py: "Download genome from Ensembl" (Download genome fasta and gtf files from Ensembl.) 
# OUTPUT output.fa.lz4
# OUTPUT output.gtf.lz4
# PARAMETER species: Species TYPE STRING DEFAULT drosophila_melanogaster ()
# PARAMETER version: Version TYPE STRING DEFAULT BDGP6.32 (Genome assembly version)
# PARAMETER release: "Ensembl release" TYPE STRING DEFAULT current (Ensembl release number or "current")
# PARAMETER ftp_site: "Ensembl ftp site" TYPE ["ensembl.org", "fungi", "plants", "bacteria"] DEFAULT "ensembl.org" ()
# RUNTIME python3

import requests
from ftplib import FTP
import os.path
from urllib.parse import urlparse
from pathlib import Path
import re
import subprocess
import sys

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils


# Calls find_files and filters away some files that are typical for gtf
def find_gtf_file(dir, file):
    return find_file(dir, file, ["abinitio", ".chr.gtf.gz", ".chr_patch"])

# Finds files that are specified by its arguments first is directory and second is file (or part of a file name), rest of the arguments are expressions for grep to exclude e.g. '[\.]', ".gz", or "name". If a file is found from a directory its path is returned, otherwise raises an error
def find_file(dir, file, ignore_names):

    try:
        url = urlparse(dir)
        print("connect to FTP server " + url.netloc)
        ftp = FTP(url.netloc)
        print("login")
        ftp.login()

        print("list files in " + url.path)
        ftp.cwd(url.path)
        files = ftp.nlst()
        ftp.quit()

    except BaseException as exc:
        raise RuntimeError("unable to get directory listing: " + dir) from exc

    filtered_files = []
    
    # Filter strings away if specified
    for ftp_file in files:
        for ignore_name in ignore_names:
            if ignore_name in ftp_file:
                break
        # for-else is executed only if the for loop exited normally (i.e. without 'break')
        else:
            if file in ftp_file:
                filtered_files.append(ftp_file)

    # check empty string first, because wc will count the empty line
    if len(filtered_files) == 0:
            raise RuntimeError("Can't find file " + file + " from " + dir)

    if len(filtered_files) != 1:
            raise RuntimeError("Searching for one file, but found many: " + str(filtered_files))

    return dir + filtered_files[0]

# Returns a path to a directory where a genome is found. First argument is TYPE.
def get_dir(host, species, release, file_type):

    print("get_dir()", host, species, release, file_type)

    dir = ""

    if host.startswith("ftp://ftp.ensembl.org/"):
        if release is "current":
            # ftp://ftp.ensembl.org/pub/current_fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
            dir = host + "current_" + file_type + "/"
        else:
            # ftp://ftp.ensembl.org/pub/release-80/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.dna.toplevel.fa.gz
            dir = host + release + "/"+ file_type + "/"
    elif host.startswith("ftp://ftp.ensemblgenomes.org/"):
        # ftp://ftp.ensemblgenomes.org/pub/plants/current/fasta/populus_trichocarpa/dna/Populus_trichocarpa.JGI2.0.26.dna.toplevel.fa.gz
        dir = host + release + "/" + file_type + "/"

    if host is "ftp://ftp.ensemblgenomes.org/pub/bacteria/":
        # iterate over collections to find the right one
        i = 0
        while True:
            collection = dir + "bacteria_" + i + "_collection/"
            response = requests.get(collection)

            if response.status_code == 200:
                species_url = dir + "bacteria_" + i + "_collection/" + species + "/"
                response = requests.get(species_url)

                if response.status_code == 200:
                    dir = species_url
                    break
            else:
                # we have checked all collections
                raise RuntimeError(species + "not found from bacteria collections 0 - " + str(i - 1))
            
            i = i + 1

    else:
        dir = dir + species + "/"

    if file_type is "fasta":
        dir = dir + "dna/"

    return dir

def check_version(species, release, version, ftp_version):

    print("checking " + species + " version " + version + " from release " + release)

    if version == ftp_version:
        print("current release is " + ftp_release)
    else:
        raise RuntimeError("""

            Version of """ + species + """ in """ + release + """ is set to """ + version + """, but the version on the server is """ + ftp_version + """
            - find out the last release which still contains """ + version + """ and replace """ + release + """ with it in genomes.tsv
            - consider disabling the indexing of the old version
            - add a new row for the new version
            
            """)

def download(host, species, release, ftp_release):
    print("downloading " + species + " version " + version + " from release " + ftp_release)

    fasta_dir = get_dir(host, species, release, "fasta")
    fasta_url = find_file(fasta_dir, ".dna.toplevel.fa.gz", [])
    download_url(fasta_url, "output.fa.lz4")

    gtf_dir = get_dir(host, species, release, "gtf")
    gtf_url = find_gtf_file(gtf_dir, ".gtf.gz")
    download_url(gtf_url, "output.gtf.lz4")

    ftp_fasta_filename = os.path.basename(fasta_url)
    ftp_gtf_filename = os.path.basename(gtf_url)

    fasta_filename = tool_utils.remove_postfix(ftp_fasta_filename, '.gz') + ".lz4"
    gtf_filename = tool_utils.remove_postfix(ftp_gtf_filename, '.gz') + ".lz4"

    output_names = {
        'output.fa.lz4': fasta_filename,
        'output.gtf.lz4': gtf_filename,
    }

    tool_utils.write_output_definitions(output_names)

    print("done " + species + " version " + ftp_release)

# convert bytes to human readable string
def to_hr(byte_size: int):
    return str(round(byte_size / 1024 / 1024, 1)) + " MiB"

def download_url(url: str, file_name: str):

    print("download url " + url  + " to " + file_name)

    bash_cmd = "curl -s " + url + " | gunzip | lz4 > " + file_name
    bash_process = subprocess.run(["bash", "-c", bash_cmd])

    if bash_process.returncode != 0:
        raise RuntimeError("download failed with exit code: " + bash_process.returncode + ", command: " + bash_cmd)

    #print("download url " + url  + " to " + file_name)
    # parsed_url = urlparse(url)
    # ftp_server = parsed_url.netloc
    # ftp_dir = os.path.dirname(parsed_url.path)
    # ftp_file = os.path.basename(parsed_url.path)
    # try:
    #     print("connect to FTP server " + str(ftp_server))
    #     with FTP(ftp_server) as ftp:
    #         print("login")
    #         ftp.login()
    #         ftp.cwd(ftp_dir)
    #         total_size = ftp.size(ftp_file)
    #         print("download")
    #         size = 0
    #         chunk_count = 0
    #         with open(file_name, "wb") as f:
    #             def write_chunk(data):
    #                 f.write(data)
    #                 nonlocal size 
    #                 nonlocal chunk_count
    #                 size += len(data)
    #                 chunk_count = chunk_count + 1
    #                 if chunk_count % 10000 == 0:
    #                     print("downloaded", to_hr(size) + " / " + to_hr(total_size))
    #             ftp.retrbinary('RETR ' + ftp_file, write_chunk)        
    # except BaseException as exc:
    #     raise RuntimeError("ftp download failed", ftp_server, ftp_dir, ftp_file) from exc

    # url = url.replace("ftp://", "http://")
    # print("download url " + url  + " to " + file_name)
    # resp = requests.get(url, stream=True)
    # total = int(resp.headers.get('content-length', 0))
    # chunk_count = 0
    # downloaded_size = 0
    # with open(file_name, 'wb') as file:
    #     for data in resp.iter_content(chunk_size=1024 * 1024 * 100):
    #         size = file.write(data)
    #         downloaded_size = downloaded_size + size            
    #         # if (chunk_count % 100 == 0):
    #         print("downloaded", to_hr(downloaded_size) + " / " + to_hr(total))
    #         chunk_count = chunk_count + 1
    #     print("download completed")


def get_host(ftp_site):
    if ftp_site == "ensembl.org":
        return "ftp://ftp.ensembl.org/pub/"
    elif ftp_site == "fungi":
        return "ftp://ftp.ensemblgenomes.org/pub/fungi/"
    elif ftp_site == "plants":
        return "ftp://ftp.ensemblgenomes.org/pub/plants/"
    elif ftp_site == "bacteria":
        return "ftp://ftp.ensemblgenomes.org/pub/bacteria/"
    else:
        raise RuntimeError("unknown ftp_site " + ftp_site)

def parse_gtf_url(gtf_url):
    # name 
    # (Drosophila_melanogaster.BDGP6.80)
    basename = os.path.basename(gtf_url).replace("\\.gtf\\.gz$", "")
    name = re.sub("\.gtf\.gz$", "", basename)

    # a string before the first period 
    # (Drosophila_melanogaster)
    species_name = name.split(".", 1)[0]

    # everything after the first period 
    #(BDGP6.80)
    ftp_v_and_r = name.split(".", 1)[1]

    # remove everything after the last period 
    # (BDGP6)
    index_of_last_dot = ftp_v_and_r.rindex(".")
    ftp_version = ftp_v_and_r[:index_of_last_dot]

    # everything after the last period 
    # (80)
    ftp_release = ftp_v_and_r[index_of_last_dot + 1:]

    return species_name, ftp_version, ftp_release


# convert parameter options to real addresses
print("get_host()", ftp_site)

host = get_host(ftp_site)

print("get_dir", host, species, release)
gtf_dir = get_dir(host, species, release, "gtf")

print("find_gtf_file()", gtf_dir)
gtf_url = find_gtf_file(gtf_dir, ".gtf.gz")


species_name, ftp_version, ftp_release = parse_gtf_url(gtf_url)

print("check_version()", species, release, version, ftp_version)
check_version(species, release, version, ftp_version)

print("download()", host, species_name, ftp_release)
download(host, species, release, ftp_release)

