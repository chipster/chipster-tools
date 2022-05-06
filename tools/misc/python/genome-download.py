# TOOL genome-download.py: "Download genome from Ensembl" (Download genome fasta and gtf files from Ensembl.) 
# OUTPUT output
# PARAMETER species: Species TYPE STRING DEFAULT drosophila_melanogaster ()
# PARAMETER version: Version TYPE STRING DEFAULT BDGP6 (Genome assembly version)
# PARAMETER release: "Ensembl release" TYPE STRING DEFAULT current (Ensembl release number or "current")
# PARAMETER ftp_site: "Ensembl ftp site" TYPE ["ensembl.org", "fungi", "plants", "bacteria"] DEFAULT "ensembl.org" ()
# RUNTIME python3

import requests
from ftplib import FTP
import os.path

# Calls find_files and filters away some files that are typical for gtf
def find_gtf_file(dir, file):
    find_file(dir, file, ["abinitio", ".chr.gtf.gz", ".chr_patch"])

# Finds files that are specified by its arguments first is directory and second is file (or part of a file name), rest of the arguments are expressions for grep to exclude e.g. '[\.]', ".gz", or "name". If a file is found from a directory its path is returned, otherwise raises an error
def find_file(dir, file, ignore_names):

    try:
        ftp = FTP(dir)
        ftp.login()
        #ftp.cwd()

        dirs=[]
        files=[]
        #for item in ftp.mlsd(dir):
        for item in ftp.mlsd():
            if item[1]['type'] == 'dir':
                dirs.append(item[0])
            else:
                files.append(item[0])

        ftp.quit()

    except BaseException as exc:
        raise RuntimeError("unable to get directory listing: " + dir) from exc

    list = []

    print("files: ", files)
    print("ignore_names: ", ignore_names)
    print("file: " + file)
    
    # Filter strings away if specified
    for ftp_file in files:
        print("ftp_file: ", ftp_file)
        for ignore_name in ignore_names:
            print("ignore_name: ", ignore_name)
            if ignore_name in ftp_file:
                continue
        if file in ftp_file:
            list.append(ftp_file)

    # check empty string first, because wc will count the empty line
    if len(list) == 0:
            raise RuntimeError("Can't find file " + file + " from " + dir)

    if len(list) != 1:
            raise RuntimeError("Searching for one file, but found many: " + list)

    return dir + list[0]

# Returns a path to a directory where a genome is found. First argument is TYPE.
def get_dir(host, species, release, file_type):

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

def check_version(host, species, release):
    gtf_dir = get_dir(host, species, release, "gtf")
    gtf_url = find_gtf_file(gtf_dir, ".gtf.gz")

    # name 
    # (Drosophila_melanogaster.BDGP6.80)
    name = os.path.basename(gtf_url).replace(".gtf.gz$", "")

    # a string before the first period 
    # (Drosophila_melanogaster)
    species_name = name.split(".", 1)[0]

    # everything after the first period 
    #(BDGP6.80)
    ftp_v_and_r = name.split(".", 1)[1]

    # remove everything after the last period 
    # (BDGP6)
    ftp_version = ftp_v_and_r[ftp_v_and_r.rindex("."):]

    # everything after the last period 
    # (80)
    ftp_release = ftp_v_and_r[:ftp_v_and_r.rindex(".")]

    if version is not ftp_version:
        raise RuntimeError("""

            Version of $SPECIES in $RELEASE is set to $VERSION, but the version on the server is $FTP_VERSION
            - find out the last release which still contains $VERSION and replace $RELEASE with it in genomes.tsv
            - consider disabling the indexing of the old version
            - add a new row for the new version
            
            """)

host = ""

if ftp_site is "ensembl.org":
    host = "ftp://ftp.ensembl.org/"
elif ftp_ist is "fungi":
    host = "ftp://ftp.ensemblgenomes.org/pub/fungi/"
elif ftp_ist is "plants":
    host = "ftp://ftp.ensemblgenomes.org/pub/plants/"
elif ftp_ist is "bacteria":
    host = "ftp://ftp.ensemblgenomes.org/pub/bacteria/"

check_version(host, species, release)

