# TOOL filter-fasta.py: "Filter fasta" (Keep only primary chromosomes in fasta. Chromosomes that have karyotype data in Ensembl are considered as the primary chromosomes.)
# INPUT input.fa.gz TYPE GENERIC
# INPUT OPTIONAL coord_system.txt TYPE GENERIC
# INPUT OPTIONAL seq_region.txt TYPE GENERIC
# INPUT OPTIONAL karyotype.txt TYPE GENERIC
# OUTPUT output.fa
# PARAMETER ftp_site: "Ensembl ftp site" TYPE ["ensembl.org", "fungi", "plants", "bacteria"] DEFAULT "ensembl.org" ()
# RUNTIME python3

# add the tools dir to path, because __main__ script cannot use relative imports
sys.path.append(os.getcwd() + "/../toolbox/tools")
from common.python import tool_utils

def main():

    coord_input = "coord_system.txt"
    seq_input = "seq_region.txt"
    karyotype_input = "karyotype.txt"

    karyotype_chr = get_karyotype_chromosomes(coord_input, seq_input, karyotype_input)

    print("karyotype chromosomes", karyotype_chr)

def get_karyotype_chromosomes(coord_input: str, seq_input: str, karyotype_input: str):

    for file in [ coord_input, seq_input, karyotype_input]:
        if not os.path.exists(file):
            raise RuntimeError("cannot filter fasta without mysql file", file)

    coord_chromosomes = []

    with open(coord_input) as lines:
        for line in lines:
            tabs = line.split("\t")
            # "primary_assembly" in drosophila_melanogaster
            if tabs[2] == "chromosome" or tabs[2] == "primary_assembly":
                # there are multiple rows for different assembly versions (human)
                # but that shouldn't matter as long as only one of them is used in karyotype data
                coord_chromosomes.append(tabs[1])

    if len(coord_chromosomes) == 0:
        raise RuntimeError("no chromosomes were found from " + coord_input)

    coord_chr_set = set(coord_chromosomes)

    chr_map = {}

    with open(seq_input) as lines:
        for line in lines:
            tabs = line.split("\t")
            if tabs[2] in coord_chr_set:
                if tabs[0] in chr_map:
                    raise RuntimeError("duplicate id", tabs[0])
                chr_map[tabs[0]] = tabs[1]

    karyotype_chr = set()

    with open(karyotype_input) as lines:
        for line in lines:
            tabs = line.split("\t")
            if tabs[1] in chr_map:
                chr_name = chr_map[tabs[1]]
                karyotype_chr.add(chr_name)

    return karyotype_chr

main()
