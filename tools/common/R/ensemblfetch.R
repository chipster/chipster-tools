# TOOL ensemblfetch.R: "Retrieve data for a given organism in Ensembl" (Retrieves the genomic, cDNA or protein dataset of the given species from the Ensembl data bases.)
# OUTPUT OPTIONAL ensemblfetch.fasta
# OUTPUT OPTIONAL ensemblfetch.gtf
# OUTPUT OPTIONAL ensemblfetch_species.tsv
# OUTPUT OPTIONAL ensemblfetch.log
# PARAMETER OPTIONAL species_menu: "Ensembl species" TYPE [species: "Use the text box below", Homo_sapiens: "Human", Mus_musculus: "Mouse", Danio_rerio: "Zebrafish", Mus_spretus: "Algerian_mouse", Ciona_intestinalis: "C.intestinalis", Caenorhabditis_elegans: "Caenorhabditis_elegans", Gallus_gallus: "Chicken", Pan_troglodytes: "Chimpanzee", Gadus_morhua: "Cod", Bos_taurus: "Cow", Canis_lupus: "Dog", Drosophila_melanogaster: "Fruitfly", Cavia_porcellus: "Guinea_Pig", Equus_caballus: "Horse", Sus_scrofa: "Pig", Oryctolagus_cuniculus: "Rabbit", Rattus_norvegicus: "Rat", Saccharomyces_cerevisiae: "Saccharomyces_cerevisiae", Ovis_aries: "Sheep", Gasterosteus_aculeatus: "Stickleback", Xenopus_tropicalis: "Xenopus" ] DEFAULT species ( Quick selection meny of commonly used Ensembl species.)
# PARAMETER OPTIONAL species: "Species name" TYPE STRING (Then latin name of the species for which the data is retrieved. Note that you should use under score: _ in stead of the space character in the species name. For example homo_sapiens)
# PARAMETER OPTIONAL type: "Data to retrieve" TYPE [dna: "Genomic DNA", cdna: "cDNA transcripts", pep: "Protein sequences", gtf: "GTF file" ] DEFAULT dna (Sequence data type to retrieve)
# PARAMETER OPTIONAL names: "List the available species names" TYPE [ yes: "All species", nonbac: "List non-bacterial species", no: No] DEFAULT "no" (List the available species names)
# PARAMETER OPTIONAL save_log: "Collect a log file about the ensemblfetch run" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the enseblfetch run.)

# K.M 28.10.2013

# enseblfecth settings
ensemblfetch.binary <- file.path(chipster.module.path, "../admin/shell/ensemblfetch.sh ")

if (species_menu != "species") {
    species <- species_menu
}


# remove spaces from sepcies name and convert all letters to lower keys
species <- gsub(" ", "_", species)
species <- tolower(species)

if (names == "yes" || names == "nonbac" || nchar(species) < 3) {
    if (names == "nonbac") {
        command.to_run <- paste(ensemblfetch.binary, " -names -bacteria no > ensemblfetch_species.tsv")
        echo.command <- paste('echo "', command.to_run, ' "> ensemblfetch.log')
        system(echo.command)
        system(command.to_run)
    } else {
        command.to_run <- paste(ensemblfetch.binary, " -names > ensemblfetch_species.tsv")
        echo.command <- paste('echo "', command.to_run, ' "> ensemblfetch.log')
        system(echo.command)
        system(command.to_run)
    }
} else {
    command.to_run <- paste(ensemblfetch.binary, " -type ", type, " -out ensemblfetch.fasta ", species, " > ensemblfetch.log")
    system(command.to_run)
    if (type == "gtf") {
        system("mv ensemblfetch.fasta ensemblfetch.gtf")

        # Handle output names
        source(file.path(chipster.common.lib.path, "tool-utils.R"))
        outputnames <- matrix(NA, nrow = 1, ncol = 2)
        outputnames[1, ] <- c("ensemblfetch.gtf", paste(species, ".gtf", sep = ""))
        # Write output definitions file
        write_output_definitions(outputnames)
    } else {
        # Handle output names
        source(file.path(chipster.common.lib.path, "tool-utils.R"))
        outputnames <- matrix(NA, nrow = 1, ncol = 2)
        outputnames[1, ] <- c("ensemblfetch.fasta", paste(species, "_", type, ".fasta", sep = ""))
        # Write output definitions file
        write_output_definitions(outputnames)
    }
}



system("ls -l >> ensemblfetch.log")

if (save_log == "no") {
    system("rm -f ensemblfetch.log")
}
