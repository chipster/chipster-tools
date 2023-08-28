# TOOL enafetch.R: "Retrieve data from ENA database" (Retrieve data from the European Nucleotide Archive, ENA. ENA contains high throughput sequencing data from many sources, including the Sequence Read Archive, SRA.)
# OUTPUT OPTIONAL {...}.gz: "gzipped files"
# OUTPUT OPTIONAL {...}.fasta: "fasta files"
# OUTPUT OPTIONAL {...}.dat: "embl files"
# OUTPUT OPTIONAL {...}.txt: "text files"
# OUTPUT OPTIONAL {...}.xml: "xml files"
# OUTPUT OPTIONAL {...}.bam: "bam files"
# OUTPUT OPTIONAL {...}.cram: "cram files"
# OUTPUT OPTIONAL {...}.crai: "cram files"
# PARAMETER entry_id: "ENA ID" TYPE STRING DEFAULT "entry" (Give the ID of the ENA dataset to be retrieved. For example SRR000021)
# PARAMETER format: "Data format" TYPE [default: "Default", fastq: "FASTQ", fasta: "FASTA", embl: "EMBL formatted sequence" ] DEFAULT default (Define the format for retrieved data. Please use the FASTQ format for reads.)
# IMAGE comp-r-4-2-3-enabrowsertools
# RUNTIME R-4.2.3
# TOOLS_BIN ""

source(file.path(chipster.common.path, "tool-utils.R"))

# use enabrowsertools from the image
ena.path <- file.path("/opt/chipster/tools", "enabrowsertools/python3")
ena.binary <- file.path(ena.path, "enaDataGet")

if (format != "default") {
    command.full <- paste(ena.binary, "-f", format, entry_id)
} else {
    command.full <- paste(ena.binary, entry_id)
}

runExternal(command.full)

ena.version <- system(paste(ena.binary, "--version"), intern = TRUE)
documentVersion("enaBrowserTools", ena.version)
documentCommand(command.full)

system("find")

file.found <- FALSE

# fastq files are stored in a subdirectory. Move them to current directory (but keep them gzipped)
# cram files are stored in a subdirectory. Move them to current diretory. Those are not gzipped.
for (file in list.files(recursive = TRUE, include.dirs = TRUE)) {
    if (dir.exists(file)) {
        print(paste("skip directory", file))
    } else {
        if (any(endsWith(file, c(".gz", ".fasta", ".dat", ".txt", ".xml", ".bam", ".cram", ".crai")))) {
            file.found <- TRUE
            filename <- basename(file)

            if (file == filename) {
                print(paste("no need to move file", file))
            } else {
                # move file from the subfolder to the job working diretory
                print(paste("move file", file, " to ", filename))
                file.rename(file, filename)
            }

            # gunzip fasta and dat files if needed
            # .fastq.gz files shouldn't be extracted
            if (any(endsWith(file, c(".fasta.gz", ".dat.gz")))) {
                # TODO find suitable ENA ID to test this
                print(paste("extract file", filename))
                runExternal(paste("gunzip", filename))
            }
        } else {
            print(paste("skip unknown file type:", file))
        }
    }
}

system("find")

if (!file.found) {
    stop("CHIPSTER-NOTE: no results, click 'Show screen output' to see details")
}

## check storage resources and needs
# freespace <- system(" df ./ | awk '{print $4}' ", intern = TRUE )
# sra.binary <- file.path(sra.path, "vdb-dumb")
# commad.full <- paste(sra.binary, ' --info ', entry_id, '| grep -w size | tr -d ","', " | awk '{print $3}'")
# srasize <-system(command.full)
# space_needed = 4*srasize

# if ( freespace < space_needed ){
# 	stop("CHIPSTER-NOTE: THe SRA entry you are trying to retrieve is too larage for the chipster server")
# }
