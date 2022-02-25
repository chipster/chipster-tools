# TOOL bbduk-trim.R: "Trim QuantSeq reads using BBDuk" (Remove the adapter contamination, polyA read through, and low quality tails.)
# INPUT reads.fq: "FASTQ file" TYPE GENERIC 
# OUTPUT trimmed.fq.gz
# PARAMETER OPTIONAL k: "Kmer length" TYPE INTEGER FROM 1 DEFAULT 13 (Kmer length used for finding contaminants. Contaminants shorter than k will not be found. k must be at least 1.)
# PARAMETER OPTIONAL ktrim: "Trim reads" TYPE [f: "Don't trim", r: "Trim to the right", l: "Trim to the left"] DEFAULT r (Trim reads to remove bases matching reference kmers.)
# PARAMETER OPTIONAL mink: "Mink" TYPE INTEGER DEFAULT 5 (Look for shorter kmers at read tips down to this length, when k-trimming or masking. 0 means disabled.)
# PARAMETER OPTIONAL qtrim: "Trim read ends" TYPE [rl: (trim both ends), f: (neither end), r: (right end only), l: (left end only), w: (sliding window)] DEFAULT rl (Trim read ends to remove bases with quality below trimq.Performed AFTER looking for kmers.)
# PARAMETER OPTIONAL trimq: "" TYPE DECIMAL DEFAULT 10 (Regions with average quality BELOW this will be trimmed, if qtrim is set to something other than f.  Can be a floating-point number like 7.3.)
# PARAMETER OPTIONAL minlength: "Minimum length" TYPE INTEGER DEFAULT 20 (Reads shorter than this after trimming will be discarded. Pairs will be discarded if both are shorter.)

source(file.path(chipster.common.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

# Uncompress input if compressed
unzipIfGZipFile("reads.fq")

# File paths
java.path <- paste("PATH=", c(file.path(chipster.tools.path,"rtg","jre","bin")),":$PATH",sep="",collapse="")
bbduk.binary <- c(file.path(chipster.tools.path, "bbmap", "bbduk.sh"))
polya.file <- c(file.path(chipster.tools.path, "bbmap", "resources", "polyA.fa.gz"))
trueseq.file <- c(file.path(chipster.tools.path, "bbmap", "resources", "truseq.fa.gz"))

# Run BBduk
system(paste(bbduk.binary,"-v 2> version.tmp"))
version <- system("grep Version version.tmp",intern = TRUE)
documentVersion("BBDuk",version)
bbduk.options <- paste("in=reads.fq ", "out=trimmed.fq ", "ref=", polya.file,",",trueseq.file, " k=",k," ktrim=",ktrim," mink=",mink, " qtrim=",qtrim, " trimq=",trimq," minlength=",minlength,sep="",collapse="")
bbduk.command <- paste(bbduk.binary, bbduk.options)
documentCommand(bbduk.command)
#system(bbduk.command)
runExternal(bbduk.command, java.path)

# Compress output
system("gzip trimmed.fq")

# Output names
inputnames <- read_input_definitions()
basename <- strip_name(inputnames$reads.fq)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=1, ncol=2)
outputnames[1,] <- c("trimmed.fq.gz", paste(basename, "_trimmed.fq.gz", sep=""))

# Write output definitions file
write_output_definitions(outputnames)