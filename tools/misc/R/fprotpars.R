# TOOL fprotpars.R: "Protein parsimony algorithm" (Protein parsimony algorithm of the PHYLIP package)
# INPUT sequence: sequence TYPE GENERIC
# OUTPUT OPTIONAL protpars.out.txt
# OUTPUT OPTIONAL protpars.tree.txt
# OUTPUT OPTIONAL fprotpars.log
# PARAMETER OPTIONAL njumble: "Number of times to randomise" TYPE INTEGER FROM 0 DEFAULT 0 (Number of times to randomise)
# PARAMETER OPTIONAL seed: "Random number seed between 1 and 32767 (must be odd\)" TYPE INTEGER FROM 1 TO 32767 DEFAULT 1 (Random number seed between 1 and 32767 (must be odd\))
# PARAMETER OPTIONAL outgrno: "Species number to use as outgroup" TYPE INTEGER FROM 0 DEFAULT 0 (Species number to use as outgroup)
# PARAMETER OPTIONAL dothreshold: "Use threshold parsimony" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Use threshold parsimony)
# PARAMETER OPTIONAL threshold: "Threshold value" TYPE DECIMAL FROM 1 DEFAULT 1 (Threshold value)
# PARAMETER OPTIONAL whichcode: "Use which genetic code" TYPE [U: Universal, M: Mitochondrial, V: "Vertebrate mitochondrial", F: "Fly mitochondrial", Y: "Yeast mitochondrial"] FROM 1 TO 1 DEFAULT U (Use which genetic code)
# PARAMETER OPTIONAL trout: "Write out trees to tree file" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Write out trees to tree file)
# PARAMETER OPTIONAL printdata: "Print data at start of run" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Print data at start of run)
# PARAMETER OPTIONAL treeprint: "Print out tree" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Print out tree)
# PARAMETER OPTIONAL stepbox: "Print steps at each site" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Print steps at each site)
# PARAMETER OPTIONAL ancseq: "Print sequences at all nodes of tree" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT N (Print sequences at all nodes of tree)
# PARAMETER OPTIONAL dotdiff: "Use dot differencing to display results" TYPE [<undefined>: " ", Y: Yes, N: No] DEFAULT Y (Use dot differencing to display results)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the analysis run.)

emboss.path <- file.path(chipster.tools.path, "emboss", "bin")

source(file.path(chipster.common.lib.path, "zip-utils.R"))
unzipIfGZipFile("sequence")

# check sequece file type
inputfile.to.check <- ("sequence")
sfcheck.binary <- file.path(chipster.module.path, "/shell/sfcheck.sh")
sfcheck.command <- paste(sfcheck.binary, emboss.path, inputfile.to.check)
str.filetype <- system(sfcheck.command, intern = TRUE)

if (str.filetype == "Not an EMBOSS compatible sequence file") {
    stop("CHIPSTER-NOTE: Your input file is not a sequence file that is compatible with the tool you try to use")
}

# count the query sequeces
seqcount.exe <- file.path(emboss.path, "seqcount sequence -filter")
str.queryseq <- system(seqcount.exe, intern = TRUE)
num.queryseq <- as.integer(str.queryseq)

# round(num.queryseq)

if (num.queryseq > 100000) {
    stop(paste("CHIPSTER-NOTE: Too many query sequences. Maximun is 100000 but your file contains ", num.queryseq))
}

emboss.binary <- file.path(emboss.path, "fprotpars")
emboss.parameters <- paste("-auto -progress N")
emboss.parameters <- paste(emboss.parameters, "-sequence sequence")
emboss.parameters <- paste(emboss.parameters, "-outfile protpars.out.txt")
emboss.parameters <- paste(emboss.parameters, "-outtreefile protpars.tree.txt")
emboss.parameters <- paste(emboss.parameters, "-njumble", njumble)
emboss.parameters <- paste(emboss.parameters, "-seed", seed)
emboss.parameters <- paste(emboss.parameters, "-outgrno", outgrno)
emboss.parameters <- paste(emboss.parameters, "-dothreshold", dothreshold)
emboss.parameters <- paste(emboss.parameters, "-threshold", threshold)
emboss.parameters <- paste(emboss.parameters, "-whichcode", whichcode)
emboss.parameters <- paste(emboss.parameters, "-trout", trout)
emboss.parameters <- paste(emboss.parameters, "-printdata", printdata)
emboss.parameters <- paste(emboss.parameters, "-treeprint", treeprint)
emboss.parameters <- paste(emboss.parameters, "-stepbox", stepbox)
emboss.parameters <- paste(emboss.parameters, "-ancseq", ancseq)
emboss.parameters <- paste(emboss.parameters, "-dotdiff", dotdiff)

command.full <- paste(emboss.binary, emboss.parameters, " >> fprotpars.log 2>&1")
echo.command <- paste('echo "', command.full, ' "> fprotpars.log')
system(echo.command)

system(command.full)

system("ls -l >> fprotpars.log")

if (save_log == "no") {
    system("rm -f fprotpars.log")
}
