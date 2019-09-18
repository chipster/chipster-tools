# TOOL check-fastq.R: "Check FASTQ file for errors" (This tool does some quick validity checks for a FASTQ file to spot some common problems.)
# INPUT reads: "FASTQ file" TYPE GENERIC
# OUTPUT OPTIONAL fail.log
# OUTPUT OPTIONAL pass.log


# AMS 2017-03-29

source(file.path(chipster.common.path, "zip-utils.R"))
source(file.path(chipster.common.path, "tool-utils.R"))

unzipIfGZipFile("reads")

inputnames <- read_input_definitions()

# Truncated files of lack last newline. this throws off other checks.
#system(paste("sed -i -e '$a\' reads"))

alltests <- TRUE

# Check number of lines
system("echo '#Test 1: number of lines is divisible by four' > check.log")
# Truncated files of lack last newline, so wc -l as such won't work
linenumber <- as.integer(system(paste("grep . reads | wc -l"), intern=TRUE))
if (linenumber%%4 == 0){
	system("echo PASS >> check.log")
}else{
	alltests <- FALSE
	system("echo FAIL >> check.log")
	system("echo >> check.log")
	system("echo 'Last line:' >> check.log")
	system("tail -1 reads >> check.log")
	system("echo >> check.log")
}
system("echo >> check.log")

# Check that read and quality lengths match
system("echo '#Test 2: read lengths and quality lengths match' >> check.log")
system(paste("cat reads | paste - - - - | awk -F\"\t\" '{ if (length($2) != length($4)) print $0 }' | tr '\t' '\n' > nonmatching"))
if (fileNotOk("nonmatching", minlines=1)){
	system("echo PASS >> check.log")
}else{
	alltests <- FALSE
	system("echo FAIL >> check.log")
	system("echo >> check.log")
	system("echo 'Mismatched reads:' >> check.log")
	system("cat nonmatching >> check.log")
}
system("echo >> check.log")

# Check for duplicated read names
system("echo '#Test 3: duplicated read names' >> check.log")
#system("grep \"^@\" reads |sort | uniq -c | awk '{if ($1 > 1) print $0}' > duplicates")
system("awk 'NR%4==1' reads |sort | uniq -c | awk '{if ($1 > 1) print $0}' > duplicates")
if (fileNotOk("duplicates", minlines=1)){
	system("echo PASS >> check.log")
}else{
	alltests <- FALSE
	system("echo FAIL >> check.log")
	system("echo >> check.log")
	system("echo `wc -l duplicates` >> check.log")	
}
system("echo >> check.log")

system("echo 'Test summary:' >> check.log")
system(paste("echo 'file name:", inputnames$reads, "' >>check.log"))
system(paste("echo 'reads:", floor(linenumber / 4), "' >>check.log"))
if (alltests){
	system("echo 'Overall result: PASS' >> check.log")
	system("mv check.log pass.log")
	
}else{
	system("echo 'Overall result: FAIL' >> check.log")
	system("echo >> check.log")
	system("echo 'Failures in test 1 and 2 often result from truncated files. Failure in test 3 can result e.g. from incorrect merging of FASTQ files.' >> check.log")
	system("mv check.log fail.log")
}

# Determine base name
basename <- strip_name(inputnames$reads)

# Make a matrix of output names
outputnames <- matrix(NA, nrow=2, ncol=2)
outputnames[1,] <- c("pass.log", paste(basename, "_PASS.log", sep =""))
outputnames[2,] <- c("fail.log", paste(basename, "_FAIL.log", sep =""))

# Write output definitions file
write_output_definitions(outputnames)