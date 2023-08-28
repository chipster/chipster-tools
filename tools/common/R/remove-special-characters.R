# TOOL remove-special-characters.R: "Remove special characters from a file" (Removes potetially problematic characters from a file. Please note that this tol removes all instances of the characters and does not do any context checking. In some case it may render the file unusable.)
# INPUT input: "File to process" TYPE GENERIC
# OUTPUT OPTIONAL output
# OUTPUT OPTIONAL output.gz
# PARAMETER remove: "Charackters to remove" TYPE [all: "all", ex:"\!", pi: "\|",hy: "-", us: "_", br: "\(\)", cb: "\{\}", sb: "\[\]", ac: "\'"] DEFAULT all (Characters to remove.)

source(file.path(chipster.common.lib.path, "tool-utils.R"))
source(file.path(chipster.common.path, "zip-utils.R"))

gzip <- FALSE
if (isGZipFile("input")) {
   gzip <- TRUE
   unzipIfGZipFile("input")
}
file.copy("input", "output")

# Remove characters
# Exlamation mark "!"
if (remove == "ex" || remove == "all") {
   system("sed -i s/\\!//g output")
}
# Pipe character "|"
if (remove == "pi" || remove == "all") {
   system("sed -i s/\\|//g output")
}
# Hyphen "-"
if (remove == "hy" || remove == "all") {
   system("sed -i s/-//g output")
}
# Underscore "_"
if (remove == "us" || remove == "all") {
   system("sed -i s/_//g output")
}
# Brackets "()"
if (remove == "br" || remove == "all") {
   system("cat output | sed s/\\(//g | sed s/\\)//g > output_tmp")
   file.copy("output_tmp", "output", overwrite = TRUE)
}
# Curly brackets "{}"
if (remove == "cb" || remove == "all") {
   system("cat output | sed s/\\{//g | sed s/\\}//g > output_tmp")
   file.copy("output_tmp", "output", overwrite = TRUE)
}
# Square brackets "[]"
if (remove == "sb" || remove == "all") {
   # system("cat output | sed s/\\\[//g | sed s/\\\]//g > output_tmp")
   system("sed -i s/[][]//g output")
}
# Single quote "'" (sometimes used as "prime" as in "3'", "5'")
if (remove == "ac" || remove == "all") {
   system("sed -i s/\\'//g output")
   # system("sed -i s/\\'//g")
}

# Gzip output if input was gzipped
if (gzip) {
   system("gzip output")
   system("ls -l")
}


# Output display name
# Determine base name
inputnames <- read_input_definitions()
name <- paste(inputnames$input)

# Make a matrix of output names
outputnames <- matrix(NA, nrow = 2, ncol = 2)
outputnames[1, ] <- c("output", name)
outputnames[2, ] <- c("output.gz", name)

# Write output definitions file
write_output_definitions(outputnames)
