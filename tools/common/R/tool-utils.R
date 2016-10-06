read_input_definitions <- function(){
	# Read in the data
	list <- scan("chipster-inputs.tsv", what="", sep="\n", comment.char="#")
	# Separate elements by one or more whitepace
	inputdef <- strsplit(list, "[[:space:]]+")
	# Extract the first vector element and set it as the list element name
	names(inputdef) <- sapply(inputdef, function(list) list[[1]])
	# Remove the first vector element from each list element
	inputdef <- lapply(inputdef, function(list) list[-1])
	
	return(inputdef)
}							

write_output_definitions <- function(output_names){
	write.table(output_names, file = "chipster-outputs.tsv", row.names = F, col.names= F, quote = F, sep = "\t")
}

# Removes user-specified postfix if present
#
remove_postfix <- function(name, postfix){
	pf <- paste(postfix, "$", sep="")
	if (grepl(pf, name)){
		basename <- substr(name, 1, (nchar(name) - nchar(postfix)))
		return(basename)
	}else{
		return(name)
	}
}

# Removes extension, i.e. everything after the last dot (including the dot)
#
remove_extension <- function(name){
	return(sub("\\.[^.]*$", "", name))
}

# Strips common file extensions from a file name
#
strip_name <- function(name){
	known_postfixes <- c(".gz", ".bam", ".sam", ".fa", ".fasta", ".fq", ".fastq", ".gtf", ".tsv", ".txt", "_trimmed", "_filtered")
	newname <- name
	while (TRUE){
		for (i in known_postfixes){
			newname <- remove_postfix(newname, i)
		}
		if (nchar(newname) == nchar(name)){
			break
		}
		name <- newname
	}
	return(name)
	
}

# If the names look like typical paired-end names: *_1, *_2, remove the ending and return the name.
# If not, return first name as-is.
#
paired_name <- function(name1, name2){
	if (grepl("_1$", name1) && grepl("_2$", name2)){
		return(remove_postfix(name1, "_1"))
	}
	if (grepl("_2$", name1) && grepl("_1$", name2)){
		return(remove_postfix(name1, "_2"))
	}
	return(name1)
}

# Takes as an argument a file name of a list of input names. It Goes through the list and compares 
# it to chipster-inputs.tsv. Returns a vector of internal input names or gives an error message, 
# if a listed input file is missing.
#
make_input_list <-function(listfile){
	
	# read list file
	name.list <- scan(listfile, what="", sep="\n")
	
	# read input names
	input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")
	
	# Check for duplicated etries
	if (anyDuplicated(name.list)){
		message <- paste("Input file list contains one or more duplicated entries:\n", name.list[duplicated(name.list)], "\nFilenames must be unique so they can be correctly assigned.")
		stop(paste('CHIPSTER-NOTE: ', message ))
	}
	
	# Check that inputs exist
	sdf <- setdiff(name.list, input.names[,2])
	if (identical(sdf, character(0))){
		input.list <- vector(mode="character", length=0)
		for (i in 1:length(name.list)) {
			input.list <- c(input.list, paste(input.names[grep(paste("^", name.list[i], "$", sep =""), input.names[,2]), 1]))
		}
		
	}else{
		message <- paste("Input file list includes one or more files that has not been selected:", sdf)
		stop(paste('CHIPSTER-NOTE: ', message))
		
	}
	return(input.list)
}