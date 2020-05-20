read_input_definitions <- function() {
  # Read in the data
  list <- scan("chipster-inputs.tsv",what = "",sep = "\n",comment.char = "#")
  # Separate elements by one or more whitepace
  inputdef <- strsplit(list,"[\t]+")
  # Extract the first vector element and set it as the list element name
  names(inputdef) <- sapply(inputdef,function(list) list[[1]])
  # Remove the first vector element from each list element
  inputdef <- lapply(inputdef,function(list) list[-1])

  return(inputdef)
}

write_output_definitions <- function(output_names) {
  write.table(output_names,file = "chipster-outputs.tsv",row.names = FALSE,col.names = FALSE,quote = FALSE,sep = "\t")
}

# Removes user-specified postfix if present
#
remove_postfix <- function(name,postfix) {
  pf <- paste(postfix,"$",sep = "")
  if (grepl(pf,name)) {
    basename <- substr(name,1,(nchar(name) - nchar(postfix)))
    return(basename)
  } else {
    return(name)
  }
}

# Removes extension, i.e. everything after the last dot (including the dot)
#
remove_extension <- function(name) {
  return(sub("\\.[^.]*$","",name))
}

# Strips common file extensions from a file name
#
strip_name <- function(name) {
  known_postfixes <- c(".gz",".bam",".sam",".fa",".fasta",".fq",".fastq",".gtf",".tsv",".txt","vcf","_trimmed","_filtered")
  newname <- name
  while (TRUE) {
    for (i in known_postfixes) {
      newname <- remove_postfix(newname,i)
    }
    if (nchar(newname) == nchar(name)) {
      break
    }
    name <- newname
  }
  return(name)

}

# If the names look like typical paired-end names: *_1, *_2, remove the ending and return the name.
# If not, return first name as-is.
#
paired_name <- function(name1,name2) {
  if (grepl("_1$",name1) && grepl("_2$",name2)) {
    return(remove_postfix(name1,"_1"))
  }
  if (grepl("_2$",name1) && grepl("_1$",name2)) {
    return(remove_postfix(name1,"_2"))
  }
  return(name1)
}

# Takes as an argument a file name of a list of input names. It Goes through the list and compares 
# it to chipster-inputs.tsv. Returns a vector of internal input names or gives an error message, 
# if a listed input file is missing.
#
make_input_list <- function(listfile) {

  # read list file
  name.list <- scan(listfile,what = "",sep = "\n")

  # read input names
  input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")

  # Check for duplicated etries
  if (anyDuplicated(name.list)) {
    message <- paste("Input file list contains one or more duplicated entries:\n",name.list[duplicated(name.list)],"\nFilenames must be unique so they can be correctly assigned.")
    stop(paste('CHIPSTER-NOTE: ',message))
  }

  # Check that inputs exist
  sdf <- setdiff(name.list,input.names[,2])
  if (identical(sdf,character(0))) {
    input.list <- vector(mode = "character",length = 0)
    for (i in 1:length(name.list)) {
      input.list <- c(input.list,paste(input.names[grep(paste("^",name.list[i],"$",sep = ""),input.names[,2]),1]))
    }

  } else {
    message <- paste("Input file list includes one or more files that has not been selected:",sdf)
    stop(paste('CHIPSTER-NOTE: ',message))

  }
  return(input.list)
}

# fileOk and fileNotOk ar functions to check that file exists and optionally that it is larger than the provided threshold.
# fileOk("myfile.txt") : returns TRUE if file myfile.txt exists
# fileOk("myfile.txt", 100) : returns TRUE if file myfile.txt exists and is min. 100 bytes
# fileOk("myfile.txt",,2) : returnd TRUE if file myfile.txt exists and is min. 2 lines
# The two functions work the same exept fileOk returns TRUE if the file meets the criteria, and fileNotOk returns TRUE
# if it does not. The idea is to make if statements more readable.

fileOk <- function(filename,minsize,minlines) {
  # Check that file exists
  if (!(file.exists(filename))) {
    return(FALSE)
  }
  # Check minimum size if provided
  if (!(missing(minsize))) {
    if (file.info(filename)$size < minsize) {
      return(FALSE)
    }
  }
  # Check line number if provided
  if (!(missing(minlines))) {
    linenumber <- as.integer(system(paste("wc -l <",filename),intern = TRUE))
    if (linenumber < minlines) {
      return(FALSE)
    }
  }
  # If all checks pass, return TRUE
  return(TRUE)
}

fileNotOk <- function(filename,minsize,minlines) {
  # Check if file exists
  if (!(file.exists(filename))) {
    return(TRUE)
  }
  # Check minimum size if provided
  if (!(missing(minsize))) {
    if (file.info(filename)$size < minsize) {
      return(TRUE)
    }
  }
  # Checkline number if provided
  if (!(missing(minlines))) {
    linenumber <- as.integer(system(paste("wc -l <",filename),intern = TRUE))
    if (linenumber < minlines) {
      return(TRUE)
    }
  }
  # If all checks pass, return FALSE
  return(FALSE)
}

fileCheck <- function(filename,minsize,minlines) {
  # Check if file exists
  if (!(file.exists(filename))) {
    stop(paste('CHIPSTER-NOTE: ',"Required file",filename,"does not exist."))
  }
  # Check minimum size if provided
  if (!(missing(minsize))) {
    filesize <- file.info(filename)$size
    if (filesize < minsize) {
      stop(paste('CHIPSTER-NOTE: ',"Required file",filename,"exists, but is smaller (",filesize,") than required size (",minsize,")."))
    }
  }
  # Checkline number if provided
  if (!(missing(minlines))) {
    linenumber <- as.integer(system(paste("wc -l <",filename),intern = TRUE))
    if (linenumber < minlines) {
      stop(paste('CHIPSTER-NOTE: ',"Required file",filename,"exists, but has fewer lines (",linenumber,") than required (",minlines,")."))
    }
  }
}

# Wrapper for system2() command. 
# Optionally captures stderr.
# Checks for exit status and gives an error message if status != 0
# Allows setting environment variables give as character vector in style "VARIABLE=value"
#
runExternal <- function(command,env = NULL,capture = TRUE,checkexit = TRUE) {

  # Split command to words
  wcom <- strsplit(command," ")[[1]]

  # if environment varaiables are set, the command has to be started with "bash -c" and encased in single quotes.
  if (!is.null(env)) {
    wcom <- c("bash","-c","\'",wcom,"\'")
  }

  # Add timestamp and command line to stderr.log
  system("echo \"# Chipster log\" >> stderr.log")
  system("echo \"# \"\`date\` >> stderr.log")
  write(paste(wcom,collapse = " "),file = "stderr.log",append = TRUE)

  # Capture of stderr is optional
  if (capture) {
    # Run command, capture stderr
    exitcode <- system2(wcom[1],wcom[2:length(wcom)],stderr = "stderr.tmp",env = env)

    # Append error messages to the log file. stderr.tmp is overwritten each time runExternal is called.
    system("echo >> stderr.log")
    system("cat stderr.tmp >> stderr.log")
    system("echo  >> stderr.log")
  } else {
    # Run command without capturing stderr
    exitcode <- system2(wcom[1],wcom[2:length(wcom)],env = env)
  }
  # Show error message if command fails
  if (checkexit) {
    if (exitcode != 0) {

      error <- scan("stderr.tmp",what = " ",sep = "\n")

      msg <- paste("External program did not complete succesfully.\n")
      msg <- paste(msg,"Please check your parameters, input files and input file assignment.\n\n")
      msg <- paste(msg,"Failed command:\n")
      msg <- paste(msg,command,"\n\n")
      msg <- paste(msg,"Error:\n")
      msg <- paste(msg,paste(error,collapse = "\n"))

      stop(paste('CHIPSTER-NOTE: ',msg))
    }
  }
}

# Changes the file names in a text file to display names according to chipster-inputs.tsv
#
displayNamesToFile <- function(input.file) {

  # Read input names
  input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")
  # Go through input names and change names
  for (i in 1:nrow(input.names)) {
    sed.command <- paste("s/",input.names[i,1],"/",input.names[i,2],"/",sep = "")
    system(paste("sed -i",sed.command,input.file))
  }
}

# Formats and prints out the command to stdout. Input names are substituted with
# display names to aid readability and help spot input assignment errors.
#
documentCommand <- function(command.string) {
  # Substitute input names
  input.names <- read.table("chipster-inputs.tsv",header = FALSE,sep = "\t")
  for (i in 1:nrow(input.names)) {
    command.string <- gsub(input.names[i,1],input.names[i,2],command.string)
  }
  cat("##","COMMAND:",command.string,"\n")
}

# Prints out version information.
#
documentVersion <- function(application,version.string) {
  cat("##","VERSION:",application,"\n")
  cat("##",version.string,"\n")
}