# TOOL share_via_object_storage.R: "Share a file" (Uploads selected file to Object Storage service of Chipster and gennerates a public link for the file)
# INPUT file: File TYPE GENERIC 
# OUTPUT link.html 
# OUTPUT OPTIONAL file.md5
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL add_md5sum: "Include md5 sum" TYPE [yes: Yes, no: No] DEFAULT no (Add md5 checksum file to the shared link)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the Mimimap2 mapping process.)

source(file.path(chipster.common.path, "tool-utils.R"))

#Check s3cmd and configuration
s3cmd_conf <- as.integer(system("s3cmd --dump-config 2> /dev/null | wc -l", intern = TRUE))
if ( s3cmd_conf < 1 ){
	stop("CHIPSTER-NOTE: Share a file tool has not been configured for this Chipster server.")
}

# 
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")


key <- stringi::stri_rand_strings(1, 12, pattern = "[a-z0-9]")
now <- as.integer(system("date +%s", intern = TRUE))
then <- now + 604800

date_command <- paste('date --date="@', then, '" +%d-%m-%y 2>&1', sep="")
#echo.command <- paste("echo '", date_command, "' >> log.txt")
#system(echo.command)

exp_date <- system(date_command, intern = TRUE)
#"system(paste("echo now: ",now ," then: ", then, " exp_date:", exp_date, " >> log.txt"))

cat(exp_date, file="log.txt", append=TRUE, sep="\n")

bname1 <- paste(key, exp_date, sep="-")
bname <- paste("s3://", bname1, sep="")
s3command <- paste("s3cmd mb -P ", bname, " 2>&1  >> log.txt" )
#echo.command <- paste ("echo ' ", s3command, "' >> log.txt" )
#system(echo.command)
system(s3command)

fname <- input.names[1,2]
system(paste("cp file ", fname ))

system("ls -l >> log.txt")

s3command <- paste("s3cmd put -P ", fname, bname, " >> log.txt" )
system(s3command)

if ( add_md5sum == "yes" ){
	md5file <- paste(fname, ".md5", sep="")
	system(paste("md5sum ", fname, " > ", md5file ))
	s3command <- paste("s3cmd put -P ", md5file, bname, " >> log.txt" )
	system(s3command)
}

cat("<html>",file="link.html",sep="\n")
cat("Temporary public link for file ",file="link.html", append=TRUE)
ns <- paste(fname, ":")
cat( ns , file="link.html", append=TRUE, sep="\n")
link <- paste('<a href="https://',bname1 ,'.object.pouta.csc.fi/', fname,'">https://', bname1, '.object.pouta.csc.fi/', fname, '</a>', sep="")
cat(link, file="link.html", append=TRUE, sep="\n")
if ( add_md5sum == "yes" ){
	link <- paste('<a href="https://',bname1 ,'.object.pouta.csc.fi/', md5file,'">https://', bname1, '.object.pouta.csc.fi/', md5file, '</a>', sep="")
	cat(link, file="link.html", append=TRUE, sep="\n")
}

cat("</html>",file="link.html", append=TRUE, sep="\n")

if ( save_log == "no") {
	system ("rm -f log.txt")
}










