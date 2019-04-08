# TOOL share_via_object_storage.R: "Share a file" (Uploads the selected file to Object Storage service of Chipster and generates a public link for the file.)
# INPUT file: File TYPE GENERIC 
# OUTPUT link.html 
# OUTPUT OPTIONAL file.md5
# OUTPUT OPTIONAL log.txt
# PARAMETER OPTIONAL add_md5sum: "Include md5 sum" TYPE [yes: Yes, no: No] DEFAULT no (Add md5 checksum file to the shared link)
# PARAMETER OPTIONAL save_log: "Collect a log file" TYPE [yes: Yes, no: No] DEFAULT no (Collect a log file about the Mimimap2 mapping process.)

source(file.path(chipster.common.path, "tool-utils.R"))

s3.conf.file<- file.path(chipster.common.path, "../../admin/shell/s3.conf")

s3cmd.binary <- file.path(chipster.tools.path, "Python-2.7.12", "bin", "s3cmd")
s3cmd.binary <- paste(s3cmd.binary, "-c", s3.conf.file)


#Check s3cmd and configuration
s3check.command <- paste(s3cmd.binary, "--dump-config 2> /dev/null | wc -l")
s3cmd_conf <- as.integer(system(s3check.command, intern = TRUE))
if ( s3cmd_conf < 1 ){
	stop(paste("CHIPSTER-NOTE: Share a file tool has not been configured for this Chipster server.", s3check.command ))
}

# 
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")

# storage time in seconds
validfor <- 432000

key <- stringi::stri_rand_strings(1, 12, pattern = "[a-z0-9]")
now <- as.integer(system("date +%s", intern = TRUE))
then <- now + validfor

date_command <- paste('date --date="@', then, '" +%d-%m-%y 2>&1', sep="")
#echo.command <- paste("echo '", date_command, "' >> log.txt")
#system(echo.command)

exp_date <- system(date_command, intern = TRUE)
#"system(paste("echo now: ",now ," then: ", then, " exp_date:", exp_date, " >> log.txt"))

cat(exp_date, file="log.txt", append=TRUE, sep="\n")

bname1 <- paste("chitemp", then, key, sep="-")
bname <- paste("s3://", bname1, sep="")
s3command <- paste( s3cmd.binary ,"mb -P ", bname, " 2>&1  >> log.txt" )
cat(s3command, file="log.txt", append=TRUE, sep="\n")
#echo.command <- paste ("echo ' ", s3command, "' >> log.txt" )
#system(echo.command)
system(s3command)

fname <- input.names[1,2]
system(paste("cp file ", fname ))

system("ls -l >> log.txt")

s3command <- paste(s3cmd.binary, "put -P ", fname, bname, " >> log.txt" )
system(s3command)
cat(s3command, file="log.txt", append=TRUE, sep="\n")

if ( add_md5sum == "yes" ){
	md5file <- paste(fname, ".md5", sep="")
	system(paste("md5sum ", fname, " > ", md5file ))
	s3command <- paste(s3cmd.binary, "put -P ", md5file, bname, " >> log.txt" )
	system(s3command)
}

cat("<html>",file="link.html",sep="\n")
cat("Temporary public link for file ",file="link.html", append=TRUE)
ns <- paste(fname, ":")
cat( ns , file="link.html", append=TRUE, sep="\n")
cat("</br>", file="link.html", append=TRUE, sep="\n")
link <- paste('<a href="https://',bname1 ,'.object.pouta.csc.fi/', fname,'">https://', bname1, '.object.pouta.csc.fi/', fname, '</a>', sep="")
cat("<li>", file="link.html", append=TRUE, sep="\n")
cat(link, file="link.html", append=TRUE, sep="\n")
cat("</li>", file="link.html", append=TRUE, sep="\n")
if ( add_md5sum == "yes" ){
	cat("<li>", file="link.html", append=TRUE, sep="\n")
	link <- paste('<a href="https://',bname1 ,'.object.pouta.csc.fi/', md5file,'">https://', bname1, '.object.pouta.csc.fi/', md5file, '</a>', sep="")
	cat(link, file="link.html", append=TRUE, sep="\n")
	cat("</li>", file="link.html", append=TRUE, sep="\n")
	cat("Links are valid until:",file="link.html", append=TRUE, sep=" " )
	
} else {	
	cat("The link is valid until:",file="link.html", append=TRUE, sep=" " )
}
cat(exp_date,file="link.html", append=TRUE, sep="\n ")

cat("</html>",file="link.html", append=TRUE, sep="\n")


#Clean old stuff
s3command <- paste(s3cmd.binary, "ls | grep 's3://chitemp-' | awk '{print $3}' | awk -F '-' '{ if ( $2 < ", now, " ) print $0 }' > rb-list.tmp " )
system(s3command)
system("echo removing outdated buckets: >> log.txt")
system("cat rb-list.tmp >> log.txt")
forloop <- paste("for f in $(cat rb-list.tmp); do ", s3cmd.binary, " rb --recursive $f >> log.txt; done")
system(forloop)

if ( save_log == "no") {
	system ("rm -f log.txt")
}











