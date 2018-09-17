# TOOL encrypt.R: "Encrypt a file" (Calculates GPG encrypted version of a file or opens encrypition)
# INPUT file: File TYPE GENERIC 
# OUTPUT OPTIONAL file.gpg 
# OUTPUT OPTIONAL file_decrypted 
# OUTPUT OPTIONAL keyfile 
# PARAMETER OPTIONAL task: "Select task" TYPE [encrypt: "Encrypt a file", decrypt: "Decrypt a file"] DEFAULT encrypt (Choose if the input file will be encrypted or decrypted)
# PARAMETER OPTIONAL pwlen: "Password legnth" TYPE INTEGER DEFAULT 16 (Length of the encrytion password to be generated. This fied is ignored in the case of decryption.)
# PARAMETER OPTIONAL dkey: "Decryption password" TYPE STRING (The password string that is used to open the encryption. This fiels is ignored in ecryption.)

source(file.path(chipster.common.path, "tool-utils.R"))

# We need to change the file name before compressing so it is preserved when uncompressing
input.names <- read.table("chipster-inputs.tsv", header=F, sep="\t")


if( task == "encrypt"){
    #rename input file
	key <- stringi::stri_rand_strings(1, pwlen)
	eccommand <- paste( 'gpg --yes --batch --passphrase="', key,'" -c file >> gpg.log', sep="" )
	system(eccommand)
	echo.command <- paste(" echo '", key, "' > keyfile", sep ="")
	system(echo.command)
	system ("ls -l >> gpg.log" )
	filename <- paste(input.names[1,2], ".gpg", sep = "")
	keyname <- paste(input.names[1,2], ".key.txt", sep = "")
	# And finally change Chipster display name to match original file name
	outputnames <- matrix(NA, nrow=2, ncol=2)
	outputnames[1,] <- c("file.gpg", filename)
	outputnames[2,] <- c("keyfile", keyname)
	# Write output definitions file
	write_output_definitions(outputnames)
	
}
if ( task == "decrypt") {
	# Gzip
	eccommand <- paste ( 'gpg --yes --batch --passphrase="', dkey,'" -d file > file_decrypted 2> dclog', sep="" )
	system(eccommand)
	dc.check <- system("grep -c 'failed:' dclog", intern = TRUE )
	if ( dc.check != 0 ){
		stop(paste ("CHIPSTER-NOTE: ERROR: The key", dkey, " does not match with the file to be decrypted." ) )		
	}
	
}


