merge_tables <- function(input.names, include.everything){
	# Read in and clean the first table
	tablename <- paste(input.names[1,1])
	merged <- read.table(tablename, sep="\t", header=T, row.names=1)
	
	for (i in 1:ncol(merged)) {
		merged[,i] <- gsub("\t+", " ", merged[,i], perl=T)
		merged[,i] <- gsub("\n+", " ", merged[,i], perl=T)
		merged[,i] <- gsub(" +", " ", merged[,i], perl=T)
		
		merged[,i] <- gsub("\"+", ",", merged[,i], perl=T)
		merged[,i] <- gsub("\'+", ",", merged[,i], perl=T)
	}
	
	# Cycle through the rest of the tables
	for (j in 2:nrow(input.names)) {
		# Read the next table
		tablename <- paste(input.names[j,1])
		table<-read.table(tablename, sep="\t", header=T, row.names=1)
		# Clean the table
		for (i in 1:ncol(table)) {
			table[,i] <- gsub("\t+", " ", table[,i], perl=T)
			table[,i] <- gsub("\n+", " ", table[,i], perl=T)
			table[,i] <- gsub(" +", " ", table[,i], perl=T)
			
			table[,i] <- gsub("\"+", ",", table[,i], perl=T)
			table[,i] <- gsub("\'+", ",", table[,i], perl=T)
		}
		# Combine tables using row names
		include=FALSE
		if( include.everything == "yes" ) include=TRUE
		merged <- merge(merged, table, by.x="row.names", by.y="row.names", all=include)
		row.names(merged)<-merged$Row.names
		merged<-merged[-1]
		
	}
	
	# Writes out the combined table
	# write.table(merged, "combined.tsv", sep="\t", row.names=T, col.names=T, quote=F)
	return(merged)
}	