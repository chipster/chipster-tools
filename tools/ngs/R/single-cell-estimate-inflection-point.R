# By Dawit Yohannes (Univ of Helsinki)
# A modified version of Macosco et al R script 
# to enable estimation of the "knee point" numerically (using R package inflection)

# Input: two Arguments:
# - cell read counts (generated using BAMTagHistogram of Drop-seq_tools), and 
# - output directory where the result files are put
# Output: two files:
# - CumCellReadCount.pdf (the cumulative cell read count plot), and 
# -inflectionPoint.txt is 1row x 1col table whose one element holds the estimated inflection point (also indicated on the plot)

# Example: Rscript /wrk/dawit/dropSeqProjects/estimateInflectionPointNumCells.R /wrk/dawit/dropSeqProjects/Thymus/cell_readcounts.txt.gz /wrk/dawit/dropSeqProjects/inflectionPointTesting/



#-------------------------------------------------------------------------------
# read file name arguments : 
args <- commandArgs(trailingOnly = TRUE)
cellReadCountsFile=args[1]
outputDir=args[2]


## Install inflection package if it is not available already
#loadPacks <- function(package.list = c("inflection")){
#	new.packages <-package.list[!(package.list %in% installed.packages()[,"Package"])]
#	if(length(new.packages)) install.packages(new.packages)
#	lapply(eval(package.list), require, character.only=TRUE);
#}
#loadPacks();

# INSTEAD: install.packages("inflection")  <--- Taavi / Petri!


# read in the cellReadCounts File
a=read.table(cellReadCountsFile, header=F, stringsAsFactors=F)
x=cumsum(a$V1)
x=x/max(x)

# inflection point finding using extreme distance estimator(ede) from inflection package
infl = ede(1:length(x),x,0)
print(infl[1,2])


# Plot Cumulative Read proportion of cells 
pdf(file=paste(outputDir,"/CumCellReadCount.pdf",sep=""))
cellIds = 1:length(x)
plot(cellIds, x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
		ylab="cumulative fraction of reads", xlim=c(1,length(x)))

abline(v=infl[1,2],col="red")
abline(h=x[infl[1,2]],col="red")

axis(1, at=infl[1,2],labels=infl[1,2],col="red")
axis(2, at=x[infl[1,2]],labels=round(x[infl[1,2]],2),col="red",las=2)

dev.off()


#write inflectionPoint.txt 
write.table(infl[1,2],paste(outputDir,"/inflectionPoint.txt",sep=""),col.names = c("InflectionPoint"),row.names=F)
