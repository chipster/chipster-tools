# TOOL single-cell-bam-tag-histogram.R: "BAM tag histogram" (Extracts the reads per cell barcode and draw a cumulative distribution plot. The plot helps to select how many cells to include in the DGE matrix.)
# INPUT input.bam: "prepared BAM" TYPE GENERIC
# OUTPUT OPTIONAL cell_readcounts.txt.gz
# OUTPUT OPTIONAL tag_histogram.pdf
# OUTPUT OPTIONAL inflectionPoint.txt
# PARAMETER OPTIONAL x_axis_max: "Max to x axis" TYPE INTEGER FROM 0 TO 500000 DEFAULT 1500 (Upper limit for x-axis in the histogram. If you cannot see the knee in the curve, try tuning this parameter.) 


# OUTPUT OPTIONAL CumCellReadCount.pdf

# ML 12.10.2016 created

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command start
command <- paste(path.dropseq, "/BAMTagHistogram I=input.bam O=cell_readcounts.txt.gz TAG=XC", sep="")

# run the tool
command <- paste(command, " 2>> log.txt")

system(command)

#a=read.table("100cells_numReads_perCell_XC_mq_10.txt.gz", header=F, stringsAsFactors=F) x=cumsum(a$V1)
a=read.table("cell_readcounts.txt.gz", header=F, stringsAsFactors=F) 
x=cumsum(a$V1)
x=x/max(x)
#pdf(file="tag_histogram.pdf")
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,x_axis_max))
#dev.off()



# Using Dawit's script:
library("inflection")

# read in the cellReadCounts File
cellReadCountsFile <- "cell_readcounts.txt.gz" # added, by Maria
a=read.table(cellReadCountsFile, header=F, stringsAsFactors=F)
x=cumsum(a$V1)
x=x/max(x)

# inflection point finding using extreme distance estimator(ede) from inflection package
infl = ede(1:length(x),x,0)
print(infl[1,2])


# Plot Cumulative Read proportion of cells 
# pdf(file=paste(outputDir,"/CumCellReadCount.pdf",sep=""))
pdf(file="tag_histogram")
cellIds = 1:length(x)
plot(cellIds, x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
		ylab="cumulative fraction of reads", xlim=c(1,length(x)))

abline(v=infl[1,2],col="red")
abline(h=x[infl[1,2]],col="red")

axis(1, at=infl[1,2],labels=infl[1,2],col="red")
axis(2, at=x[infl[1,2]],labels=round(x[infl[1,2]],2),col="red",las=2)

dev.off()

#write inflectionPoint.txt 
write.table(infl[1,2],"inflectionPoint.txt",col.names = c("InflectionPoint"),row.names=F)



# stop(paste('CHIPSTER-NOTE: ', command))

#EOF