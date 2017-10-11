# TOOL single-cell-bam-tag-histogram.R: "Estimate number of usable cells" (Extracts the reads per cell barcode and draw a cumulative distribution plot. The plot helps to select how many cells to include in the DGE matrix.)
# INPUT input.bam: "prepared BAM" TYPE GENERIC
# OUTPUT OPTIONAL cell_readcounts.txt.gz
# OUTPUT OPTIONAL inflectionPoint.pdf
# PARAMETER OPTIONAL select_cells: "Select only a number of top cell barcodes for plotting" TYPE [TRUE, FALSE] DEFAULT FALSE (Switch to TRUE if you want to select the number of the top cell barcodes, and set the number in the parameter below. You might want to do this if you suspect that there are lots of nonsense cell barcodes in your data.)
# PARAMETER OPTIONAL number_cells: "Number of cell barcodes to use" TYPE INTEGER FROM 0 TO 500000 DEFAULT 0 (How many cell barcodes to use when drawing the knee plot. If you suspect that there are lots of nonsense cell barcodes in your data, you can try tuning this parameter to a number closer to expected number of actual beads in your data. NOTE that if you did not select TRUE in the parameter above, this parameter does nothing.) 


# OUTPUT OPTIONAL CumCellReadCount.pdf
# OUTPUT OPTIONAL inflectionPoint.txt

# ML 12.10.2016 created
# ML 19.7.2017 Change tool name and add the inflection point information to the pdf

path.dropseq <- c(file.path(chipster.tools.path, "drop-seq_tools"))

# command start
command <- paste(path.dropseq, "/BAMTagHistogram I=input.bam O=cell_readcounts.txt.gz TAG=XC", sep="")

# run the tool
command <- paste(command, " 2>> log.txt")

system(command)

#a=read.table("100cells_numReads_perCell_XC_mq_10.txt.gz", header=F, stringsAsFactors=F) x=cumsum(a$V1)
#a=read.table("cell_readcounts.txt.gz", header=F, stringsAsFactors=F) 
#x=cumsum(a$V1)
#x=x/max(x)
#pdf(file="tag_histogram.pdf")
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,x_axis_max))
#dev.off()


# Using Dawit's script:
library("inflection")

# read in the cellReadCounts File
cellReadCountsFile <- "cell_readcounts.txt.gz" # added, by Maria
a=read.table(cellReadCountsFile, header=F, stringsAsFactors=F)
# total number of barcodes:
num_barcodes <- dim(a)[1]

if (select_cells == TRUE) { 
	if (number_cells < 10) { 
		stop(paste('CHIPSTER-NOTE: ', "You wanted to select the number of cells to use -please set a reasonable (bigger) number to the parameter."))	
	}else if (number_cells > num_barcodes) { 
		stop(paste('CHIPSTER-NOTE: ', "You wanted to select the number of cells to use, but the number you chose exceeded the number of barcodes in the dataset."))	
	}else {
		a2 <- a[1:number_cells,]
		x=cumsum(a2$V1)
		x=x/max(x)
		x_axis_max <- number_cells
	}
} else{ 
	x=cumsum(a$V1)
	x=x/max(x)
	x_axis_max <- num_barcodes
}

# inflection point finding using extreme distance estimator(ede) from inflection package
infl = ede(1:length(x),x,0)
# print(infl[1,2])

# percentage of reads covered by cells in the inflection point:
per_reads <- round(x[infl[1,2]], digits=2)*100

pdf(file="inflectionPoint.pdf")
cellIds = 1:length(x)
plot(cellIds, x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]",
		ylab="cumulative fraction of reads", xlim=c(1,length(x)), 
		main=paste("Inflection point: ", infl[1,2], ", covers ", per_reads, " % of reads. \n Total number of cell barcodes detected: ", 
				num_barcodes, "\n  Number of barcodes used to compute the inflection point: ",x_axis_max  ))
#     sub=paste("Total number of cell barcodes detected: ", num_barcodes), cex.sub=1.5 )

abline(v=infl[1,2],col="red")
abline(h=x[infl[1,2]],col="red")

axis(1, at=infl[1,2],labels=infl[1,2],col="red")
axis(2, at=x[infl[1,2]],labels=round(x[infl[1,2]],2),col="red",las=2)

dev.off()

## write inflectionPoint.txt 
# write.table(infl[1,2],"inflectionPoint.txt",col.names = c("InflectionPoint"),row.names=F)



# stop(paste('CHIPSTER-NOTE: ', command))

#EOF