# TOOL classification-knn.R: "KNN classification" (K-nearest neighbor classification. If you have a separate test data set, you can validate the prediction with it by setting the validation type to predict. This function does not perform any gene selection, and the analysis is run for the selected data set.)
# INPUT normalized.tsv: normalized.tsv TYPE GENE_EXPRS
# INPUT META phenodata.tsv: phenodata.tsv TYPE GENERIC
# OUTPUT knn-cross-validation.tsv: knn-cross-validation.tsv
# PARAMETER groups: "Groups" TYPE METACOLUMN_SEL DEFAULT group (Phenodata column describing the groups to test)
# PARAMETER training: "Training" TYPE METACOLUMN_SEL DEFAULT EMPTY (Phenodata column describing training chips)
# PARAMETER number.of.nearest.neighbors: "Number of nearest neighbors" TYPE INTEGER FROM 1 TO 1000 DEFAULT 2 (Number of nearest neighbors)
# PARAMETER number.of.votes: "Number of votes" TYPE INTEGER FROM 0 TO 1000 DEFAULT 2 (Number of votes needed to get a definite answer)
# PARAMETER validation.type: "Validation type" TYPE [crossvalidate: "cross-validate", predict: predict] DEFAULT crossvalidate (Type of analysis)

# JTT 26.6.2006: KNN classification created
# MK 02.07.2013: Bugs fixed, added group and traning column selection to the header section

# Parameter settings (default) for testing purposes
# number.of.nearest.neighbors<-2
# number.of.votes<-2
# validation.type<-"crossvalidate"

# Renaming the variables
knn.type <- validation.type
k.no <- number.of.nearest.neighbors
k.vote <- number.of.votes

# Load the libraries
library(class)

# Loads the data
file <- c("normalized.tsv")
dat <- read.table(file, sep = "\t", header = T, row.names = 1)

# Reads the phenodata table
phenodata <- read.table("phenodata.tsv", header = T, sep = "\t")

# Checks whether the training variable is present if phenodata
# If training part of phenodata has not been filled, but the column (header) is present,
# all the chips belong to the training set
if (training != "EMPTY") {
    if (training %in% colnames(phenodata)) {
        tr <- phenodata[, training]

        if (length(unique(tr)) != 2) {
            stop("CHIPSTER-NOTE: The training column in your phenodata must contain exactly two groups!")
        }

        if (tr[1] == " " | tr[1] == "" | any(is.na(tr)) == T) {
            tr <- rep(1, nrow(phenodata))
        }
    }
}

# Separates expression values and flags
calls <- dat[, grep("flag", names(dat))]
dat2 <- dat[, grep("chip", names(dat))]

# Are the parameter values sensical?
if (k.no > length(dat2)) {
    stop("CHIPSTER-NOTE: The number of neighbors is larger than the number of chips!")
}
if (k.vote > k.no) {
    stop("CHIPSTER-NOTE: The number of votes needed to give a definitive answer is larger than the number of neighbors!")
}

if (validation.type == "predict") {
    # Which parts of the data are training and test sets?
    dat3 <- split(as.data.frame(t(dat2)), tr)
    train <- dat3$"1"
    test <- dat3$"2"
    # Which part of the phenodata are training and test set
    phenodata_training <- phenodata[tr == unique(tr)[1], ]
    if (knn.type == "predict") {
        phenodata_test <- phenodata[tr == unique(tr)[2], ]
    }
} else {
    train <- t(dat2)
}

# Defines the true classification of the training set
if (validation.type == "predict") {
    cl <- phenodata_training[, groups]
} else {
    cl <- phenodata[, groups]
}

# Runs the KNN analysis and reports the results
if (knn.type == "crossvalidate") {
    knn.cross <- knn.cv(train = train, cl = cl, k = k.no, l = k.vote)
    # Writes a table of known versus predicted classes
    write.table(data.frame(sample = rownames(train), known.classes = cl, prediction = knn.cross), file = "knn-cross-validation.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
}

if (knn.type == "predict") {
    cl_test <- phenodata_test[, groups]
    knn.predict <- knn(train = train, test = test, cl = cl, k = k.no, l = k.vote)
    # Writes a table of known versus predicted classes
    write.table(data.frame(sample = rownames(test), known.classes = cl_test, prediction = knn.predict), file = "knn-cross-validation.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
}
