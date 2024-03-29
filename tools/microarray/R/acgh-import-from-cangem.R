# TOOL acgh-import-from-cangem.R: "Import from CanGEM" (Load a microarray data set from the CanGEM database, perform background-correction and normalization, and append chromosomal locations of the microarray probes.)
# OUTPUT normalized.tsv: normalized.tsv
# OUTPUT META phenodata.tsv: phenodata.tsv
# PARAMETER accession: Accession TYPE STRING DEFAULT CG- (Accession of either a data set, an experiment, a series, or single microarray results.)
# PARAMETER username: Username TYPE STRING DEFAULT empty (Username, in case the data is password-protected. WARNING: This will store your username password in the Chipster history files. To avoid this, use the session parameter.)
# PARAMETER password: Password TYPE STRING DEFAULT empty (Password, in case the data is password-protected. WARNING: This will store your username password in the Chipster history files. To avoid this, use the session parameter.)
# PARAMETER session: "Session ID" TYPE STRING DEFAULT empty (Session ID. To avoid saving your username password in Chipster history files, log in at http: www.cangem.org using a web browser, then copy&paste your session ID from the lower right corner of the CanGEM website. This will allow Chipster to access your password-protected data until you log out of the web site (or the session times out\).)
# PARAMETER agilent.filtering: "Agilent filtering" TYPE [yes: yes, no: no] DEFAULT yes (Whether to filter outliers from Agilent 2-color arrays. Will be ignored, if downloaded files are 1-color arrays, or not in Agilent file format. Check the help file for details about the filtering function.)
# PARAMETER background.treatment: "Background treatment" TYPE [none: none, subtract: subtract, normexp: normexp, rma: rma] DEFAULT normexp (Background treatment method. RMA is available only for one-color arrays.)
# PARAMETER background.offset: "Background offset" TYPE [0: 0, 50: 50] DEFAULT 50 (Background offset.)
# PARAMETER intra.array.normalization: "Intra array normalization" TYPE [none: none, median: median, loess: loess] DEFAULT loess (Intra-array normalization method for Agilent arrays. Will be ignored, if downloaded files are not in Agilent file format.)
# PARAMETER inter.array.normalization: "Inter array normalization" TYPE [none: none, quantile: quantile, scale: scale] DEFAULT none (Inter-array normalization method for Agilent arrays. Will be ignored, if downloaded files are not in Agilent file format.)
# PARAMETER affymetrix.normalization: "Affymetrix normalization" TYPE [gcrma: gcrma, rma: rma, mas5: mas5] DEFAULT gcrma (Normalization method for Affymetrix arrays. Will be ignored, if downloaded files are not in Affymetrix file format.)
# PARAMETER genome.build: "Genome build" TYPE [none: none, GRCh37: GRCh37, NCBI36: NCBI36, NCBI35: NCBI35, NCBI34: NCBI34] DEFAULT GRCh37 (The genome build to use. GRCh37 = hg19, NCBI36 = hg18, NCBI35 = hg17, NCBI34 = hg16.)

# Ilari Scheinin <firstname.lastname@gmail.com>
# 2014-03-22

source(file.path(chipster.common.path, "library-Chipster.R"))

# check for valid accession
accession <- toupper(accession)
if (length(grep("^CG-(SET|EXP|SER|RES|SAM|PLM)-[0-9]+$", accession)) == 0) {
  stop("CHIPSTER-NOTE: Not a valid accession: ", accession)
}

# construct the string used in authenticating
if (session != "empty" && session != "") {
  auth <- paste("&PHPSESSID=", session, sep = "")
} else if (username != "empty" && username != "" && password != "empty" && password != "") {
  auth <- paste("&username=", username, "&password=", password, sep = "")
} else {
  auth <- ""
}

# fetch list of arrays from CanGEM
cangem.samples <- read.table(paste("http://www.cangem.org/scripts/listhybs.php?accession=", accession, auth, sep = ""), header = TRUE, sep = "\t", quote = "", as.is = TRUE, comment.char = "")
cangem.samples$GUID <- NULL

# check that we did get some results
if (nrow(cangem.samples) == 0) {
  stop("CHIPSTER-NOTE: No array data found.")
}

# check for single platform, array type and file format
if (length(unique(cangem.samples$PlatformAccession)) > 1) {
  stop("CHIPSTER-NOTE: Multiple platforms used: ", paste(unique(cangem.samples$Platform), collapse = ", "))
}
if (length(unique(cangem.samples$Type)) > 1) {
  stop("CHIPSTER-NOTE: Multiple array types used: ", paste(unique(cangem.samples$Type), collapse = ", "))
}
if (length(unique(cangem.samples$Format)) > 1) {
  stop("CHIPSTER-NOTE: Multiple file formats used: ", paste(unique(cangem.samples$Format), collapse = ", "))
}

microarrays <- sprintf("microarray%.3i", 1:nrow(cangem.samples)) # or .cel/.tsv ???
chips <- paste("chip.", microarrays, sep = "")
arrays <- list()
probes.to.genes <- data.frame()
onecolor <- FALSE
if (cangem.samples$Format[1] == "agilent") {
  library(limma)
  # define filtering function
  if (agilent.filtering == "yes") {
    cangem.filter <- function(x) {
      # required fields and their values, field.name=required.value
      reqs <- c(
        ControlType = 0,
        gIsBGNonUnifOL = 0,
        gIsBGPopnOL = 0,
        gIsFeatNonUnifOL = 0,
        gIsFeatPopnOL = 0,
        gIsFound = 1,
        gIsPosAndSignif = 1,
        gIsSaturated = 0,
        gIsWellAboveBG = 1,
        gSurrogateUsed = 0,
        IsManualFlag = 0,
        rIsBGNonUnifOL = 0,
        rIsBGPopnOL = 0,
        rIsFeatNonUnifOL = 0,
        rIsFeatPopnOL = 0,
        rIsFound = 1,
        rIsPosAndSignif = 1,
        rIsSaturated = 0,
        rIsWellAboveBG = 1,
        rSurrogateUsed = 0
      )
      # start with a weight of 1 for every spot
      weights <- rep(1, length = nrow(x))
      # check which fields are present, and check if they have the required value
      for (req in names(reqs)) {
        if (req %in% colnames(x)) {
          weights <- weights & x[, req] == reqs[req]
        }
      }
      as.numeric(weights)
    }
  } else {
    cangem.filter <- function(x) {
      1
    }
  }
  # load arrays one by one to be able to handle with files with different number of rows
  for (i in 1:nrow(cangem.samples)) {
    prob <- TRUE
    if (!onecolor) {
      # first try to read the 'normal' columns (means), but if it fails, try medians
      try(
        {
          array.raw <- read.maimages(paste("http://www.cangem.org/download.php?hybridization=", cangem.samples[i, "Accession"], auth, sep = ""),
            source = cangem.samples[i, "Format"], wt.fun = cangem.filter
          )
          prob <- FALSE
        },
        silent = TRUE
      )
      if (prob) {
        try(
          {
            array.raw <- read.maimages(paste("http://www.cangem.org/download.php?hybridization=", cangem.samples[i, "Accession"], auth, sep = ""),
              source = cangem.samples[i, "Format"], columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal", R = "rMedianSignal", Rb = "rBGMedianSignal"), wt.fun = cangem.filter
            )
            prob <- FALSE
          },
          silent = TRUE
        )
      }
    }
    if (prob) {
      try(
        {
          array.raw <- read.maimages(paste("http://www.cangem.org/download.php?hybridization=", cangem.samples[i, "Accession"], auth, sep = ""),
            source = cangem.samples[i, "Format"], columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal", R = "gMedianSignal", Rb = "gBGMedianSignal"), wt.fun = cangem.filter
          )
          prob <- FALSE
          onecolor <- TRUE
        },
        silent = TRUE
      )
    }
    if (prob) {
      stop("CHIPSTER-NOTE: Could not read file ", cangem.samples[i, "Name"])
    }
    if (background.treatment == "rma") {
      array.bg <- backgroundCorrect(array.raw, method = "none")
    } else {
      array.bg <- backgroundCorrect(array.raw, method = background.treatment, normexp.method = "mle", offset = as.numeric(background.offset))
    }
    if (onecolor) {
      array <- as.vector(array.bg$R)
    } else {
      array.ma <- normalizeWithinArrays(array.bg, method = intra.array.normalization)
      array <- as.vector(array.ma$M)
      # take inverse for dye swaps
      if (cangem.samples[i, "SampleChannel"] == "Cy3") {
        array <- -array
      }
    }
    # names(array) <- array.bg$genes[,cangem.samples[i,'ProbeNames']]
    names(array) <- array.bg$genes[, "ProbeName"]
    # average replicate probes
    replicate.probes <- unique(names(array)[duplicated(names(array))])
    array.uniques <- array[!names(array) %in% replicate.probes]
    array.replicates <- array[names(array) %in% replicate.probes]
    array.replicates.avg <- aggregate(array.replicates, list(probe = names(array.replicates)), median, na.rm = TRUE)
    array.replicates <- array.replicates.avg$x
    names(array.replicates) <- array.replicates.avg$probe
    array <- c(array.uniques, array.replicates)
    arrays[[i]] <- array
    if (cangem.samples[i, "ProbeNames"] != "ProbeName") {
      probes.to.genes <- unique(rbind(probes.to.genes, array.bg$genes[, c("ProbeName", cangem.samples[i, "ProbeNames"])]))
    }
  }
  # go through all the loaded arrays and build a list of all probes found.
  all.probes <- character()
  for (i in 1:length(arrays)) {
    all.probes <- union(all.probes, names(arrays[[i]]))
  }
  # build a matrix of the measurement values
  dat <- matrix(nrow = length(all.probes), ncol = length(chips), dimnames = list(all.probes, chips))
  for (i in 1:length(arrays)) {
    dat[, chips[i]] <- arrays[[i]][all.probes]
  }
  # if the background treatment method is RMA, perform it now.
  # otherwise it has been performed already.
  if (background.treatment == "rma") {
    if (!onecolor) {
      stop("CHIPSTER-NOTE: RMA background treatment method is only available for one-color arrays.")
    }
    library(preprocessCore)
    dat <- rma.background.correct(dat, copy = TRUE)
    dat <- dat + abs(min(dat)) + 2
    rownames(dat) <- all.probes
    colnames(dat) <- chips
  }
  # in case of one-color arrays perform log2 transformation.
  # otherwise it has been perfomed as part of normalizeWithinArrays
  if (onecolor) {
    dat <- log2(dat)
  }
  # normalize between arrays.
  dat <- normalizeBetweenArrays(dat, method = inter.array.normalization)
  # average from probes to genes if needed
  if (cangem.samples[1, "ProbeNames"] != "ProbeName") {
    library(affy)
    dat <- 2^dat
    rownames(probes.to.genes) <- probes.to.genes$ProbeName
    pNList <- probes.to.genes[rownames(dat), cangem.samples[1, "ProbeNames"]]
    ngenes <- length(unique(pNList))
    pNList <- split(0:(length(pNList) - 1), pNList)
    dat <- .Call("rma_c_complete_copy", dat, pNList, ngenes, normalize = FALSE, background = FALSE, bgversion = 2, verbose = TRUE, PACKAGE = "affy")
    colnames(dat) <- chips
  }
  dat <- signif(dat, digits = 3)
  dat <- as.data.frame(dat)
  chiptype <- cangem.samples$BioconductorPackage[1]
  if (is.na(chiptype)) {
    if (cangem.samples$Type[1] == "miRNA") {
      chiptype <- "miRNA"
    } else {
      chiptype <- "cDNA"
    }
  }
} else if (cangem.samples$Format[1] == "affymetrix") {
  library(affy)
  library(gcrma)
  onecolor <- TRUE
  for (i in 1:nrow(cangem.samples)) {
    download.file(paste("http://www.cangem.org/download.php?hybridization=", cangem.samples[i, "Accession"], auth, sep = ""), cangem.samples[i, "FileName"], quiet = TRUE)
  }
  raw <- ReadAffy(filenames = cangem.samples$FileName)
  chiptype <- paste(raw@annotation, ".db", sep = "")
  calls <- as.data.frame(exprs(mas5calls(raw)))
  names(calls) <- paste("flag.", microarrays, sep = "") # or names(calls) ???
  if (affymetrix.normalization == "mas5") {
    dat <- as.data.frame(signif(log2(exprs(mas5(raw))), digits = 3))
  }
  rm(raw)
  gc()
  if (affymetrix.normalization == "rma") {
    dat <- as.data.frame(signif(exprs(justRMA(filenames = cangem.samples$FileName)), digits = 3))
  }
  if (affymetrix.normalization == "gcrma") {
    dat <- as.data.frame(signif(exprs(justGCRMA(filenames = cangem.samples$FileName, type = "fullmodel", fast = TRUE, optimize.by = "speed")), digits = 3))
  }
  names(dat) <- chips # or names(dat) ???
  dat <- data.frame(dat, calls)
  unlink(cangem.samples$FileName)
} else {
  stop("CHIPSTER-NOTE: Unsupported file format: ", cangem.samples$Format[1])
}

phenodata <- data.frame(sample = microarrays)
phenodata$original_name <- cangem.samples$FileName
phenodata$chiptype <- chiptype
phenodata$group <- ""
phenodata$description <- cangem.samples$Name
cangem.samples$ProbeNames <- NULL
cangem.samples$BioconductorPackage <- NULL
if (onecolor) {
  cangem.samples$SampleChannel <- NULL
  cangem.samples$ReferenceSample <- NULL
  cangem.samples$ReferenceSex <- NULL
  cangem.samples$ReferenceAccession <- NULL
}

for (col in colnames(cangem.samples)) {
  if (all(is.na(cangem.samples[, col]))) {
    cangem.samples[, col] <- NULL
  }
}

phenodata <- cbind(phenodata, cangem.samples)

dat2 <- data.frame(probe = rownames(dat), stringsAsFactors = FALSE)
rownames(dat2) <- dat2$probe

if (genome.build != "none") {
  # load platform
  platform <- read.table(paste("http://www.cangem.org/download.php?platform=", cangem.samples$PlatformAccession[1], "&flag=", genome.build, auth, sep = ""), sep = "\t", header = TRUE, as.is = TRUE)
  colnames(platform) <- tolower(colnames(platform))
  colnames(platform)[colnames(platform) == "chr"] <- "chromosome"
  rownames(platform) <- platform[, 1]
  platform$chromosome <- factor(platform$chromosome, levels = c(1:22, "X", "Y", "MT"), ordered = TRUE)
  dat2 <- cbind(dat2, platform[dat2$probe, c("chromosome", "start", "end")])
  dat2 <- dat2[order(dat2$chromosome, dat2$start), ]
}

if (chiptype != "cDNA" && chiptype != "miRNA") {
  # including gene names to data
  library(chiptype, character.only = T)
  lib2 <- sub(".db", "", chiptype)
  symbol <- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "SYMBOL", sep = "")))))[dat2$probe, ])
  genename <- gsub("\'", "", data.frame(unlist(as.list(get(paste(lib2, "GENENAME", sep = "")))))[dat2$probe, ])
  symbol <- gsub("#", "", symbol)
  genename <- gsub("#", "", genename)
  dat2 <- cbind(dat2, symbol, description = genename)
}

dat2 <- cbind(dat2, dat[dat2$probe, ])
if (ncol(dat) == 1) {
  colnames(dat2)[ncol(dat2)] <- chips[1]
}
dat2$probe <- NULL

writeData(dat2, "normalized.tsv")
writePhenodata(phenodata, "phenodata.tsv")

# EOF
