# Prints out version information.
#
printVersion <- function(application, version.string) {
  cat("##", "VERSION:", application, "\n")
  cat("##", version.string, "\n")
}

documentVersion <- function(application, version.string) {
  write(version.string, file = file.path(chipster.versions.path, paste(application, ".txt", sep = "")))
}

documentRVersion <- function() {
  documentVersion("R", R.version.string)
}

documentRSessionInfo <- function() {
  documentVersion("R", R.version.string)
  documentVersion("R sessionInfo", capture.output(sessionInfo()))
}


# This is called from java code RCompJob.java
#
documentR <- function() {
  documentRVersion()
  documentRSessionInfo()
}
