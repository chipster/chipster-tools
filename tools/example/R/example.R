# TOOL example.R: "R template for tools" (R template for tools)
# PARAMETER OPTIONAL int: Integer TYPE INTEGER FROM -100 TO 100 DEFAULT 0
# PARAMETER OPTIONAL dec: Decimal TYPE DECIMAL FROM -0.1 TO 1.0 DEFAULT 0
# PARAMETER OPTIONAL string: String TYPE STRING
# PARAMETER OPTIONAL unchecked: "Unchecked string" TYPE UNCHECKED_STRING
# PARAMETER OPTIONAL enum: Enum TYPE [first: First, second: Second, third: Third] DEFAULT first
# PARAMETER OPTIONAL column: Column TYPE COLUMN_SEL
# PARAMETER OPTIONAL metacolumn: "Phenodata column" TYPE METACOLUMN_SEL
# PARAMETER OPTIONAL file: File TYPE INPUT_SEL

# Chipster variables
# ------------------------------------------------------------------------------
# Variables created by the R job. Use these to for example constract paths to
# applications and other R scripts.
chipster.tools.path # "/opt/chipster-web-server/tools"
chipster.common.path # "../toolbox/tools/common/R"
chipster.common.lib.path # "../toolbox/tools/common/R/lib"
chipster.module.path # "../toolbox/tools/misc"
chipster.module.lib.path # "../toolbox/tools/misc/lib"
chipster.java.libs.path # "/opt/chipster-web-server/lib"
chipster.threads.max # "2"
chipster.memory.max # "8192"


# Use a library script file
# ------------------------------------------------------------------------------
# Library R files should go to tools/common/R/lib.
# Load the library with:
source(file.path(chipster.common.lib.path, "tool-utils.R"))

# Use a function from the library script (no need for prefix or anything):
read_input_definitions()


# Document application version
# ------------------------------------------------------------------------------
# Store the version of the application used by a script. Application version
# information is stored in with the job in the session db and it's visible for
# the users in the job details.
#
# Version functions are implemented in the version-utils.R, which is always
# sourced at the start of a R job, so no need to explicitly source
# version-utils.R.
#
# First argument is the name of the application and the second is the string
# containing version information. The string can have multiple lines.
documentVersion("bowtie", version.string)


# Notify user with Chipster note
# ------------------------------------------------------------------------------
# Use Chipster note to notify the user in a situation where running a script
# can't continue, but the user can do something to fix the problem.
#

if (param.count < 2) {
    stop(paste("CHIPSTER-NOTE:", "Count needs to greater than 2"))
}




# inputs
# - multi
# - all types
# - same display name?
# - empty description?


# output rename


# unzip
