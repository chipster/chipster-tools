# TOOL example.R: "Provide template for tools" (Provide template for tools)
# PARAMETER OPTIONAL int: Integer TYPE INTEGER FROM -100 TO 100 DEFAULT 0
# PARAMETER OPTIONAL dec: Decimal TYPE DECIMAL FROM -0.1 TO 1.0 DEFAULT 0
# PARAMETER OPTIONAL string: String TYPE STRING
# PARAMETER OPTIONAL unchecked: "Unchecked string" TYPE UNCHECKED_STRING
# PARAMETER OPTIONAL enum: Enum TYPE [first: First, second: Second, third: Third] DEFAULT first
# PARAMETER OPTIONAL column: Column TYPE COLUMN_SEL
# PARAMETER OPTIONAL metacolumn: "Phenodata column" TYPE METACOLUMN_SEL
# PARAMETER OPTIONAL file: File TYPE INPUT_SEL

# VARIABLES, examples values
chipster.tools.path # "/opt/chipster-web-server/tools"
chipster.common.path # "../toolbox/tools/common/R"
chipster.common.lib.path # "../toolbox/tools/common/R/lib"
chipster.module.path # "../toolbox/tools/misc"
chipster.module.lib.path # "../toolbox/tools/misc/lib"
chipster.java.libs.path # "/opt/chipster-web-server/lib"
chipster.threads.max # "2"
chipster.memory.max # "8192"




# include utils
source(file.path(chipster.common.lib.path, "tool-utils.R"))


# source(file.path(chipster.common.lib.path, "version-utils.R"))
# documentVersion("R", R.version.string)





# inputs
# - multi
# - all types
# - same display name?
# - empty description?

# how to call a function from own lib file


# chipster-note

# output rename

# tool versions

# unzip
