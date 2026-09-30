# TOOL test-r-in-container-without-tools-bin.R: "Test R in container image without tools-bin" (Data input output test.)
# INPUT input: "input" TYPE GENERIC (Just some test file)
# OUTPUT output
# RUNTIME R-4.4.3
# TOOLS_BIN ""

system("cp input output")
print("hello! hei!")
