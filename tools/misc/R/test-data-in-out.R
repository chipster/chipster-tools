# TOOL test-data-in-out.R: "Test data input and output in R" (Data input output test.)
# INPUT input: "input" TYPE GENERIC (Just some test file)
# OUTPUT output
# PARAMETER delay: Delay TYPE INTEGER FROM 0 TO 10000 DEFAULT 1 (Delay in seconds)

system("cp input output")

print("Sleeping now!")

for (i in 1:delay) {
    print(paste("Sleeping ", i))
    Sys.sleep(1)
}

print("Waking up")
