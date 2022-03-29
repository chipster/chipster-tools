# TOOL test-data-in-out-without-tools-bin.py: "Test data input and output without tools-bin" (Data input output test.) 
# INPUT input TYPE GENERIC
# OUTPUT output
# PARAMETER delay: Delay TYPE INTEGER FROM 0 TO 10000 DEFAULT 1 (Delay in seconds)
# TOOLS_BIN ""

import shutil
import time

time.sleep(delay)
shutil.copyfile('input', 'output')

