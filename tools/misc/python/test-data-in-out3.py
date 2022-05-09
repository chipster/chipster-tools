# TOOL test-data-in-out3.py: "Test data input and output in Python3" (Data input output test.) 
# INPUT input TYPE GENERIC
# OUTPUT output
# PARAMETER delay: Delay TYPE INTEGER FROM 0 TO 10000 DEFAULT 1 (Delay in seconds)
# RUNTIME python3

import shutil
import time

time.sleep(delay)
shutil.copyfile('input', 'output')

