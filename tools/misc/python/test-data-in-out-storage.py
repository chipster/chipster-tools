# TOOL test-data-in-out-storage.py: "Test storage volume for working directory" (Data input output test.) 
# INPUT input TYPE GENERIC
# OUTPUT output
# PARAMETER delay: Delay TYPE INTEGER FROM 0 TO 10000 DEFAULT 1 (Delay in seconds)
# STORAGE 10

import shutil
import time

time.sleep(delay)
shutil.copyfile('input', 'output')

