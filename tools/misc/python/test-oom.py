# TOOL test-oom.py: "Test out of memory handling" ()
# INPUT input TYPE GENERIC
# OUTPUT output
# PARAMETER size: "Size in GiB" TYPE INTEGER FROM 0 TO 100 DEFAULT 1 ()
# RUNTIME python3

import shutil
import array
import time

arrays = []

for i in range(size):
    print("reserve memory " + str(i))
    for j in range(1024):
        arrays.append(array.array('b', [0])*1024*1024)

print("sleep to wait for memory monitoring")
time.sleep(10)

print("chipster_threads_max: " + chipster_threads_max)
print("chipster_memory_max: " + chipster_memory_max)
shutil.copyfile("input", "output")
