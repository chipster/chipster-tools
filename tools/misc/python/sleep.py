# TOOL sleep.py: "Sleep in Python" ()
# PARAMETER delay: Delay TYPE INTEGER FROM 0 TO 1000 DEFAULT 10 (Delay in seconds)

import time
# cannot use tool_utils in python2

print("Start : %s" % time.ctime())
time.sleep(delay)
print("End : %s" % time.ctime())
