# TOOL test-screen-output.py: "Test screen output in Python" (Screen output test) 

import time

i = 20
print "Counting to %d" % (i)
for x in range(0, i ):
    print "Counting %d" % (x)
    time.sleep(1)
