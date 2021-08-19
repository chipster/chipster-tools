# TOOL test-screen-output-image.py: "Test screen output in Python and container image" (Screen output test)
# IMAGE single-shot-comp-16.04

import time

i = 20
print "Counting to %d" % (i)
for x in range(0, i ):
    print "Counting %d" % (x)
    time.sleep(1)
