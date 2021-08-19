# TOOL test-data-in-out-image.py: "Test data input and output in Python and container image" (Data input output test.) 
# INPUT input TYPE GENERIC
# OUTPUT output
# OUTPUT OPTIONAL missing_output.txt
# IMAGE single-shot-comp-16.04

import shutil

shutil.copyfile('input', 'output')

