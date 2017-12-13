# TOOL test-data-in-out.py: "Test data input and output in Python" (Data input output test.) 
# INPUT input TYPE GENERIC
# OUTPUT output
# OUTPUT OPTIONAL missing_output.txt

import shutil

shutil.copyfile('input', 'output')

