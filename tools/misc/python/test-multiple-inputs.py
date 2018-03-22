# TOOL test-multiple-inputs.py: "Test multiple inputs in Python" (Multiple input output test.) 
# INPUT input{...} TYPE GENERIC
# INPUT other TYPE GENERIC
# OUTPUT output
# OUTPUT OPTIONAL missing_output.txt

import shutil

shutil.copyfile('other', 'output')

