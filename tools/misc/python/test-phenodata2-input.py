# TOOL test-phenodata-input2.py: "Test phenodata2 input in Python" () 
# INPUT META phenodata2.tsv TYPE GENERIC
# OUTPUT output.tsv

import shutil

shutil.copyfile('phenodata2.tsv', 'output.tsv')

