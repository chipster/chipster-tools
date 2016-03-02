# TOOL test-phenodata-input.py: "Test phenodata input in Python" () 
# INPUT META phenodata.tsv TYPE GENERIC
# OUTPUT output.tsv

import shutil

shutil.copyfile('phenodata.tsv', 'output.tsv')

