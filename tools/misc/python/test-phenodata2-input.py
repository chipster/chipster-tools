# TOOL test-phenodata2-input.py: "Test phenodata2 input in Python" ()
# INPUT input.tsv TYPE GENERIC
# INPUT META phenodata2.tsv TYPE GENERIC
# OUTPUT output.tsv

import shutil

shutil.copyfile('phenodata2.tsv', 'output.tsv')

