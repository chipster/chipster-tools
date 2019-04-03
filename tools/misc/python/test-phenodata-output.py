# TOOL test-phenodata-output.py: "Test phenodata output in Python" ()
# INPUT input{...}.tsv TYPE GENERIC
# OUTPUT output.tsv
# OUTPUT META phenodata.tsv

with open('output.tsv', 'w') as f:
	f.write('identifier	chip.sample1\n')
	f.write('test	output\n')

with open('phenodata.tsv', 'w') as f:
	f.write('sample	chiptype	group\n')
	f.write('sample001.tsv	chiptype	1\n')
