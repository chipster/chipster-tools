# TOOL test-phenodata-output.py: "Test phenodata output in Python" ()
# INPUT input{...}.tsv TYPE GENERIC
# OUTPUT phenodata.tsv
# OUTPUT output.tsv

with open('output.tsv', 'w') as f:
	f.write('identifier	chip.sample1\n')
	f.write('test	output\n')

with open('phenodata.tsv', 'w') as f:
	f.write('dataset	column	sample	chiptype	experiment	group	library_size\n')
	f.write('ngs-data-table.tsv	chip.sample000.tsv	sample000.tsv	not applicable	rna_seq		\n')
	f.write('ngs-data-table.tsv	chip.sample001.tsv	sample001.tsv	not applicable	rna_seq		\n')
	