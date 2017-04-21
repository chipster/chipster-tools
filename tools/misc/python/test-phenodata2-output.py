# TOOL test-phenodata2-output.py: "Test phenodata2 output in Python" () 
# OUTPUT phenodata2.tsv
# OUTPUT output.tsv

with open('output.tsv', 'w') as f:
	f.write('test	output\n')
	f.write('test	output\n')

with open('phenodata2.tsv', 'w') as f:
	f.write('dataset	column	sample	original_name	chiptype	group\n')
	f.write('output.tsv	chip.microarray001.cel	microarray001.cel	cancerGSM11814.cel	hgu133ahsentrezg.db	1\n')
	f.write('output.tsv	chip.microarray002.cel	microarray002.cel	cancerGSM11830.cel	hgu133ahsentrezg.db	2\n')
