# TOOL test-phenodata-output.py: "Test phenodata output in Python" () 
# OUTPUT phenodata.tsv
# OUTPUT output.tsv

with open('output.tsv', 'w') as f:
	f.write('test	output\n')
	f.write('test	output\n')

with open('phenodata.tsv', 'w') as f:
	f.write('sample	original_name	chiptype	group\n')
	f.write('microarray001.cel	cancerGSM11814.cel	hgu133ahsentrezg.db	1\n')
	f.write('microarray002.cel	cancerGSM11830.cel	hgu133ahsentrezg.db	2\n')
