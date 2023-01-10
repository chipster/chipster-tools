# TOOL remove-duplicate-genes.py: "Remove duplicate genes" (This tool takes a list of identifiers and removes duplicates. The output file is compatible with downstream tools used for example for pathway analysis.)
# INPUT genelist.tsv: "Table with identifiers" TYPE GENERIC 
# OUTPUT unique-genes.tsv: "Table listing the unique genes." 
# PARAMETER id_col_title: "ID column" TYPE COLUMN_SEL DEFAULT ensembl_id (If multiple rows have the same value in this column, only one of them is kept.)
# RUNTIME python3
# TOOLS_BIN ""

import csv
import requests
import json

print("parse input file")

column_titles = None
rows = []

with open("genelist.tsv") as fd:
	rd = csv.reader(fd, delimiter="\t", quotechar='"')
	# skip column title row
	column_titles = next(rd)

	for row in rd:
		rows.append(row)

print("found " + str(len(rows)) + " rows")

no_first_column_title = len(column_titles) + 1 == len(rows[0])

def get_column_index(column_to_find):

	print("column titles: " + str(column_titles))

	if column_to_find == "":
		print("column is empty string, using rownames")
		return 0

	if column_to_find == " ":
		print("column is a space character, using rownames")
		return 0

	print("trying to find column " + column_to_find)
	col_index = column_titles.index(column_to_find)
	print("column index: " + str(col_index))

	if no_first_column_title:
		# the first column title is missing i.e. R rowname
		col_index = col_index + 1
	
	return col_index

id_col_index = 0
loc_col_index = None

try:
	id_col_index = get_column_index(id_col_title)

except ValueError:
	# this shouldn't be possible, because the column comes from the COLUMN_SEL parameter
	raise RuntimeError("CHIPSTER-NOTE: column " + id_col_title + " not found")
	
try:
	loc_col_index = get_column_index("location")
	
except ValueError:
	print("column 'location' not found, ignore")

	
print("sort")

# priorities of the location values
loc_prio = [
	"upstream",
	"overlapStart",
	"inside",
	"overlapEnd",
	"downstream",
]

def compareLocation(row):
	location = row[loc_col_index]
	try:
		return loc_prio.index(location)
	except ValueError:
		# lower priority than all known values
		return len(loc_prio)

def compareId(row):
	return row[id_col_index]

if loc_col_index is not None:
	print("sort by location")
	rows.sort(key=compareLocation)

print("sort by ID")
# preserves the order by location among same IDs
rows.sort(key=compareId)

print("remove duplicates")
unique_ids = set()
unique_rows = []

for row in rows:
	id = row[id_col_index]

	if id not in unique_ids:
		unique_ids.add(id)
		unique_rows.append(row)

print("open the results file")

with open("unique-genes.tsv", "w") as out_file:

	print("write column titles")

	# write original column titles
	out_file.write("\t".join(column_titles) + "\n")

	for row in unique_rows:
		out_file.write("\t".join(row) + "\n")
		