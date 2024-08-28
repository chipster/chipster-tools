# TOOL annotate-ensembl-ids.py: "Annotate Ensembl identifiers" (Annotates Ensembl IDs with gene symbols and descriptions, creates a new table containing these and the values in the original input file. The Ensembl IDs can be either the rownames or in the first column of the input table. They can also be in the middle if the name of the column is ensembl_id.)
# INPUT genelist.tsv: genelist.tsv TYPE GENERIC
# OUTPUT annotated.tsv: annotated.tsv
# RUNTIME python3
# TOOLS_BIN ""

import csv
import requests
import json


# https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    # for i in range(0, len(lst), n):
    #     yield lst[i:i + n]
    return [lst[i : i + n] for i in range(0, len(lst), n)]


def annotate(id_list):
    req_body = {"ids": id_list}

    request_json = json.dumps(req_body)

    server = "http://rest.ensembl.org"
    ext = "/lookup/id"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    url = server + ext

    response = requests.post(url, headers=headers, data=request_json)

    if not response.ok:
        if response.status_code == 429:
            # rate limiting, print some info from headers
            print("Retry-After: " + str(response.headers.get("Retry-After")))
            print(
                "X-RateLimit-Limit: " + str(response.headers.get("X-RateLimit-Limit"))
            )
            print(
                "X-RateLimit-Reset: " + str(response.headers.get("X-RateLimit-Reset"))
            )
            print(
                "X-RateLimit-Period: " + str(response.headers.get("X-RateLimit-Period"))
            )
            print(
                "X-RateLimit-Remaining: "
                + str(response.headers.get("X-RateLimit-Remaining"))
            )

        response.raise_for_status()

    return response.json()


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

# use column "ensembl_id" if it exists, otherwise use the first column
id_col_index = 0
id_col_title = None
other_col_titles = None

try:
    column_to_find = "ensembl_id"

    id_col_index = column_titles.index(column_to_find)
    id_col_title = column_titles[id_col_index]
    other_col_titles = column_titles.copy()
    other_col_titles.pop(id_col_index)

    if no_first_column_title:
        # the first column title is missing i.e. R rowname
        id_col_index = id_col_index + 1
        print(
            "the first column does not have title, new column index: "
            + str(id_col_index)
        )

except ValueError:
    print("column '" + column_to_find + "' not found, using first")

    if no_first_column_title:
        other_col_titles = column_titles
    else:
        id_col_title = column_titles[0]
        other_col_titles = column_titles.copy()
        other_col_titles.pop(0)

# there is a limit how many IDs can be queried in one request
max_ids_per_query = 1000

# limit how many request are allowed in one job to avoid Rest API rate limits
# max_chunks = 100 # increased to 300 to allow customer to run with bigger inputs
max_chunks = 300

row_chunks = chunks(rows, max_ids_per_query)

print("going to make " + str(len(row_chunks)) + " request(s)")

if len(row_chunks) > max_chunks:
    raise RuntimeError(
        "CHIPSTER-NOTE: Max "
        + str(max_chunks * max_ids_per_query)
        + " genes are allowed in one job. The input file has "
        + str(len(rows))
        + " genes."
    )

print("open the results file")

with open("annotated.tsv", "w") as annotated:
    print("write column titles")

    if no_first_column_title:
        # the first column title is missing i.e. R rowname
        annotated.write("symbol\tdescription\t" + "\t".join(column_titles) + "\n")
    else:
        print("keep the column title of first column: " + str(id_col_title))
        # keep the colum title of first column
        annotated.write(
            id_col_title
            + "\tsymbol\tdescription\t"
            + "\t".join(other_col_titles)
            + "\n"
        )

    for row_chunk in row_chunks:
        print("annotate " + str(len(row_chunk)) + " gene(s)")

        not_found_ids = []

        # pick the id column
        id_list = [row[id_col_index] for row in row_chunk]

        response = annotate(id_list)

        for row in row_chunk:
            ensemblid = row[id_col_index]
            other_columns = row.copy()
            other_columns.pop(id_col_index)

            id_response = response.get(ensemblid)

            if not id_response:
                not_found_ids.append(ensemblid)
                # deprecated ID
                print("ensemblid not found: " + ensemblid)
                # write the existing row anyway
                id_response = {"symbol": "", "description": ""}

            # empty string when there is a response, but it doesn't have display_name or description
            # (e.g. novel_transcripts)
            symbol = id_response.get("display_name", "")
            description = id_response.get("description", "")

            # print(ensemblid, symbol, description)
            annotated.write(
                ensemblid
                + "\t"
                + symbol
                + "\t"
                + description
                + "\t"
                + "\t".join(other_columns)
                + "\n"
            )

        # give up if most of the IDs in this chunk were not found
        if len(not_found_ids) > 0.9 * len(row_chunk):
            raise RuntimeError(
                "CHIPSTER-NOTE: "
                + str(len(not_found_ids))
                + "/"
                + str(len(row_chunk))
                + " IDs not found. You can only annotate Ensembl IDs with this tool. The IDs need to be in the first column of the table."
            )
