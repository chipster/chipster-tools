# TOOL gunzip.py: "Extract .gz file" (Extract a gzip file, which usually has a file extension .gz)
# INPUT input_file: "Gzip file" TYPE GENERIC (Gzip compressed file)
# OUTPUT output_file: "Extracted file"
# RUNTIME python3

import gzip
import shutil
import tool_utils


def main():
    infile = gzip.open("input_file", "rb")

    with open("output_file", "wb") as outfile:
        # copy in chunks
        shutil.copyfileobj(infile, outfile)
    infile.close()

    # set dataset name
    input_name = tool_utils.read_input_definitions()["input_file"]
    output_names = {"output_file": tool_utils.remove_postfix(input_name, ".gz")}

    tool_utils.write_output_definitions(output_names)


if __name__ == "__main__":
    main()
