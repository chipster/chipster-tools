# TOOL list_zip.py: "List contents of a zip file" (List the contents of a zip file.)
# INPUT input_file: ".zip file" TYPE GENERIC (Zip file.)
# OUTPUT output_file: "List of files in the zip package"
# RUNTIME python3

import os
import zipfile
from tool_utils import *

def main():

    document_python_version()

    input_name = read_input_definitions()['input_file']
    input_basename = remove_postfix(input_name, ".zip")
    output_name = input_basename + '_list.txt'

    with zipfile.ZipFile("input_file", "r") as zip_file:
        with open('output_file', 'w') as list_file:
            for member in zip_file.infolist():
                # member.is_dir() after python 3.6
                if member.filename[-1] == '/':
                    print('skipping directory: ' + member.filename)
                    continue

                # remove paths from dataset names, because those aren't supported in client
                list_file.write(os.path.basename(member.filename) + '\n')

    # set dataset names
    write_output_definitions({
        'output_file': output_name
    })

if __name__ == "__main__":
    main()
