# TOOL extract_zip.py: "Extract .zip file" (Extract a zip file (usually has a file extension .zip.\)
# All files are extracted by default, but you can choose to extract only one or some of the files by giving the names of the files to be extracted as a parameter. Alternatively you can provide a list of filenames as a text file, one name per line. 
# To see the contents of a .zip file, use tool \"Misc / Utilities / List contents of a zip file\".)
# INPUT input_file: ".zip file" TYPE GENERIC (Zip file.)
# INPUT OPTIONAL file_list: "List of files to extract" TYPE GENERIC (List of files to extract. One filename per line.)
# OUTPUT output_file{...}: "Extracted file(s)"
# PARAMETER OPTIONAL names: "Extract by filename" TYPE STRING (Filenames of the files to extract. If more than one, separate names with a comma (e.g. abc123_1.fq,abc123_2.fq\). Alternatively you can provide a list of filenames as a text file, one name per line.)
# PARAMETER OPTIONAL extensions: "Extract by extension" TYPE STRING (Filename extension of the files to extract (e.g .html\). If more than one, separate with a comma (e.g .txt,.log\).)
# PARAMETER OPTIONAL full_paths: "Keep directories" TYPE [yes: Yes, no: No] DEFAULT no (Use the whole file path for the filename.)
# RUNTIME python3

import os
import zipfile
import shutil
import subprocess
from tool_utils import *

def main():

    document_python_version()
    # discard stderr. zipinfo in Ubuntu 16.04 expects to get the zip file 
    # always and prints an error when we don't give one
    document_version("zipinfo", subprocess.check_output('zipinfo 2>/dev/null | head -n 1', shell=True).decode('utf-8'))
    document_version("grep", subprocess.check_output('grep --version | head -n 1', shell=True).decode('utf-8'))

    output_names = {}
    include_list = []
    include_list_raw = []
    # by default, extract all
    someonly = False

    # Filename list
    if names:
        include_list_raw = names.split(',')
        someonly = True
        
    # Extension list
    if extensions:
        someonly = True
        os.system("zipinfo -1 input_file > toc")
        # Remove any * characters as they are treated literally and not as a wildcard.
        clean_ext = extensions.translate(str.maketrans('','', '*'))
        extension_list = clean_ext.split(',')
        for ext in extension_list:
            os.system("grep -i " + ext +"$ toc >> selected_list")
        with open('selected_list','r') as f:
            include_list_raw = [l.strip() for l in f]  

    # File of filenames
    if os.path.isfile('file_list'):
        someonly = True
        with open('file_list','r') as f:
            include_list_raw = [l.strip() for l in f]
            
   
    for name in include_list_raw:
        include_list.append(os.path.basename(name.strip()))

    # it's nice to do this with a python library, but it is also   
    # faster than the unzip command, at least in macOS 
    # (decompression of 4 GB session, 22 vs. 32 seconds)

    with zipfile.ZipFile("input_file", "r") as f:        
        for member in f.infolist():
            # member.is_dir() after python 3.6
            if member.filename[-1] == '/':
                print('skipping directory: ' + member.filename)
                continue

            if someonly:
                if not os.path.basename(member.filename) in include_list:
                    print('skipping, not in list: ' + member.filename)
                    continue

            # fixed names for output files
            output_file = 'output_file' + str(len(output_names) + 1)

            # dataset name for the client
            dataset_name = member.filename
            if full_paths == 'no':
                # remove paths from dataset names, because those aren't supported in client
                dataset_name = os.path.basename(dataset_name)
            output_names[output_file] = dataset_name

            # extract without paths, because those could point to parent directories
            with open(output_file, 'wb') as outfile, f.open(member) as infile:
                shutil.copyfileobj(infile, outfile)

    # set dataset names
    write_output_definitions(output_names)

if __name__ == "__main__":
    main()
