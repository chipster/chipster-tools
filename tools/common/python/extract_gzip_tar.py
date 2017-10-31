# TOOL extract_gzip_tar.py: "Extract .tar or .tar.gz file" (Extract a tar file. The file can be gzipped (usually has a file extension .tar.gz.\)
# All files are extracted by default, but you can choose to extract only one or some of the files by giving the names of the files to be extracted as a parameter. Alternatively you can provide a list of filenames as a text file, one name per line. 
# To see the contents of a .tar file, use tool \"Utilities: List contents of a tar file\".)
# INPUT input_file: ".tar.gz file" TYPE GENERIC (Tar file. Can be gzip compressed.)
# INPUT OPTIONAL file_list: "List of files to extract" TYPE GENERIC (List of files to extract. One Filename per line.)
# OUTPUT output_file{...}: "Extracted file(s)"
# PARAMETER OPTIONAL names: "Extract by filename" TYPE STRING (Filenames of the files to extract. If more than one, separate names with a comma (e.g. abc123_1.fq,abc123_2.fq\). Alternatively you can provide a list of filenames as a text file, one name per line.)
# PARAMETER OPTIONAL extensions: "Extract by extension" TYPE STRING (Filename extension of the files to extract (e.g .html\). If more than one, separate with a comma (e.g .txt,.log\).)

import os
import sys
import tarfile
import posixpath
from tool_utils import *

def main():  

    input_names = read_input_definitions()
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
        os.system("tar tf input_file > toc")
        extension_list = extensions.split(',')
        for ext in extension_list:
            os.system("grep " + ext +"$ toc >> selected_list")
        with open('selected_list','r') as f:
            include_list_raw = [l.strip() for l in f]  

    # File of filenames
    if os.path.isfile('file_list'):
        someonly = True
        with open('file_list','r') as f:
            include_list_raw = [l.strip() for l in f]
            
   
    for name in include_list_raw:
        include_list.append(os.path.basename(name.strip()))
        
    # benchmark this in comparison to gunzip and pigz
    infile = tarfile.open('input_file', 'r')

    for member in infile.getmembers():
        if not member.isfile():
            print 'skipping, not a file: ' + member.name
            continue
        if someonly:
            if not os.path.basename(member.name) in include_list:
                print 'skipping, not in list: ' + member.name
                continue
        # fixed names for output files
        output_file = 'output_file' + str(len(output_names) + 1)
        # remove paths from dataset names, because those aren't supported in client
        output_names[output_file] = posixpath.basename(member.name)
        # extract without paths, because those could point to parent directories
        member.name = output_file
        infile.extract(member, '.')

    infile.close()

    # set dataset names
    write_output_definitions(output_names)

if __name__ == "__main__":
    main()
