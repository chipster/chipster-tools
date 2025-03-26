# TOOL make_zip.py: "Make a zip package" (Create a zip package)
# INPUT input_file{...}: "Files to add" TYPE GENERIC (Files to add to the package.)
# OUTPUT package.zip: "Zip package"
# PARAMETER OPTIONAL compression: "Compression" TYPE [yes: Yes, no: No] DEFAULT no (Compressed package is smaller, but is slower to create and extract.)
# PARAMETER OPTIONAL rename: "Rename duplicate files" TYPE [yes: Yes, no: No] DEFAULT no (Rename files with unique names by appending a number when necessary.)
# RUNTIME python3

import os
import zipfile
import tool_utils

def main():
    
    if compression == "yes":
        comp_mode = zipfile.ZIP_DEFLATED
    elif compression == "no":
        comp_mode = zipfile.ZIP_STORED
    else:
        raise Exception("unknown value for compression: " + compresion)
    
    input_files = []
    client_names = {}
    
    for file in os.listdir("."):
        if file.startswith("input_file"):
            client_name = tool_utils.read_input_definitions()[file]
                
            if client_name in client_names.values():
                if rename == "no":
                    raise Exception("Duplicate file name " + client_name + ". Please rename files with unique names.")
                else:
                    print("filename " + client_name + " exists already")
                    basename, file_extension = os.path.splitext(client_name)
                    i = 2
                    while client_name in client_names.values():
                        client_name = basename + "_" + str(i) + file_extension
                        i += 1
                    print("rename to " + client_name)
            
            client_names[file] = client_name
            input_files.append(file)
            
    with zipfile.ZipFile("package.zip", 'w', comp_mode) as zipf:
        for file in input_files:        
            client_name = client_names[file]
            print("package " + file  + " as " + client_name)
            zipf.write(file, client_name)

if __name__ == "__main__":
    main()
