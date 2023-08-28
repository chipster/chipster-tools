# TOOL verify-tools-bin.py: "Verify tools-bin" ()
# OUTPUT files.md5
# RUNTIME python3

import subprocess
from typing import Iterable
from datetime import datetime
import os


def get_time():
    return str(datetime.now().strftime("%H:%M:%S"))


def run_bash(cmd: str):
    run_process(["bash", "-c", cmd])


def run_process(cmd: Iterable[str]):
    process = subprocess.run(cmd)
    if process.returncode != 0:
        raise RuntimeError(
            "process failed with return code: "
            + str(process.returncode)
            + ", command: "
            + str(cmd)
        )


work_dir = os.getcwd()
os.chdir(chipster_tools_path)

file_list_path = work_dir + "/files.txt"
checksum_path = work_dir + "/files.md5"

print(get_time() + " list files")

run_bash("find . -type f | sort > " + file_list_path)

print(get_time() + " count files")

file_count = 0

with open(file_list_path) as file_list:
    for file in file_list:
        file_count += 1

print(get_time() + " file count: " + str(file_count))

print(get_time() + " calculate checksums")

checksums = 0

with open(file_list_path) as file_list:
    for file in file_list:
        # remove new line character with rstrip()
        run_bash("md5sum '" + file.rstrip() + "' >> " + checksum_path)
        checksums += 1
        if checksums % 1000 == 0:
            print(get_time() + " checksums: " + str(checksums))
