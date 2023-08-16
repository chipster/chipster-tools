import sys
import os
import chipster_variables  # this file is created by PythonCompJob.java when starting the job


def document_version(application, version_string):
    # print(application, version_string)
    file = open(
        os.path.join(chipster_variables.chipster_versions_path, application + ".txt"),
        "w",
    )
    file.write(version_string)
    file.close()


def document_version_old(application, version_string):
    print("## VERSION: " + application)
    for line in version_string.splitlines():
        print("## " + line)


def document_python_version():
    python_version = (
        str(sys.version_info.major)
        + "."
        + str(sys.version_info.minor)
        + "."
        + str(sys.version_info.micro)
    )
    document_version("python", python_version)
