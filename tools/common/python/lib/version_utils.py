import sys

def document_version(application, version_string):
    print('## VERSION: ' + application)
    for line in version_string.splitlines():
        print('## ' + line)

def document_python_version():
    python_version = str(sys.version_info.major) + '.' + str(sys.version_info.minor) + '.' + str(sys.version_info.micro)
    document_version("python", python_version)

