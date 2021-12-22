#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# export these to environment variables so that these are available in ../scripts/run-in-pod.bash
export TOOLS_PATH="/opt/chipster/tools"
export TMPDIR_PATH="/opt/chipster/tools/.tmp"
