#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Run a script in tool installation pod
# 
# This is the stable interface that all build scripts should call. We should be able to change where and how the 
# containers are run without changing all build scripts. 

source $(dirname "$0")/vm-utils.bash

JOB_NAME="$1"
BUILD_NUMBER="$2"
user="$3"
command="$4"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

cat /dev/stdin | ssh $K3S_BUILD_HOST "\
    TOOLS_PATH="$TOOLS_PATH" \
    TMPDIR_PATH="$TMPDIR_PATH" \
    bash tmp/$name/run-in-pod.bash \"$JOB_NAME\" \"$BUILD_NUMBER\" \"$user\" \"$command\""
