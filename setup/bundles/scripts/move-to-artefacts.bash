#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Move folder from the tool installation pod to artefacts
# 
# This is the stable interface that all build scripts should call. 
#
# See start-pod.bash for longer explanation.

source $(dirname "$0")/vm-utils.bash

source="$1"
JOB_NAME="$2"
BUILD_NUMBER="$3"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

ssh $K3S_BUILD_HOST "\
    TOOLS_PATH="$TOOLS_PATH" \
    TMPDIR_PATH="$TMPDIR_PATH" \
    bash tmp/$name/move-to-artefacts.bash \"$source\" \"$JOB_NAME\" \"$BUILD_NUMBER\""
