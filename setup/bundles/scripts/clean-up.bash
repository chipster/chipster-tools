#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Clean up after the tool installation
# 
# This is the stable interface that all build scripts should call. 
#
# See start-pod.bash for longer explanation.

source $(dirname "$0")/vm-utils.bash

JOB_NAME="$1"
BUILD_NUMBER="$2"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

ssh $K3S_BUILD_HOST "bash tmp/$name/clean-up.bash \"$JOB_NAME\" \"$BUILD_NUMBER\""

ssh $K3S_BUILD_HOST rm -rf tmp/$name
