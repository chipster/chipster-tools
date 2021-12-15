#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/vm-utils.bash

source="$1"
JOB_NAME="$2"
BUILD_NUMBER="$3"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

ssh $K3S_BUILD_HOST "bash tmp/$name/move-to-artefacts.bash \"$source\" \"$JOB_NAME\" \"$BUILD_NUMBER\""
