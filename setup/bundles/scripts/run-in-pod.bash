#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/vm-utils.bash

JOB_NAME="$1"
BUILD_NUMBER="$2"
user="$3"
command="$4"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

cat /dev/stdin | ssh $K3S_BUILD_HOST "bash tmp/$name/run-in-pod.bash \"$JOB_NAME\" \"$BUILD_NUMBER\" \"$user\" \"$command\""
