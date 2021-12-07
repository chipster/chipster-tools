#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/vm-utils.bash

JOB_NAME="$1"
BUILD_NUMBER="$2"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

ssh $K3S_BUILD_HOST "bash tmp/$name/clean-up.bash \"$JOB_NAME\" \"$BUILD_NUMBER\""

ssh $K3S_BUILD_HOST rm -rf tmp/$name
