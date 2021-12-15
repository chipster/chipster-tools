#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/vm-utils.bash

JOB_NAME="$1"
BUILD_NUMBER="$2"
image="$3"
BUNDLE_COLLECTION_VERSION="$4"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

ssh $K3S_BUILD_HOST mkdir -p tmp/$name

rsync -r $BUNDLE_SCRIPTS_DIR/k3s-build-host/ $K3S_BUILD_HOST:tmp/$name

ssh $K3S_BUILD_HOST "bash tmp/$name/start-pod.bash \"$JOB_NAME\" \"$BUILD_NUMBER\" \"$image\" \"$BUNDLE_COLLECTION_VERSION\""
