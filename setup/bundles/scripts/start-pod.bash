#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Start a tool installation pod
# 
# This is the stable interface that all build scripts should call. 
# 
# We should be able to change where and how the containers are run without changing all build scripts. 
# At the moment we have a separate host for running containers. Let's call it K3s host. However, it's easier to 
# manage all scripts in one place on this Jenkins host. Thus we first copy the scripts to the K3s host and then
# run another script there which does the actual work.

source $(dirname "$0")/vm-utils.bash

JOB_NAME="$1"
BUILD_NUMBER="$2"
image="$3"
BUNDLE_COLLECTION_VERSION="$4"

name="$(get_name $JOB_NAME $BUILD_NUMBER)"

ssh $K3S_BUILD_HOST mkdir -p tmp/$name

rsync -r $BUNDLE_SCRIPTS_DIR/k3s-build-host/ $K3S_BUILD_HOST:tmp/$name

ssh $K3S_BUILD_HOST "bash tmp/$name/start-pod.bash \"$JOB_NAME\" \"$BUILD_NUMBER\" \"$image\" \"$BUNDLE_COLLECTION_VERSION\""
