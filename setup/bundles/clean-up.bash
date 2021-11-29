#!/bin/bash

set -e

source ~/jenkins-env.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 2 ]; then
  echo "Clean-up tool installation pod"
  echo "Usage: $0 JOB_NAME BUILD_NUNMBER"
  exit 0
fi

JOB_NAME="$1"
BUILD_NUMBER="$2"

adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' )"
name="tool-install-$adjusted_job_name-$BUILD_NUMBER"

kubectl delete deployment $name

echo "umount tools-bin"
sudo umount /mnt/data/$name/tools-bin

echo "delete /mnt/data/$name"
sudo rm -rf /mnt/data/$name
