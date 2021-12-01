#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source ~/jenkins-env.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 2 ]; then
  echo "Clean-up tool installation pod"
  echo "Usage: $0 JOB_NAME BUILD_NUNMBER"
  exit 0
fi

JOB_NAME="$1"
BUILD_NUMBER="$2"

adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' | tr '.' '-')"
name="tool-install-$adjusted_job_name-$BUILD_NUMBER"

kubectl delete deployment $name

if grep -qs "/mnt/data/$name/tools-bin " /proc/mouns; then
  echo "umount tools-bin"
  sudo umount /mnt/data/$name/tools-bin
else
  echo "/mnt/data/$name/tools-bin is not a mount"
fi

echo "delete /mnt/data/$name"
sudo rm -rf /mnt/data/$name
