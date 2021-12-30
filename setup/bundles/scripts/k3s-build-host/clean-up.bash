#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source ~/jenkins-env.bash
source $(dirname "$0")/bundle-utils.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 2 ]; then
  echo "Clean-up tool installation pod"
  echo "Usage: $0 JOB_NAME BUILD_NUNMBER"
  exit 0
fi

JOB_NAME="$1"
BUILD_NUMBER="$2"

name=$(get_name $JOB_NAME $BUILD_NUMBER)

kubectl delete deployment $name

if grep -qs "/mnt/data/$name/tools-bin " /proc/mounts; then
  echo "umount tools-bin"
  sudo umount /mnt/data/$name/tools-bin
else
  #echo "/mnt/data/$name/tools-bin is not a mount"
  :
fi

#echo "delete /mnt/data/$name"
sudo rm -rf /mnt/data/$name
