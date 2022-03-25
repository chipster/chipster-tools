#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source ~/jenkins-env.bash

get_name () {
  JOB_NAME=$1
  BUILD_NUMBER=$2

  adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' | tr '.' '-' | tr '[:upper:]' '[:lower:]')"
  echo "tool-install-$adjusted_job_name-$BUILD_NUMBER"
}

get_time () {
  date +%s
}

show_elapsed () {
  start_time="$1"
  elapsed=$(( $(get_time) - start_time ))
  eval "echo Elapsed time: $(date -ud "@$elapsed" +'$((%s/3600/24)) days %H hr %M min %S sec')"
}

