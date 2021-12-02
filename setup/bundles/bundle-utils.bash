#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source ~/jenkins-env.bash

get_name () {
  JOB_NAME=$1
  BUILD_NUMBER=$2

  adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' | tr '.' '-' )"
  echo "tool-install-$adjusted_job_name-$BUILD_NUMBER"
}


