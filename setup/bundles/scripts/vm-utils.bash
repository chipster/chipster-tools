#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

get_name () {
  JOB_NAME=$1
  BUILD_NUMBER=$2

  adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' | tr '.' '-' | tr '[:upper:]' '[:lower:]')"
  echo "tool-install-$adjusted_job_name-$BUILD_NUMBER"
}


