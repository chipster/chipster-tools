#!/bin/bash

set -e

source ~/jenkins-env.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 4 ]; then
  echo "Run command in the tool installation pod"
  echo "Usage: $0 JOB_NAME BUILD_NUNMBER USER COMMAND"
  exit 0
fi

JOB_NAME="$1"
BUILD_NUMBER="$2"
user="$3"
command="$4"
kubectl_exec_opts="-it"

adjusted_job_name="$(echo "$JOB_NAME" | tr '_' '-' )"
name="tool-install-$adjusted_job_name-$BUILD_NUMBER"

pod_name=$(kubectl get pod | grep $name | grep Running | cut -d " " -f 1)

if kubectl exec $pod_name -- id -u $user >> /dev/null 2>&1; then
  echo "user $user exists already"
else
  echo "create user $user"
  kubectl exec $pod_name -- groupadd -g 1000 ubuntu
  kubectl exec $pod_name -- useradd ubuntu -u 1000 -g 1000 --create-home --shell /bin/bash
fi

if [ "$command" == "-" ]; then
  echo "read command from stdin"
  command=$(</dev/stdin)
  kubectl_exec_opts=""
fi

echo "** running in $image container as $user"
kubectl exec $kubectl_exec_opts $pod_name -- su - -c "set -e; $command"
