#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source ~/jenkins-env.bash
source $(dirname "$0")/bundle-utils.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 4 ]; then
  echo "Run command in the tool installation pod"
  echo "Usage: $0 JOB_NAME BUILD_NUNMBER USER COMMAND"
  exit 0
fi

JOB_NAME="$1"
BUILD_NUMBER="$2"
user="$3"
command="$4"
# for interactive use
#kubectl_exec_opts="-it"
kubectl_exec_opts=""

name=$(get_name $JOB_NAME $BUILD_NUMBER)

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
  command="$(</dev/stdin)"
fi

# create temp dir for tool installations
TEMP_DIR="/opt/chipster/tools/tmp"
kubectl exec $kubectl_exec_opts $pod_name -- su - root -c "echo create $TEMP_DIR; mkdir -p $TEMP_DIR; chown $user:$user $TEMP_DIR"

echo "** run as $user"

# run strict mode to catch errors early
kubectl exec $kubectl_exec_opts $pod_name -- su - $user -c "
  set -euo pipefail
  IFS=$'\n\t'
  cd $TEMP_DIR
  $command"

kubectl exec $kubectl_exec_opts $pod_name -- su - root -c "echo delete $TEMP_DIR; rm -rf $TEMP_DIR"
