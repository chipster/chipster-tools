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
  #echo "user $user exists already"
  :
else
  echo "create user $user"
  kubectl exec $pod_name -- groupadd -g 1000 ubuntu
  kubectl exec $pod_name -- useradd ubuntu -u 1000 -g 1000 --create-home --shell /bin/bash
fi

if [ "$command" == "-" ]; then
  #echo "read command from stdin"
  command="$(</dev/stdin)"
fi

# create temp dir for tool installations (start with "." to omit it for example when copying artefacts in binaries.bash with wildcard "*")
kubectl exec $kubectl_exec_opts $pod_name -- su - root -c "mkdir -p $TMPDIR_PATH; chown $user:$user $TMPDIR_PATH"

echo "** run as $user in $pod_name"

init_and_command="
  set -euo pipefail
  IFS=$'\n\t'
  cd $TMPDIR_PATH
  TOOLS_PATH="$TOOLS_PATH"
  TMPDIR_PATH="$TMPDIR_PATH"
  $command"

# run strict mode to catch errors early
kubectl exec $kubectl_exec_opts $pod_name -- su - $user -c "$init_and_command"

# delete TEMP_DIR after each run, otherwise we should think about its file owners more carefully
kubectl exec $kubectl_exec_opts $pod_name -- su - root -c "rm -rf $TMPDIR_PATH"
