#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source ~/jenkins-env.bash
source $(dirname "$0")/bundle-utils.bash

if [ "$1" == "-h" ] || [ "$1" == "--help" ] || [ $# -ne 3 ]; then
  echo "Move directory to artefacts"
  echo "Usage: $0 SOURCE JOB_NAME BUILD_NUNMBER"
  exit 0
fi

source="$1"
JOB_NAME="$2"
BUILD_NUMBER="$3"

# different path on nfs clients and server
nfs_client_bundle="/mnt/artefacts/$JOB_NAME/$BUILD_NUMBER"
nfs_server_bundle="/export/artefacts/$JOB_NAME/$BUILD_NUMBER"

# copying 10k files to nfs takes about 6 minutes
#time mv /opt/chipster/tools/umi-tools $bundle

echo "package $source"
start_time=$(get_time)
    
# packaging and extracting takes less than 10 seconds
TOOLS_PATH="$TOOLS_PATH" \
TMPDIR_PATH="$TMPDIR_PATH" \
bash $(dirname "$0")/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF

  mkdir -p $nfs_client_bundle
  source_parent="$(dirname "$source")"
  source_basename="$(basename "$source")"
  cd "\$source_parent"

  # do not quote \$source_basename to allow wildcards to be evaluated (needed in binaries.bash)
  tar -cf $nfs_client_bundle/bundle.tar \$source_basename

  #echo "chown to 1000:1000"
  chown -R 1000:1000 $nfs_client_bundle
  chmod -R go-rwxst,u=rwX,go=rX $nfs_client_bundle
EOF

show_elapsed $start_time

echo "extract files"
start_time=$(get_time)
ssh $ARTEFACTS_HOST tar -xf $nfs_server_bundle/bundle.tar -C $nfs_server_bundle
ssh $ARTEFACTS_HOST rm $nfs_server_bundle/bundle.tar
show_elapsed $start_time

echo "delete $source"
rm -rf $source
