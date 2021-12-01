#!/bin/bash

set -e

source $OPENRC_PATH

image="comp-20.04-r-deps"
BUNDLE_COLLECTION_VERSION=""

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\""
  
ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu -" <<EOF

  # current working directory is a temp dir
  wget https://sourceforge.net/projects/bbmap/files/BBMap_38.94.tar.gz/download -O BBMap_38.94.tar.gz
  tar -xzf BBMap_38.94.tar.gz
  rm BBMap_38.94.tar.gz 
  mv bbmap /opt/chipster/tools/BBMap_38.94
  cd /opt/chipster/tools
  ln -s BBMap_38.94 bbmap

  # test
  bbmap/stats.sh
EOF

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/move-to-artefacts.bash /opt/chipster/tools/BBMap_38.94 $JOB_NAME $BUILD_NUMBER"
ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/move-to-artefacts.bash /opt/chipster/tools/bbmap $JOB_NAME $BUILD_NUMBER"

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER"
