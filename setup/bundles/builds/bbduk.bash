#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

image="comp-20-04-r-deps"
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # current working directory is a temp dir
  wget https://sourceforge.net/projects/bbmap/files/BBMap_38.94.tar.gz/download -O BBMap_38.94.tar.gz
  tar -xzf BBMap_38.94.tar.gz
  rm BBMap_38.94.tar.gz 
  mv bbmap $TOOLS_PATH/BBMap_38.94
  cd $TOOLS_PATH
  ln -s BBMap_38.94 bbmap

  # test
  bbmap/stats.sh
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/BBMap_38.94 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bbmap $JOB_NAME $BUILD_NUMBER
