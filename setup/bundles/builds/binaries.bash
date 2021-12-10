#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

# we don't need a python image to install python
# most likely the tool wrappers are going to be written in R and hence this r-deps image will be used to run this python eventually
image="comp-20.04-r-deps"

# this installation doesn't need anythin from tools-bin
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # simply copy the old files (built in Ubuntu 16.04), because we cannot build these anymore in Ubuntu 20.04
  f="binaries-Ubuntu-16.04_2021-08-25.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f
  
  ls -lah $TOOLS_PATH/

  exit $?

  # used to be:
  # checkout https://github.com/chipster/chipster-tools.git
  # cd build && sudo bash run_install_chipster.bash
  
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/'*' $JOB_NAME $BUILD_NUMBER
