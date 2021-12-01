#!/bin/bash

set -e

source $OPENRC_PATH

# wrapper script is going to be written in R
image="comp-20.04-r-deps"

# this needs python from tools-bin
BUNDLE_COLLECTION_VERSION="58"

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image $BUNDLE_COLLECTION_VERSION"

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root -" <<EOF

  apt-get update
  apt-get install -y libssl-dev
  apt-get install -y python3-venv
  apt-get install -y build-essential
  apt-get install -y python3-dev
EOF
  
ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu -" <<EOF

  /opt/chipster/tools/python-3.8.11/bin/python3 -m venv /opt/chipster/tools/umi-tools/venv
  source /opt/chipster/tools/umi-tools/venv/bin/activate
  
  pip3 install wheel
  pip3 install umi_tools
  
  umi_tools --help
EOF

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/move-to-artefacts.bash /opt/chipster/tools/umi-tools $JOB_NAME $BUILD_NUMBER"

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER"
