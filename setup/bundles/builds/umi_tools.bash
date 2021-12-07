#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# wrapper script is going to be written in R
image="comp-20.04-r-deps"

# this needs python from tools-bin
BUNDLE_COLLECTION_VERSION="58"

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image $BUNDLE_COLLECTION_VERSION

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF

  apt-get update
  apt-get install -y libssl-dev
  apt-get install -y python3-venv
  apt-get install -y build-essential
  apt-get install -y python3-dev
EOF
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  /opt/chipster/tools/python-3.8.11/bin/python3 -m venv /opt/chipster/tools/umi-tools/venv
  source /opt/chipster/tools/umi-tools/venv/bin/activate
  
  pip3 install wheel
  pip3 install umi_tools
  
  umi_tools --help
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash /opt/chipster/tools/umi-tools $JOB_NAME $BUILD_NUMBER
