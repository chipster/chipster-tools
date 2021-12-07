#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

image="comp-20.04-r-deps"

# no tools-bin dependencies
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"

tools_dir="/opt/chipster/tools"

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  wget https://a3s.fi/bundle-builds/R-3.3.2_2020-10-15.tar.lz4
  lz4 -d R-3.3.2_2020-10-15.tar.lz4 -c | tar x -C /opt/chipster/tools
  rm R-3.3.2_2020-10-15.tar.lz4

  ls -lah $tools_dir
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash /opt/chipster/tools/R-3.3.2 $JOB_NAME $BUILD_NUMBER
