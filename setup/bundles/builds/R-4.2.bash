#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

image="comp-20-04-r-deps"

# no tools-bin dependencies
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # variable f needs to be escaped to be evaluated only later on the k3s host
  f="R-4.2.0-single-cell_2023-01-05.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f
  f="R-4.2.0-phyloseq_2023-02-17.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f

  ls -lah $TOOLS_PATH/
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-4.2.0-single-cell $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-4.2.0-phyloseq $JOB_NAME $BUILD_NUMBER
