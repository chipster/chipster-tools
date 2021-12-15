#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

image="comp-20.04-r-deps"

# no tools-bin dependencies
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # variable f needs to be escaped to be evaluated only later on the k3s host
  # recent packages are downloaded directly from Allas
  f="R-3.6.1-phyloseq_2020-11-04.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f
  f="R-3.6.1-single-cell_2020-10-14.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f

  # old packages had their primary copy in nic
  lz4 -d /mnt/artefacts/downloads/R-3.4.3_ubuntu-16.04_2018-08-29.tar.lz4 | tar x -C $TOOLS_PATH
  lz4 -d /mnt/artefacts/downloads/R-3.6.1_ubuntu-16.04_2019-07-31.tar.lz4 | tar x -C $TOOLS_PATH
  lz4 -d /mnt/artefacts/downloads/R-3.6.1-plain_ubuntu-16.04_2019-07-31.tar.lz4 | tar x -C $TOOLS_PATH

  ls -lah $TOOLS_PATH/
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-3.6.1-phyloseq $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-3.6.1-single-cell $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-3.4.3 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-3.6.1 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/R-3.6.1-plain $JOB_NAME $BUILD_NUMBER
