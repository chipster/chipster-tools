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

  # variable f needs to be escaped to be evaluated only later on the k3s host
  f="R-4.1.1_2021-09-10.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $tools_dir; rm \$f
  f="R-4.1.0-single-cell_2021-09-27.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $tools_dir; rm \$f
  f="R-4.1.1-statistics_2021-12-02.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $tools_dir; rm \$f

  ls -lah $tools_dir
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash /opt/chipster/tools/R-4.1.1 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash /opt/chipster/tools/R-4.1.0-single-cell $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash /opt/chipster/tools/R-4.1.1-statistics $JOB_NAME $BUILD_NUMBER