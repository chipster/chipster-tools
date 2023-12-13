#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

# we don't need a python image to install python
# most likely the tool wrappers are going to be written in R and hence this r-deps image will be used to run this python eventually
image="comp-20-04-r-deps"

# this installation doesn't need anythin from tools-bin
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF

  apt-get update
  apt-get install -y libssl-dev libffi-dev libbz2-dev liblzma-dev

  echo "TOOLS_PATH: $TOOLS_PATH"

  # tools-bin used to be mounted in /opt/chipster/tools in old runtimes, but when R is in image, tools-bin is 
  # mounted in /opt/chipster/tools-bin instead. cellbender will write the python path in the shebang line, so 
  # the python must be in the final path when cellbender is installed (or the shebang line can be fixed afterwards).
  cd /opt/chipster
  ln -s $TOOLS_PATH tools-bin
  ls -lah

EOF
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF  
  
  # pwd is a temp dir
  # cellbender works only in Python 3.7
  wget https://www.python.org/ftp/python/3.7.17/Python-3.7.17.tgz
  tar -xzf Python-3.7.17.tgz 
  cd Python-3.7.17
  ./configure --prefix=/opt/chipster/tools-bin/python-3.7.17
  make
  make install

  /opt/chipster/tools-bin/python-3.7.17/bin/python3 --version

  cd /opt/chipster/tools-bin/python-3.7.17

  bin/pip3 install cellbender==v0.3.0
  bin/cellbender --version

EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/python-3.7.17 $JOB_NAME $BUILD_NUMBER
