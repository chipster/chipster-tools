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
  apt-get install -y build-essential
  apt-get install -y python3-dev
  apt-get install -y libssl-dev zlib1g-dev libncurses5-dev libncursesw5-dev libreadline-dev libsqlite3-dev
  apt-get install -y libgdbm-dev libdb5.3-dev libbz2-dev libexpat1-dev liblzma-dev libffi-dev uuid-dev

  # for RSeQC
  apt-get install -y gcc libz-dev
EOF
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  mkdir $TOOLS_PATH/python-3.8.11
  
  # pwd is a temp dir
  wget https://www.python.org/ftp/python/3.8.11/Python-3.8.11.tgz
  tar -xzf Python-3.8.11.tgz 
  cd Python-3.8.11
  ./configure --prefix=$TOOLS_PATH/python-3.8.11
  make
  make install

  $TOOLS_PATH/python-3.8.11/bin/python3 --version

  cd $TOOLS_PATH/python-3.8.11

  bin/pip3 install cutadapt
  bin/cutadapt --version

  bin/pip3 install multiqc
  mkdir ${TOOLS_PATH}/multiqc
  ln -s ../python-3.8.11/bin/multiqc ${TOOLS_PATH}/multiqc/multiqc

  bin/pip3 install RSeQC==5.0.1
  ln -s python-3.8.11/bin/ ${TOOLS_PATH}/rseqc

EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/python-3.8.11 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/multiqc $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/rseqc $JOB_NAME $BUILD_NUMBER
