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

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF
  apt-get install -q -y unzip
EOF
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # FastQC, GPL v3 or later

  wget -nv https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
  unzip -q fastqc_v0.11.9.zip
  chmod a+x FastQC/fastqc
  mv FastQC fastqc-0.11.9
  mv fastqc-0.11.9 ${TOOLS_PATH}/
  rm fastqc_v0.11.9.zip
  
  cd $TOOLS_PATH
  ln -s fastqc-0.11.9 fastqc

  # Enabrowsertools

  cd $TMPDIR_PATH
  # download this specific commit, because the latest release version 1.5.4 from 2020 tries to use deprecated ssl v1
  wget https://github.com/enasequence/enaBrowserTools/archive/7075a896f822e3ea3d3fac8bc10bcfeeb2506685.tar.gz
  tar xzf 7075a896f822e3ea3d3fac8bc10bcfeeb2506685.tar.gz -C ${TOOLS_PATH}/
  cd ${TOOLS_PATH}
  ln -s enaBrowserTools-7075a896f822e3ea3d3fac8bc10bcfeeb2506685 enabrowsertools


  ls -lah $TOOLS_PATH/
  
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/fastqc-0.11.9 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/fastqc $JOB_NAME $BUILD_NUMBER
