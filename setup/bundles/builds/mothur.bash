#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

image="comp-20.04-r-deps"
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

  # mothur GPLv3

  # Retain version 1.41.3 as backup
  cd ${TMPDIR_PATH}/
  wget -nv https://github.com/mothur/mothur/releases/download/v1.41.3/Mothur.linux_64.zip
  unzip -q Mothur.linux_64.zip
  mv mothur ${TOOLS_PATH}/mothur-1.41.3
  rm -rf  __MACOSX
  #cd ${TOOLS_PATH}
  #ln -s mothur-1.41.3 mothur

  # Make version 1.44.3 the default
  cd ${TMPDIR_PATH}/
  wget -nv https://github.com/mothur/mothur/releases/download/v1.44.3/Mothur.linux.zip
  unzip -q Mothur.linux.zip
  mv mothur ${TOOLS_PATH}/mothur-1.44.3
  rm -rf  __MACOSX
  cd ${TOOLS_PATH}
  ln -s mothur-1.44.3 mothur

  mkdir -p ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v102.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/silva-gold.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  curl -s https://a3s.fi/bundle-builds/mothur-silva-v138.1_2021-04-26.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/

  ln -s v138.1 ${TOOLS_PATH}/mothur-silva-reference/silva

  mkdir -p ${TOOLS_PATH}/mothur-unite-reference/
  curl -s https://a3s.fi/bundle-builds/mothur-UNITEv8_2020-12-15.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-unite-reference/

  ls -lah
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur-1.41.3 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur-1.44.3 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur-silva-reference $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur-unite-reference $JOB_NAME $BUILD_NUMBER
