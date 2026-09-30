#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

image="comp-20-04-r-deps"
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

  mkdir -p ${TOOLS_PATH}/mothur-silva-reference/

  curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/silva-gold.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
  mkdir ${TOOLS_PATH}/mothur-silva-reference/v138.2
  curl -s https://a3s.fi/bundle-builds/mothur-silva.nr_v138_2_2025-03-17.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/v138.2


  # downloaded from https://unite.ut.ee/repository.php (two small "Fungi" packages from the table, not the large "download1" etc. packages at the top)

  mkdir -p ${TOOLS_PATH}/mothur-unite-reference/  
  curl -s https://a3s.fi/bundle-builds/mothur-UNITEv10_2025-03-17.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-unite-reference/

  rm ${TOOLS_PATH}/mothur-unite-reference/UNITEv10_sh_97.fasta
  rm ${TOOLS_PATH}/mothur-unite-reference/UNITEv10_sh_97.tax
  rm ${TOOLS_PATH}/mothur-unite-reference/UNITEv10_sh_97_s.fasta
  rm ${TOOLS_PATH}/mothur-unite-reference/UNITEv10_sh_97_s.tax
  ls -lah
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur-silva-reference $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/mothur-unite-reference $JOB_NAME $BUILD_NUMBER
