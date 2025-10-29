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

# dada2-silva-reference_2025-08-14.tar.lz4 created from https://zenodo.org/records/14169026
# wget https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1
# wget https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1
# wget https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1
# mkdir dada2-silva dada2-silva-reference
# cd dada2-silva dada2-silva-reference
# mv 'silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1' silva_nr99_v138.2_toGenus_trainset.fa.gz
# mv 'silva_v138.2_assignSpecies.fa.gz?download=1' silva_v138.2_assignSpecies.fa.gz
# mv 'silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1' silva_nr99_v138.2_toSpecies_trainset.fa.gz
# gunzip silva_nr99_v138.2_toGenus_trainset.fa.gz 
# gunzip silva_nr99_v138.2_toSpecies_trainset.fa.gz 
# gunzip silva_v138.2_assignSpecies.fa.gz 
# cd ..
# tar -c dada2-silva-reference | lz4 > dada2-silva-reference.tar.lz4

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # silva 138.2
  f="dada2-silva-reference_2025-08-14.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f

  ls -lah $TOOLS_PATH/
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/dada2-silva-reference $JOB_NAME $BUILD_NUMBER
