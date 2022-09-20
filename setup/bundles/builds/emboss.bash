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

# bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF
#   # emboss needs this to create png images
#   # TODO we would probably need this (or without -dev) in comp-20.04-r-deps too
#   # TODO this was installed in Ubuntu 16.04, but does not exists in ubuntu 20.04. Would libbd-dev work?
#   apt-get -y install libgd2-noxpm-dev 
# EOF

# there is the same Emboss version for Ubuntu 16.04 already
EMBOSS_PATH=${TOOLS_PATH}/EMBOSS-6.5.7-20.04
EMBOSS_OPTIONS="--prefix=${EMBOSS_PATH} --without-x"
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # current working directory is a temp dir

  # EMBOSS, GPL  
 	 	 
  wget ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz

  tar -xzvf EMBOSS-6.5.7.tar.gz

  cd EMBOSS-6.5.7

  ./configure ${EMBOSS_OPTIONS}
  make
  make install
	
  cd $TOOLS_PATH

  # there is already a symlink "emboss" which points to same Emboss version compiled in old Ubuntu 16.04
  ln -s EMBOSS-6.5.7-20.04 ${TOOLS_PATH}/emboss-20.04

  # EMBOSS extras

  cd $TMPDIR_PATH
  wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/MEME-4.7.650.tar.gz
  tar -xzvf MEME-4.7.650.tar.gz
  cd MEME-4.7.650
  ./configure ${EMBOSS_OPTIONS}
  make
  make install
  cd ..

  # REBASE reference data and indeces	
  cd ${EMBOSS_PATH}/share/EMBOSS/data/REBASE
  #wget ftp://ftp.neb.com/pub/rebase/withrefm.txt
  #wget ftp://ftp.neb.com/pub/rebase/proto.txt
  wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/withrefm.txt
  wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/EMBOSS/proto.txt

  ../../../../bin/rebaseextract -infile withrefm.txt -protofile proto.txt -equivalences Y

EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/EMBOSS-6.5.7-20.04 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/emboss-20.04 $JOB_NAME $BUILD_NUMBER
