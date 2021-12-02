#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $OPENRC_PATH

# we don't need a python image to install python
# most likely the tool wrappers are going to be written in R and hence this r-deps image will be used to run this python eventually
image="comp-20.04-r-deps"

# this installation doesn't need anythin from tools-bin
BUNDLE_COLLECTION_VERSION=""

function finish {
  ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER"
}
trap finish EXIT

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\""

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root -" <<EOF

  apt-get update
  apt-get install -y build-essential
  apt-get install -y python3-dev
  apt-get install -y libssl-dev zlib1g-dev libncurses5-dev libncursesw5-dev libreadline-dev libsqlite3-dev
  apt-get install -y libgdbm-dev libdb5.3-dev libbz2-dev libexpat1-dev liblzma-dev libffi-dev uuid-dev
EOF
  
ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu -" <<EOF

  mkdir /opt/chipster/tools/python-3.8.11
  
  # pwd is a temp dir
  wget https://www.python.org/ftp/python/3.8.11/Python-3.8.11.tgz
  tar -xzf Python-3.8.11.tgz 
  cd Python-3.8.11
  ./configure --prefix=/opt/chipster/tools/python-3.8.11
  make
  make install

  /opt/chipster/tools/python-3.8.11/bin/python3 --version
EOF

ssh $K3S_BUILD_HOST "bash $K3S_BUNDLE_DIR/move-to-artefacts.bash /opt/chipster/tools/python-3.8.11 $JOB_NAME $BUILD_NUMBER"
