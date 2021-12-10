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
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # simply copy the old files (built in Ubuntu 16.04), because we cannot build these anymore in Ubuntu 20.04
  f="python2-Ubuntu-16.04_2021-08-25.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f
  
  ls -lah $TOOLS_PATH/

EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/Python-2.7.12 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/htseq $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/macs $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/multiqc $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/rseqc $JOB_NAME $BUILD_NUMBER

exit $?

# Python
cd ${TMPDIR_PATH}/
wget -nv https://www.python.org/ftp/python/2.7.12/Python-2.7.12.tgz
tar xf Python-2.7.12.tgz
cd Python-2.7.12
mkdir ${TOOLS_PATH}/Python-2.7.12
./configure --prefix ${TOOLS_PATH}/Python-2.7.12
make
make install
PYTHON=Python-2.7.12

# install pip
  wget -nv https://bootstrap.pypa.io/pip/2.7/get-pip.py
  ${TOOLS_PATH}/${PYTHON}/bin/python get-pip.py

# Python tools
#
# S3
  ${TOOLS_PATH}/${PYTHON}/bin/pip install s3cmd

# HTSeq
  ${TOOLS_PATH}/${PYTHON}/bin/pip install matplotlib
  ${TOOLS_PATH}/${PYTHON}/bin/pip install HTSeq==0.6.1
  # htseq-count_chr is not part of distribution
  #wget -O ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr http://$NIC_MIRROR/pub/sci/molbio/chipster/dist/tools_extras/htseq/htseq-count_chr 
  cp ${TOOLS_PATH}/${PYTHON}/bin/htseq-count ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr
  sed -i 's/HTSeq.scripts.count/HTSeq.scripts.count_chr/' ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr
  chmod 755 ${TOOLS_PATH}/${PYTHON}/bin/htseq-count_chr
  wget -O ${TOOLS_PATH}/${PYTHON}/lib/python2.7/site-packages/HTSeq/scripts/count_chr.py http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/htseq/count_chr_v2.py
  # make links
  mkdir -p ${TOOLS_PATH}/htseq/
  ln -s ../${PYTHON}/bin/htseq-qa ${TOOLS_PATH}/htseq/htseq-qa
  ln -s ../${PYTHON}/bin/htseq-count ${TOOLS_PATH}/htseq/htseq-count
  ln -s ../${PYTHON}/bin/htseq-count_chr ${TOOLS_PATH}/htseq/htseq-count_chr
  
# MACS
  ${TOOLS_PATH}/${PYTHON}/bin/pip install MACS==1.4.2
  mkdir ${TOOLS_PATH}/macs/
  ln -s ../${PYTHON}/bin/macs14 ${TOOLS_PATH}/macs/macs14
  
# MACS2
  ${TOOLS_PATH}/${PYTHON}/bin/pip install numpy
  ${TOOLS_PATH}/${PYTHON}/bin/pip install MACS2==2.1.1.20160309
  ln -s ../${PYTHON}/bin/macs2 ${TOOLS_PATH}/macs/macs2
  
# MultiQC
  ${TOOLS_PATH}/${PYTHON}/bin/pip install multiqc
  mkdir ${TOOLS_PATH}/multiqc
  ln -s ../${PYTHON}/bin/multiqc ${TOOLS_PATH}/multiqc/multiqc
    
# RSeQC
  ${TOOLS_PATH}/${PYTHON}/bin/pip install RSeQC==2.6.4
  ln -s ${PYTHON}/bin ${TOOLS_PATH}/rseqc

# ZIFA
# Requires: scipy, pandas, numpy,  scikits.learn, matplotlib, pandas
# Some of these modules are already installed, but lets still list all the requirements for ZIFA
${TOOLS_PATH}/${PYTHON}/bin/pip install numpy scipy scikits.learn matplotlib
# pandas >0.21 needs currently some extra care
${TOOLS_PATH}/${PYTHON}/bin/pip install --no-build-isolation pandas

# ZIFA
cd ${TMPDIR_PATH}/
git clone https://github.com/epierson9/ZIFA
cd ZIFA
# TODO: on commandline this requires sudo
${TOOLS_PATH}/Python-2.7.12/bin/python setup.py install
cd ..
rm -rf ZIFA
