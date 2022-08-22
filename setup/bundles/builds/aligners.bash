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
  apt-get install -q -y unzip make build-essential libz-dev gcc zlib1g-dev
EOF

  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # Hisat2

  # current working directory is a temp dir
  wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
  unzip download
  rm download
  cd hisat2-2.2.1
  cd ..
  mv hisat2-2.2.1 $TOOLS_PATH/
  cd $TOOLS_PATH
  ln -s hisat2-2.2.1 hisat2

  # STAR
  
  cd $TMPDIR_PATH
  wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
  tar -xzf 2.7.10a.tar.gz
  rm 2.7.10a.tar.gz 
  cd STAR-2.7.10a/source
  make STAR
  cd ../..
  mv STAR-2.7.10a STAR-2.7.10a_build
  mkdir STAR-2.7.10a
  mv STAR-2.7.10a_build/bin/Linux_x86_64_static/STAR STAR-2.7.10a/
  rm -rf STAR-2.7.10a_build
  mv STAR-2.7.10a $TOOLS_PATH/
  cd $TOOLS_PATH
  ln -s STAR-2.7.10a STAR

  # Bowtie2

  cd $TMPDIR_PATH
  wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.5/bowtie2-2.4.5-linux-x86_64.zip/download
  unzip download
  rm download
  mv bowtie2-2.4.5-linux-x86_64 $TOOLS_PATH/
  cd $TOOLS_PATH
  ln -s bowtie2-2.4.5-linux-x86_64 bowtie2

  # Bowtie

  cd $TMPDIR_PATH
  wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip/download
  unzip download
  rm download
  mv bowtie-1.3.1-linux-x86_64 $TOOLS_PATH
  cd $TOOLS_PATH
  ln -s bowtie-1.3.1-linux-x86_64 bowtie

  # BWA
  
  cd $TMPDIR_PATH
  wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download
  tar -xvf download
  rm download
  cd bwa-0.7.17/
  make
  cd ..
  mv bwa-0.7.17 $TOOLS_PATH/
  cd $TOOLS_PATH
  ln -s bwa-0.7.17 bwa
EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/hisat2-2.2.1 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/hisat2 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/STAR-2.7.10a $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/STAR $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bowtie2-2.4.5-linux-x86_64 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bowtie2 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bowtie-1.3.1-linux-x86_64 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bowtie $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bwa-0.7.17 $JOB_NAME $BUILD_NUMBER
bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/bwa $JOB_NAME $BUILD_NUMBER
