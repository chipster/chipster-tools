#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/../builds/build-env.bash

# any image would do, now used only for chown
image="base"

# this installation doesn't need anything from tools-bin
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"

# bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF
#   apt-get install -q -y unzip
# EOF

# run this locally on the Jenkins server for now, because moving to k3s would require more work (how to transfer scripts and genomes.txt, and allow network access to the Chipster VM)

source $(dirname "$0")/genome-utils.bash
source $(dirname "$0")/index-functions.bash

filter_genomes="/mnt/artefacts/filter_genomes/$FILTER_GENOMES_BUILD"

chipster login $CHIPSTER_SERVER -u $USERNAME -p $PASSWORD

start_jobs $filter_genomes

wait_for_jobs $filter_genomes/genomes.txt

echo "starting to download results at $(date --iso-8601=seconds)"

job_dir="/mnt/artefacts/$JOB_NAME/$BUILD_NUMBER"

mkdir -p $job_dir
pushd $job_dir

packages_file=$job_dir/genomes/packages.tsv

download_results $filter_genomes
download_mirna

popd

echo "finished indexing at $(date --iso-8601=seconds)"

bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF
  echo "chown files"
  # this is the job_dir path, but I don't want to chown everything by accident
  chown -R 1000:1000 /mnt/artefacts/index_genomes/$BUILD_NUMBER
EOF
