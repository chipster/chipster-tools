#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

chipster login $CHIPSTER_SERVER -u $USERNAME -p $PASSWORD

chipster session create download-genomes-version-check-$BUILD_NUMBER
chipster session open download-genomes-version-check-$BUILD_NUMBER

echo "starting to check versions at $(date --iso-8601=seconds)"

while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT
do
   echo "check version of $SITE $SPECIES $VERSION $RELEASE"
   chipster --quiet job run genome-download.py \
   	-p species=$SPECIES \
    -p version=$VERSION \
    -p release=$RELEASE \
    -p site=$SITE \
    -p assembly=primary_assembly \
    -p action=check_version_only
done < ~/git/chipster-tools/setup/bundles/genomes/$GENOMES

chipster session delete download-genomes-version-check-$BUILD_NUMBER

job_dir="/mnt/artefacts/$JOB_NAME/$BUILD_NUMBER"
mkdir -p $job_dir
pushd $job_dir

echo "starting to download genomes at $(date --iso-8601=seconds)"

while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT
do
   session_name="download-genome-$SPECIES-$VERSION-$RELEASE-$BUILD_NUMBER"
   chipster --quiet session create $session_name
   chipster --quiet session open $session_name
   echo "download from Ensembl to Chipster $SITE $SPECIES $VERSION $RELEASE"
   chipster --quiet job run genome-download.py \
   	-p species=$SPECIES \
    -p version=$VERSION \
    -p release=$RELEASE \
    -p site=$SITE \
    -p assembly=primary_assembly \
    -p action=download
   
   genome_dir="$SPECIES-$VERSION-$RELEASE"
   mkdir $genome_dir
   pushd $genome_dir
   
   echo "download from Chipster to artefacts"
   datasets="$(chipster -o json dataset list)"
   for j in $(seq $(echo "$datasets" | jq length)); do     
     dataset_name="$(echo $datasets | jq .[$j-1].name -r)"
     dataset_id="$(echo $datasets | jq .[$j-1].datasetId -r)"
     
     chipster --quiet dataset download --file $dataset_name $dataset_id
   done
   
   chipster --quiet session delete $session_name
   popd
   
done < ~/git/chipster-tools/setup/bundles/genomes/$GENOMES

cp ~/git/chipster-tools/setup/bundles/genomes/$GENOMES genomes.txt

echo "finished downloads at $(date --iso-8601=seconds)"
