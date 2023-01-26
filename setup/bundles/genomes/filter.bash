#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/genome-utils.bash

chipster login $CHIPSTER_SERVER -u $USERNAME -p $PASSWORD

echo "starting to filter genomes at $(date --iso-8601=seconds)"

download_dir="/mnt/artefacts/download_genomes/$DOWNLOAD_GENOMES_BUILD"

while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT
do
  echo "filter genome $SITE $SPECIES $VERSION $RELEASE $INDEX $DEFAULT at $(date --iso-8601=seconds)"
    
  session_name="$JOB_NAME-$BUILD_NUMBER-$SPECIES-$VERSION-$RELEASE"
  genome_dir="$SPECIES-$VERSION-$RELEASE"
  
  echo "find input files"
  pushd $download_dir/$genome_dir > /dev/null
  
  fasta="$(ls | grep .fa$ )"
  gtf="$(ls | grep .gtf$ )"
  
  chipster --quiet session create $session_name
  chipster --quiet session open $session_name
  
  echo "upload fasta"
    
  chipster --quiet dataset upload $fasta
        
  echo "filter fasta"
            
  chipster --quiet job run filter-fasta.py \
 	  -i input.fa=$fasta \
      --background
      
  popd > /dev/null
done < $download_dir/genomes.txt

wait_for_jobs $download_dir/genomes.txt

job_dir="/mnt/artefacts/$JOB_NAME/$BUILD_NUMBER"

mkdir -p $job_dir
pushd $job_dir

echo "starting to download results at $(date --iso-8601=seconds)"

while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT
do  
  session_name="$JOB_NAME-$BUILD_NUMBER-$SPECIES-$VERSION-$RELEASE"

  echo "download results from session $session_name"
  
  genome_dir="$SPECIES-$VERSION-$RELEASE"
  filtered_dir="$job_dir/$genome_dir"
  
  mkdir $filtered_dir
  pushd $filtered_dir > /dev/null
  
  chipster --quiet session open $session_name
  
  download_result_datasets filter-fasta.py "" "none"
  
  for file in $download_dir/$genome_dir/*.gtf; do
    ln -s $file $(basename $file)
  done
  
  popd > /dev/null
  
  if [ "$DELETE_SESSIONS" = "true" ]; then
    chipster --quiet session delete $session_name
  fi

done < $download_dir/genomes.txt

cp $download_dir/genomes.txt .

popd

echo "finished filtering at $(date --iso-8601=seconds)"
