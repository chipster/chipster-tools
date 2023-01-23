
function start_jobs {
  
  filter_genomes="$1"
  
  pushd $filter_genomes

  echo "starting to index genomes at $(date --iso-8601=seconds)"

  while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT
  do
    echo "index genome $SITE $SPECIES $VERSION $RELEASE $INDEX $DEFAULT at $(date --iso-8601=seconds)"
    
    genome_dir="$SPECIES-$VERSION-$RELEASE"
    session_name="$JOB_NAME-$BUILD_NUMBER-$SPECIES-$VERSION-$RELEASE"

    chipster --quiet session create $session_name
    chipster --quiet session open $session_name
    
    # let's use local file names to fill in the input names
    pushd $genome_dir
    
    fasta="$(ls | grep .fa$ || true)"
    gtf="$(ls | grep .gtf$ || true)"
    
    if [ -z $fasta ]; then
      echo "fasta not found in $genome_dir"
      exit 1
    fi
    
    if [ -z $gtf ]; then
      echo "gtf not found in $genome_dir"
      exit 1
    fi
          
    echo "upload $gtf"
    chipster --quiet dataset upload $gtf

    popd
          
    echo "create job: gtf-to-bed"
    chipster --quiet job run gtf-to-bed.py \
          -i input.gtf=$gtf \
          --background
    
    echo "create job: index-dexseq"
    chipster --quiet job run index-dexseq.py \
          -i input.gtf=$gtf \
          --background
          
    if [ "$INDEX" = "index-small" ] || [ "$INDEX" = "index-hisat2" ] || [ "$INDEX" = "index-all" ]; then
    
      # let's use local file names to fill in the input names
      pushd $genome_dir
    
      echo "upload $fasta"
      chipster --quiet dataset upload $fasta
      
      popd
      
    fi
          
    if [ "$INDEX" = "index-hisat2" ] || [ "$INDEX" = "index-all" ]; then
    
      if [ "$SPECIES" = "homo_sapiens" ] || [ "$SPECIES" = "mus_musculus" ] || [ "$SPECIES" = "canis_lupus_familiaris" ] || [ "$SPECIES" = "bos_taurus" ] || [ "$SPECIES" = "sus_scrofa" ] || [ "$SPECIES" = "ovis_aries" ] || [ "$SPECIES" = "felis_catus" ]; then

        chipster --quiet job run index-hisat2-large.py \
          -i input.fa=$fasta \
          -i input.gtf=$gtf \
          --background
      
      else
          
        chipster --quiet job run index-hisat2.py \
          -i input.fa=$fasta \
          -i input.gtf=$gtf \
          --background  
        
      fi
    fi
    
    if [ "$INDEX" = "index-all" ]; then
          
      chipster --quiet job run index-star.py \
          -i input.fa=$fasta \
          --background
    fi
    
    if [ "$INDEX" = "index-small" ] || [ "$INDEX" = "index-hisat2" ] || [ "$INDEX" = "index-all" ]; then
          
      echo "create jobs"
          
      chipster --quiet job run index-bwa.py \
          -i input.fa=$fasta \
          --background
          
  # Bowtie1 only for mirna indexes
  #    chipster --quiet job run index-bowtie.py \
  #   	    -i input.fa=$fasta \
  #        --background
          
      chipster --quiet job run index-bowtie2.py \
          -i input.fa=$fasta \
          --background
    fi
        
  done < $filter_genomes/genomes.txt
}

function download_results {

  filter_genomes="$1"

  mkdir -p genomes/bed
  mkdir -p genomes/dexseq
  mkdir -p genomes/gtf
  mkdir -p genomes/fasta

  while read -r SITE SPECIES VERSION RELEASE INDEX DEFAULT
  do  
    session_name="$JOB_NAME-$BUILD_NUMBER-$SPECIES-$VERSION-$RELEASE"
    genome_dir="$SPECIES-$VERSION-$RELEASE"
    
    echo "find gtf and fasta"
    gtf=$(ls $filter_genomes/$genome_dir | grep .gtf$ )
    fasta=$(ls $filter_genomes/$genome_dir | grep .fa$ )
    gtf_basename=$(basename $gtf .gtf)
    fasta_basename=$(basename $fasta .fa)

    echo "copy gtf and fasta"
    cp $filter_genomes/$genome_dir/$gtf $job_dir/genomes/gtf/
    cp $filter_genomes/$genome_dir/$fasta $job_dir/genomes/fasta/
    
    add_file_to_package genomes/gtf $gtf $gtf_basename
    add_file_to_package genomes/fasta $fasta $fasta_basename
    
    chipster --quiet session open $session_name
    
    echo "download bed file"  
    download_result_datasets gtf-to-bed.py genomes/bed $gtf_basename  
    
    echo "download dexseq index"
    download_result_datasets index-dexseq.py genomes/dexseq $gtf_basename
    
    if [ $INDEX == "no-index" ]; then
      echo "no other indexes to download"
    fi
      
    if [ "$INDEX" = "index-small" ] || [ "$INDEX" = "index-hisat2" ] || [ "$INDEX" = "index-all" ]; then
      echo "dowload bwa index"
      download_result_datasets index-bwa.py genomes/indexes/bwa $fasta_basename
      # download_result_datasets index-bowtie.py genomes/indexes/bowtie $fasta_basename

      echo "dowload bowtie2 index"
      download_result_datasets index-bowtie2.py genomes/indexes/bowtie2 $fasta_basename
      
      echo "symlink fasta for bwa"
      echo create_symlink ../../fasta/$fasta genomes/indexes/bwa $fasta $fasta_basename
      create_symlink ../../fasta/$fasta genomes/indexes/bwa $fasta $fasta_basename
      # create_symlink ../../fasta/$fasta genomes/indexes/bowtie $fasta $fasta_basename

      echo "symlink fasta for bowtie2"
      create_symlink ../../fasta/$fasta genomes/indexes/bowtie2 $fasta_basename.fa $fasta_basename
    fi
    
    if [ "$INDEX" = "index-hisat2" ] || [ "$INDEX" = "index-all" ]; then

      echo "dowload hisat2 index"
      if [ "$SPECIES" = "homo_sapiens" ] || [ "$SPECIES" = "mus_musculus" ] || [ "$SPECIES" = "canis_lupus_familiaris" ] || [ "$SPECIES" = "bos_taurus" ] || [ "$SPECIES" = "sus_scrofa" ] || [ "$SPECIES" = "ovis_aries" ] || [ "$SPECIES" = "felis_catus" ]; then
        download_result_datasets index-hisat2-large.py genomes/indexes/hisat2 $gtf_basename
      else
        download_result_datasets index-hisat2.py genomes/indexes/hisat2 $gtf_basename
      fi
      
      create_symlink ../../fasta/$fasta genomes/indexes/hisat2 $gtf_basename.fa $gtf_basename
    fi
    
    if [ "$INDEX" = "index-all" ]; then
    
      mkdir -p genomes/indexes/star/$fasta_basename
      
      echo "download star results"
      # no need for separate package name, because Star has fasta name already in the directory path
      download_result_datasets index-star.py genomes/indexes/star/$fasta_basename ""
      
      echo "symlink fasta for star"
      create_symlink ../../../fasta/$fasta genomes/indexes/star/$fasta_basename $fasta $fasta_basename
    fi
    
    if [ "$DELETE_SESSIONS" = "true" ]; then
      chipster --quiet session delete $session_name
    fi
    
    if [ "$DEFAULT" = "default" ]; then

      create_symlink $gtf genomes/gtf default $gtf_basename
      create_symlink $gtf_basename.DEXSeq.gtf genomes/dexseq default $gtf_basename
      create_symlink $gtf_basename.bed genomes/bed default $gtf_basename
      create_symlink $fasta genomes/fasta default $fasta_basename
      
      create_symlink $fasta genomes/indexes/bwa default $fasta_basename
      create_symlink $gtf_basename.fa genomes/indexes/bowtie2 default $fasta_basename
      create_symlink $gtf_basename.fa genomes/indexes/hisat2 default $gtf_basename
      create_symlink $fasta_basename genomes/indexes/star default ""
    fi

  done < $filter_genomes/genomes.txt

  echo "delete Chipster CLI logs"

  rm -rf logs/ */logs/ */*/logs/ */*/*/logs/ */*/*/*/logs/

  echo "Wait child processes to end. There shouldn't be any. If the job stops here, some Chipster CLI client process is stuck. "
  wait
}

function download_mirna {
  echo "download miRNA"

  wget http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/miRNA_mappings/mirna_genomes_mirbase-22_flat.tar.gz

  tar -xf mirna_genomes_mirbase-22_flat.tar.gz -C genomes

  # create separate packages for each mirna index like we have for the other genomes

  for species in Homo_sapiens Mus_musculus Rattus_norvegicus; do
    for aligner in bowtie bowtie2 bwa; do
      for file in $(tar -tf mirna_genomes_mirbase-22_flat.tar.gz | grep indexes/$aligner/${species}_mirna); do 
        add_file_to_package genomes $file indexes_${aligner}_${species}_mirna
      done
    done
    add_file_to_package genomes/fasta ${species}_mirna.fa ${species}_mirna
  done

  # Bowtie has only mirna genomes, so let's use that as a default
  create_symlink Homo_sapiens_mirna.fa genomes/indexes/bowtie default Homo_sapiens_mirna

  rm mirna_genomes_mirbase-22_flat.tar.gz
}