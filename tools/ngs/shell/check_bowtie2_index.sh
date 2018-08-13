#!/bin/bash 
#
#Automatic index cheking and indexing toold for bwa indexes
#NOTE! Does not do color space indexes
#K.M. 29.3. 2018



if [[ "$1" == "" ]]
then
   echo "check_bowtie2_index.sh command line syntax:"
   echo "  check_bowtie2_index.sh genome_file.fasta -bowtie2_path <path> -index_path <path> -compress -tar"
   exit 1
fi


#set defaults
compress=(0)
tar=(0)
bowtie2_path=("/opt/chipster/tools/bowtie2")

## No longer using the option to store precalcultead stuff to TMP
#set defaults
#if [ -d "/tmp" ]; then
#   mkdir -p "/tmp/bowtie2_indexes/tmp"
#   index_path=("/tmp/bowtie2_indexes/tmp")
#else
#   index_path=("./")
#fi

index_path=("./")

while [[ $# -ge 1 ]]
do
  case "$1" in
              '-bowtie2_path')
              bowtie2_path=$2
              shift
              shift
              ;; 
              '-index_path')
              index_path=$2
              shift
              shift
              ;; 
              '-compress')
              compress=(1)
              shift
              ;;
              '-tar')
              tar=(1)
              shift
              ;;
              *)
              genome=($1)
              shift
              ;;
  esac
done


#size=(`ls -l $genome | awk '{print $5}' `)
#checksum=(`md5sum $genome | awk '{print $1}'`)
location=$(pwd)

genome_file_type=$(file -b $genome | cut -d ' ' -f2)
if [[ $genome_file_type == "tar" ]]
then
  index_file_count=$(tar -tf $genome | grep -c -E ".bt2$")
  if [[ index_file_count -gt 4 ]]
  then
      gen_name=$(tar -tf $genome | grep -E ".rev.1.bt2$")
      gen_name=$(basename $gen_name .rev.1.bt2 )
      tar xvf $genome
      index_path=$(pwd)
      echo "The location of bowtie2_indexes:"
      echo "$index_path/$gen_name"
      exit 0
   else
      echo "The tar file does not contain bowtie2 index files"
      tar -tf $genome
      exit 1
   fi
fi 
  

#if [ ! -d $index_path/genome_0 ]; then
#  mkdir $index_path/genome_0
#fi


##look for matching size and md5sum
#genome_dir=(`grep -h $size $index_path/genome_*/size_and_md5 | grep $checksum | awk '{print $1}' | tail -1`)

###
#if [ ! $genome_dir == "" ]; then
#   echo "Pre-indexed genome found"
#else

echo "Calculating indexes"
cd $index_path
##   genome_dir=(`ls -d genome_* | awk -F "_" '{ print $NF }' | sort -n | tail -1 | awk '{ a = ( $1 + 1 ) }{print "genome_"a}'`)

genome_dir=("genome_dir")
mkdir $genome_dir
 cd $genome_dir
 cp "$location"/"$genome" ./$genome
 $bowtie2_path/bowtie2-build $genome $genome 

   #check that idexes are found
   echo "$genome"
   for f in 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2; do
       if [ -e "$genome"."$f" ]; then
         echo "$genome.$f OK"
       else
         echo "Indexing failed"
         echo "Index file $genome.$f not found"
         exit 1
       fi
   done
   if [[ tar -eq 1 ]]
   then 
     tar cvf bowtie2_index.tar $genome*  
     mv bowtie2_index.tar $location/bowtie2_index.tar
   fi

   cd $location

echo "The bowtie2_indexes are in directory:"
echo "$index_path/$genome_dir/$genome"
exit 

