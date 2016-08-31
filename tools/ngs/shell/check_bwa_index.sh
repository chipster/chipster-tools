#!/bin/bash 
#
#Automatic index cheking and indexing tool for bwa indexes
#NOTE! Does not do color space indexes
#K.M. 29.8. 2016


if [[ "$1" == "" ]]
then
   echo "check_bwa_index.sh command line syntax:"
   echo "  check_bwa_index.sh genome_file.fasta -bwa_path <path> -index_path <path> -compress"
   exit 1
fi

#set defaults
compress=(0)
bwa_path=("/opt/chipster/tools/bwa")
if [ -d "/tmp" ]; then
  mkdir -p "/tmp/bwa_indexes/tmp"
   index_path=("/tmp/bwa_indexes/tmp")
else
   index_path=("./")
fi


while [[ $# -ge 1 ]]
do
  case "$1" in
              '-bwa_path')
              bwa_path=$2
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
              *)
              genome=($1)
              shift
              ;;
       esac
done



size=(`ls -l $genome | awk '{print $5}' `)
checksum=(`md5sum $genome | awk '{print $1}'`)
location=(`pwd`)

if [ ! -d $index_path/genome_0 ]; then
  mkdir $index_path/genome_0
fi


#look for matching size and md5sum
genome_dir=(`grep -h $size $index_path/genome_*/size_and_md5 | grep $checksum | awk '{print $1}' | tail -1`)


#
if [ ! $genome_dir == "" ]; then
   echo "Pre-indexed genome found"
else
   echo "Calculating indexes"
   cd $index_path
   genome_dir=(`ls -d genome_* | awk -F "_" '{ print $NF }' | sort -n | tail -1 | awk '{ a = ( $1 + 1 ) }{print "genome_"a}'`)
   mkdir $genome_dir
   cd $genome_dir
   cp "$location"/"$genome" ./$genome
  
   size_mb=(` expr $size / 1000000 ` )     
   #Choose the indextype based on the genome size
   if [ $size_mb -gt 2000 ] 
   then
     indextype=("bwtsw")
   else
     indextype=("is")
   fi
   $bwa_path/bwa index -a $indextype $genome 
   
   #check that idexes are found
   for f in ann amb bwt pac sa ; do
       echo "$genome.$f"
       if [ -e "$genome"."$f" ]; then
         echo "$genome.$f OK"
       else
         echo "Indexing failed"
         echo "Index file $genome.$f not found"
         exit 1
       fi
   done
   if [[ compress -eq 1 ]]
   then 
      tar zcvf ${genome}_bwa_index.tgz $genome*  
   fi
   echo "$genome_dir $size $checksum" > size_and_md5
   cd $location
fi

echo "The bwa_indexes are in directory:"
echo "$index_path/$genome_dir"
exit 





