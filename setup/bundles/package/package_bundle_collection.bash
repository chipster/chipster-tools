#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/../builds/build-env.bash

image="comp-20.04-r-deps"

if [ -z "$BUNDLE_COLLECTION_VERSION" ]; then
  echo "$BUNDLE_COLLECTION_VERSION not set"
  exit 1
fi

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

# no need to mount the whole bundle collection, because we want to package bundles one by one
bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"\"

bundle_collection_dir=/mnt/artefacts/collect_bundles/$BUNDLE_COLLECTION_VERSION

# bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER root - <<EOF
#     apt-get install -y liblz4-tool
# EOF

build_dir=/mnt/artefacts/$JOB_NAME/$BUILD_NUMBER
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

    mkdir -p $build_dir/parts

    pushd $bundle_collection_dir

    for bundle in *; do
        cd \$bundle

        # .tar.lz4 parts

        echo "\$(date --rfc-3339=seconds) Creating tar.lz4 parts of \$bundle"

        # list recursively all files and symlinks and split it to chunks of 1k lines

        # add large and small files to separate packages
        # this allows small files to be extracted in parallel (which is otherwise super slow on shared file systems)
        # and limits the size of the largest packges
        # store the file lists in build_dir to hide them from the following find commands  
        find . -type f -size +1M > $build_dir/parts/\${bundle}_tools_part_large  
        find . -type f -size 1M > $build_dir/parts/\${bundle}_tools_part_small
        find . -type f -size -1M >> $build_dir/parts/\${bundle}_tools_part_small
        find . -type l >> $build_dir/parts/\${bundle}_tools_part_small

        cat $build_dir/parts/\${bundle}_tools_part_small | split -a 3 -l 10000 -d - \${bundle}_small_
        cat $build_dir/parts/\${bundle}_tools_part_large | split -a 3 -l 100 -d - \${bundle}_large_
        rm $build_dir/parts/\${bundle}_tools_part_small
        rm $build_dir/parts/\${bundle}_tools_part_large  

        # create a tar.lz4 package of the files in each chunk
        # ls: 		list the file names of the file list chunks
        # parallel:	execute 8 following jobs ins parallel
        # tar:		create tar archive of the files in the chunk
        # lz4c:		compress it with lz4
        # tee:		create copy of the stream
        # md5sum: 	calculate a md5 of the stream
        # sed:		replace the '-' in the md5sum output with the real filename
        # >:			redirect the compressed stream to a file in $build_dir
        ls \${bundle}_tools_part_* | parallel -j 4 "tar -c -T {} | lz4c | tee >(md5sum | sed s/-/{}.tar.lz4/ > $build_dir/parts/{}.tar.lz4.md5) > $build_dir/parts/{}.tar.lz4"  
        rm \${bundle}_tools_part_*        

        cd ..
    done
    popd

    pushd $build_dir/parts
    ls > files.txt
    popd

    echo "\$(date --rfc-3339=seconds) Done"

    sync

EOF

