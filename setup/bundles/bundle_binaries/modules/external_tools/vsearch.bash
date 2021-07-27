##depends:none

source ../installation_files/functions.bash

cd ${TMPDIR_PATH}/
wget https://github.com/torognes/vsearch/releases/download/v2.17.1/vsearch-2.17.1-linux-x86_64
chmod u+x vsearch-2.17.1-linux-x86_64
mkdir vsearch-2.17.1
mv vsearch-2.17.1-linux-x86_64 vsearch-2.17.1
ln -s vsearch-2.17.1 vsearch

cd vsearch
ln -s vsearch-2.17.1-linux-x86_64 vsearch
cd ..

mv vsearch-2.17.1 vsearch ${TOOLS_PATH}/
