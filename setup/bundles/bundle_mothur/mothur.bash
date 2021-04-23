set -e

# mothur GPLv3
cd ${TMPDIR_PATH}/

# Retain version 1.41.3 as backup
cd ${TMPDIR_PATH}/
wget -nv https://github.com/mothur/mothur/releases/download/v1.41.3/Mothur.linux_64.zip
unzip -q Mothur.linux_64.zip
mv mothur ${TOOLS_PATH}/mothur-1.41.3
rm -rf  __MACOSX
#cd ${TOOLS_PATH}
#ln -s mothur-1.41.3 mothur

# Make version 1.44.3 the default
cd ${TMPDIR_PATH}/
wget -nv https://github.com/mothur/mothur/releases/download/v1.44.3/Mothur.linux.zip
unzip -q Mothur.linux.zip
mv mothur ${TOOLS_PATH}/mothur-1.44.3
rm -rf  __MACOSX
cd ${TOOLS_PATH}
ln -s mothur-1.44.3 mothur

mkdir -p ${TOOLS_PATH}/mothur-silva-reference/
curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v102.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/v132.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/
curl -s http://www.nic.funet.fi/pub/sci/molbio/chipster/dist/tools_extras/mothur/silva/silva-gold.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-silva-reference/

ln -s v132 ${TOOLS_PATH}/mothur-silva-reference/silva

mkdir -p ${TOOLS_PATH}/mothur-unite-reference/
curl -s https://a3s.fi/bundle-builds/mothur-UNITEv8_2020-12-15.tar.lz4  | lz4 -d | tar -x -C ${TOOLS_PATH}/mothur-unite-reference/
