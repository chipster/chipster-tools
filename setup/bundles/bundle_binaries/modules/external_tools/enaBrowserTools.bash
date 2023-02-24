##depends:none

source ../installation_files/functions.bash

cd ${TMPDIR_PATH}/
wget https://github.com/enasequence/enaBrowserTools/archive/7075a896f822e3ea3d3fac8bc10bcfeeb2506685.tar.gz
tar xzf 7075a896f822e3ea3d3fac8bc10bcfeeb2506685.tar.gz -C ${TOOLS_PATH}/
cd ${TOOLS_PATH}
ln -s enaBrowserTools-7075a896f822e3ea3d3fac8bc10bcfeeb2506685 enabrowsertools
