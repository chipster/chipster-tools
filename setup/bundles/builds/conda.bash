#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

source $(dirname "$0")/build-env.bash

# we don't need a python image to install python
# most likely the tool wrappers are going to be written in R and hence this r-deps image will be used to run this python eventually
image="comp-20.04-r-deps"

# this installation doesn't need anythin from tools-bin
BUNDLE_COLLECTION_VERSION=""

function finish {
  bash $BUNDLE_SCRIPTS_DIR/clean-up.bash $JOB_NAME $BUILD_NUMBER
}
trap finish EXIT

bash $BUNDLE_SCRIPTS_DIR/start-pod.bash $JOB_NAME $BUILD_NUMBER $image \"$BUNDLE_COLLECTION_VERSION\"
  
bash $BUNDLE_SCRIPTS_DIR/run-in-pod.bash $JOB_NAME $BUILD_NUMBER ubuntu - <<EOF

  # simply copy the old files (built in Ubuntu 16.04), because we cannot build these anymore in Ubuntu 20.04
  f="conda-Ubuntu-16.04_2021-08-25.tar.lz4"; wget https://a3s.fi/bundle-builds/\$f; lz4 -d \$f -c | tar x -C $TOOLS_PATH; rm \$f
  
  ls -lah $TOOLS_PATH/

EOF

bash $BUNDLE_SCRIPTS_DIR/move-to-artefacts.bash $TOOLS_PATH/miniconda3 $JOB_NAME $BUILD_NUMBER

exit $?

# old build build setup for reference when building new tools in Ubuntu 20.04

#Miniconda3 installation:

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
TOOLS_PATH="/opt/chipster/tools"
CONDA_PATH="/opt/chipster/tools/miniconda3"
sudo mkdir -p $TOOLS_PATH
sudo chown $(whoami) $TOOLS_PATH
bash Miniconda3-latest-Linux-x86_64.sh -b -p $CONDA_PATH
export PATH=${PATH}:${CONDA_PATH}/bin

# update conda to hide warnings
conda update -n base -c defaults conda -y

conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

mv $HOME/.condarc $CONDA_PATH
conda create -n chipster_tools -y

cat << EOF > $CONDA_PATH/conda_execute
#!/bin/bash

conda_path=("$CONDA_PATH")

c_env=\$(echo \$1 | awk -F "/" '{print \$1}')
c_tool=\$(echo \$1 | awk -F "/" '{print \$2}')
shift
export PATH=\${PATH}:\${conda_path}/bin
conda activate \$c_env
\$c_tool \$@
EOF

chmod u+x $CONDA_PATH/conda_execute


# Activate conda
CONDA_PATH=("/opt/chipster/tools/miniconda3")
export PATH=${PATH}:${CONDA_PATH}/bin

# fix conda activate for scripts https://github.com/conda/conda/issues/7980
source $(conda info --base)/etc/profile.d/conda.sh

conda activate chipster_tools

# minimap
conda install minimap2=2.9 -y

# weasyprint
conda install -c jlmenut weasyprint -y
pip install -U html5lib=="0.999999999"

# umap for Seurat
# umap-learn refuses to install with the default python version 3.9
conda install -c conda-forge umap-learn python=3.8 -y
pip install umap-learn

