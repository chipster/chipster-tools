#!/bin/bash
#This tool does a temporary installation of a bioconda package
#and launches a commend from there
#
#Syntax:
#  conda_otf conda_package_name/command -with -options
#
#For example:
#  conda_otf art/art_illumina -help
#

conda_path=("/opt/chipster/tools/miniconda3")
export PATH=${conda_path}/bin:$PATH

job_dir=$(pwd)
conda_tmp_path=$job_dir/conda_tmp_$$
mkdir $conda_tmp_path 

c_env=$(echo $1 | awk -F "/" '{print $1}')
c_tool=$(echo $1 | awk -F "/" '{print $2}')
shift

conda create -y -p $conda_tmp_path $c_env  
source activate $conda_tmp_path

$c_tool $@

source deactivate
#rm -rf $conda_tmp_path
