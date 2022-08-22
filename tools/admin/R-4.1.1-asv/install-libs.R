# wget https://ftp.acc.umu.se/mirror/CRAN/src/base/R-4/R-4.1.1.tar.gz
# tar -xzf R-4.1.1.tar.gz
# cd R-4.1.1
# ./configure --with-x=no --with-pcre1 --prefix=/opt/chipster/tools/R-4.1.1-asv
# make
# make install

install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")


# mkdir dada2-silva-reference
# cd dada2-silva-reference

# wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=

# wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1

# wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1

# wget https://zenodo.org/record/4587955/files/SILVA_LICENSE.txt?download=1
 
# mv 'silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1' silva_nr99_v138.1_wSpecies_train_set.fa.gz

# gunzip silva_nr99_v138.1_wSpecies_train_set.fa.gz

# mv 'SILVA_LICENSE.txt?download=1'  SILVA_LICENSE.txt  
   
# mv 'silva_nr99_v138.1_train_set.fa.gz?download=1' silva_nr99_v138.1_train_set.fa.gz  
 
# gunzip silva_nr99_v138.1_train_set.fa.gz

# mv 'silva_species_assignment_v138.1.fa.gz?download=1' silva_species_assignment_v138.1.fa.gz
 
# gunzip silva_species_assignment_v138.1.fa.gz


