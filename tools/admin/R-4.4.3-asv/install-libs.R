# started from comp-r-4-4-3 image 

# apt update -y && apt install -y libjpeg-turbo8-dev

install.packages("BiocManager")
# ~ 15 min
BiocManager::install("dada2", version = "3.20")

