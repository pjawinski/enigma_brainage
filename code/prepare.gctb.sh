#!/bin/bash

# =================================
# === download and install gctb ===
# =================================

# set directory and download download gctb
mkdir -p /fast/software/gctb
cd /fast/software/gctb
wget https://cnsgenomics.com/software/gctb/download/gctb_2.5.2_Linux.zip
unzip gctb_2.5.2_Linux.zip
ln -s /fast/software/gctb/gctb_2.5.2_Linux/gctb /fast/software/bin/gctb

# download resources
mkdir -p /fast/software/gctb/resources
cd /fast/software/gctb/resources
wget https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/LD/Imputed/ukbEUR_Imputed.zip
wget https://sbayes.pctgplots.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/annot_baseline2.2.zip
unzip ukbEUR_Imputed.zip
unzip annot_baseline2.2.zip 
