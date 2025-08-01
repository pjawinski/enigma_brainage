#!/bin/bash

# =================================
# === prepare phenome-wide scan ===
# =================================

# clone phesant from GitHub
cd /fast/software
git clone https://github.com/MRCIEU/PHESANT.git
cd PHESANT

# run testWAS (https://github.com/MRCIEU/PHESANT/tree/master/testWAS)
cd /fast/software/PHESANT/WAS
testDir="../testWAS/"

Rscript phenomeScan.r \
--phenofile="${testDir}data/phenotypes.csv" \
--traitofinterestfile="${testDir}data/exposure.csv" \
--variablelistfile="${testDir}variable-lists/outcome-info.tsv" \
--datacodingfile="${testDir}variable-lists/data-coding-ordinal-info.txt" \
--traitofinterest="exposure" \
--resDir="${testDir}results/" \
--userId="userId"

# shortcut
Rscript phenomeScan.r --test

# run part 1 of 3
Rscript phenomeScan.r \
--phenofile="${testDir}data/phenotypes.csv" \
--traitofinterestfile="${testDir}data/exposure.csv" \
--variablelistfile="${testDir}variable-lists/outcome-info.tsv" \
--datacodingfile="${testDir}variable-lists/data-coding-ordinal-info.txt" \
--traitofinterest="exposure" \
--resDir="${testDir}results/" \
--userId="userId" \
--partIdx=1 \
--numParts=3

# run results processing
cd ../resultsProcessing/

Rscript mainCombineResults.r \
--resDir="../testWAS/results/" \
--variablelistfile="../testWAS/variable-lists/outcome-info.tsv"

