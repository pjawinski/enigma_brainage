#!/bin/bash

# =============
# === setup ===
# =============

# create some symbolic links
cd /slow/projects/enigma_brainage/
ln -s /fast/UK_Biobank/02_data_standard data/basket
ln -s /fast/UK_Biobank/04_data_genetics_linux/ data/genetics

# install conda environments
for env in default phesant; do
	if [ ! -d "envs/${env}" ]; then
		mamba env create --file envs/${env}.yml -p envs/${env}
	fi
done
conda activate envs/default

#!/usr/bin/env Rscript
# ==========================================================
# === Extract FreeSurfer variables from UKB tabular data === 
# ==========================================================
message('\n--- Extracting FreeSurfer variables from UKB tabular data ---')

# load required packages
library(dplyr)

# load 2024 basket data
message(' - loading dnanexus data (release 31/03/2025).')
bd = data.frame(data.table::fread(cmd = 'gzip -dc data/basket/dnanexus/fields.imaging.output.txt.csv.gz', header = T, sep = ','))
names(bd) = gsub("^X", "f.", names(bd))
names(bd)[1] = "f.eid"

# set variables of interest (voi)
# participant ID: f.eid
# sex (baseline): f.31.0.0 (coding: https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=9)
# age at initial imaging visit (2014+): f.21003.2.0
# scanner site at initial imaging visit (2014+): f.54.2.0
# ancestry: f.22006.0.0 
# ethnicity: f.2100.0.0
# handedness: f.1707.0.0
# iq: NA
# dx (length of illness, disease category, sub-type): NA
# ses (the following could perhaps be used as proxy): Townsend deprivation index at recruitment (field 22189) | category 76 (education/employment/housing scores)
# genetic principal components: f.22009.0.1 - f.22009.0.40
voi_main = c('f.eid','f.31.0.0','f.21003.2.0','f.54.2.0','f.22006.0.0','f.21000.0.0','f.1707.0.0',paste0('f.22009.0.',1:10))
message(sprintf(' - %d/%d requested variables available', sum(voi_main %in% names(bd)), length(voi_main)))
voi_main_renamed = c('SUBJID','SEX','AGE','SITE','ANCESTRY','ETHNICITY','HAND',paste0('PC',1:10))
rename_map = setNames(voi_main_renamed, voi_main)
names(bd) = ifelse(names(bd) %in% names(rename_map), rename_map[names(bd)], names(bd))

# recode some variables
bd$SEX = factor(bd$SEX, levels = c(0,1), labels = c(1,0)) %>% as.character %>% as.numeric # bd$SEX = factor(bd$SEX, levels = c('Male','Female'), labels = c('0','1')) %>% as.character %>% as.numeric
bd$HAND = factor(bd$HAND, levels = c(1,2,3,-3), labels = c(0,1,2,NA)) %>% as.character %>% as.numeric # bd$HAND = factor(bd$HAND, levels = c('Right-handed','Left-handed','Use both right and left hands equally','Prefer not to answer'), labels = c('0','1','2','NA')) %>% as.character %>% as.numeric
bd$ANCESTRY = factor(bd$ANCESTRY, levels = c(1), labels = c('European')) %>% as.character # bd$ANCESTRY = factor(bd$ANCESTRY, levels = c('Caucasian'), labels = c('European')) %>% as.character
bd$SITE = bd$SITE %>% as.numeric %>% scales::rescale(to = c(1, length(table(bd$SITE))))
bd$ETHNICITY = factor(bd$ETHNICITY, levels = c(1,1001,2001,3001,4001,2,1002,2002,3002,4002,3,1003,2003,3003,4003,4,2004,3004,5,6,-1,-3), labels = c('White','British','White and Black Caribbean','Indian','Caribbean','Mixed','Irish','White and Black African','Pakistani','African','Asian or Asian British','Any other white background','White and Asian','Bangladeshi','Any other Black background','Black or Black British','Any other mixed background','Any other Asian background','Chinese','Other ethnic group','Do not know','Prefer not to answer')) %>% as.character

# surface area (SA): https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=192
# left hemisphere: f.26721.2.0 to f.26754.2.0
# right hemisphere: f.26822.2.0 to f.26855.2.0
# note that temporal poles are missing
voi_sa = c(paste0('f.',26721:26754,'.2.0'),paste0('f.',26822:26855,'.2.0'))
message(sprintf(' - %d/%d surface area variables available', sum(voi_sa %in% names(bd)), length(voi_sa)))
voi_sa_renamed = c('LSurfArea','L_bankssts_surfavg','L_caudalanteriorcingulate_surfavg','L_caudalmiddlefrontal_surfavg','L_cuneus_surfavg','L_entorhinal_surfavg','L_fusiform_surfavg','L_inferiorparietal_surfavg','L_inferiortemporal_surfavg','L_isthmuscingulate_surfavg','L_lateraloccipital_surfavg','L_lateralorbitofrontal_surfavg','L_lingual_surfavg','L_medialorbitofrontal_surfavg','L_middletemporal_surfavg','L_parahippocampal_surfavg','L_paracentral_surfavg','L_parsopercularis_surfavg','L_parsorbitalis_surfavg','L_parstriangularis_surfavg','L_pericalcarine_surfavg','L_postcentral_surfavg','L_posteriorcingulate_surfavg','L_precentral_surfavg','L_precuneus_surfavg','L_rostralanteriorcingulate_surfavg','L_rostralmiddlefrontal_surfavg','L_superiorfrontal_surfavg','L_superiorparietal_surfavg','L_superiortemporal_surfavg','L_supramarginal_surfavg','L_frontalpole_surfavg','L_transversetemporal_surfavg','L_insula_surfavg','RSurfArea','R_bankssts_surfavg','R_caudalanteriorcingulate_surfavg','R_caudalmiddlefrontal_surfavg','R_cuneus_surfavg','R_entorhinal_surfavg','R_fusiform_surfavg','R_inferiorparietal_surfavg','R_inferiortemporal_surfavg','R_isthmuscingulate_surfavg','R_lateraloccipital_surfavg','R_lateralorbitofrontal_surfavg','R_lingual_surfavg','R_medialorbitofrontal_surfavg','R_middletemporal_surfavg','R_parahippocampal_surfavg','R_paracentral_surfavg','R_parsopercularis_surfavg','R_parsorbitalis_surfavg','R_parstriangularis_surfavg','R_pericalcarine_surfavg','R_postcentral_surfavg','R_posteriorcingulate_surfavg','R_precentral_surfavg','R_precuneus_surfavg','R_rostralanteriorcingulate_surfavg','R_rostralmiddlefrontal_surfavg','R_superiorfrontal_surfavg','R_superiorparietal_surfavg','R_superiortemporal_surfavg','R_supramarginal_surfavg','R_frontalpole_surfavg','R_transversetemporal_surfavg','R_insula_surfavg')
rename_map = setNames(voi_sa_renamed, voi_sa)
names(bd) = ifelse(names(bd) %in% names(rename_map), rename_map[names(bd)], names(bd))

# cortical thickness: https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=192
# left hemisphere: f.26755.2.0 to f.26788.2.0
# right hemisphere: f.26856.2.0 to f.26889.2.0
# note that temporal poles are missing
voi_ct = c(paste0('f.',26755:26788,'.2.0'),paste0('f.',26856:26889,'.2.0'))
message(sprintf(' - %d/%d cortical thickness variables available', sum(voi_ct %in% names(bd)), length(voi_ct)))
voi_ct_renamed = c('LThickness','L_bankssts_thickavg','L_caudalanteriorcingulate_thickavg','L_caudalmiddlefrontal_thickavg','L_cuneus_thickavg','L_entorhinal_thickavg','L_fusiform_thickavg','L_inferiorparietal_thickavg','L_inferiortemporal_thickavg','L_isthmuscingulate_thickavg','L_lateraloccipital_thickavg','L_lateralorbitofrontal_thickavg','L_lingual_thickavg','L_medialorbitofrontal_thickavg','L_middletemporal_thickavg','L_parahippocampal_thickavg','L_paracentral_thickavg','L_parsopercularis_thickavg','L_parsorbitalis_thickavg','L_parstriangularis_thickavg','L_pericalcarine_thickavg','L_postcentral_thickavg','L_posteriorcingulate_thickavg','L_precentral_thickavg','L_precuneus_thickavg','L_rostralanteriorcingulate_thickavg','L_rostralmiddlefrontal_thickavg','L_superiorfrontal_thickavg','L_superiorparietal_thickavg','L_superiortemporal_thickavg','L_supramarginal_thickavg','L_frontalpole_thickavg','L_transversetemporal_thickavg','L_insula_thickavg','RThickness','R_bankssts_thickavg','R_caudalanteriorcingulate_thickavg','R_caudalmiddlefrontal_thickavg','R_cuneus_thickavg','R_entorhinal_thickavg','R_fusiform_thickavg','R_inferiorparietal_thickavg','R_inferiortemporal_thickavg','R_isthmuscingulate_thickavg','R_lateraloccipital_thickavg','R_lateralorbitofrontal_thickavg','R_lingual_thickavg','R_medialorbitofrontal_thickavg','R_middletemporal_thickavg','R_parahippocampal_thickavg','R_paracentral_thickavg','R_parsopercularis_thickavg','R_parsorbitalis_thickavg','R_parstriangularis_thickavg','R_pericalcarine_thickavg','R_postcentral_thickavg','R_posteriorcingulate_thickavg','R_precentral_thickavg','R_precuneus_thickavg','R_rostralanteriorcingulate_thickavg','R_rostralmiddlefrontal_thickavg','R_superiorfrontal_thickavg','R_superiorparietal_thickavg','R_superiortemporal_thickavg','R_supramarginal_thickavg','R_frontalpole_thickavg','R_transversetemporal_thickavg','R_insula_thickavg')
rename_map = setNames(voi_ct_renamed, voi_ct)
names(bd) = ifelse(names(bd) %in% names(rename_map), rename_map[names(bd)], names(bd))

# subcortical volumes and ICV (ASEG): https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=192
# left hemisphere: see names below
# right hemisphere: see names below
# extract columns with desirable order
voi_sv = c('f.26554.2.0','f.26585.2.0','f.26558.2.0','f.26589.2.0','f.26559.2.0','f.26590.2.0','f.26560.2.0','f.26591.2.0','f.26561.2.0','f.26592.2.0','f.26562.2.0','f.26593.2.0','f.26563.2.0','f.26594.2.0','f.26564.2.0','f.26595.2.0','f.26521.2.0')
message(sprintf(' - %d/%d subcortical volumes + ICV variables available', sum(voi_sv %in% names(bd)), length(voi_sv)))
voi_sv_renamed = c('LLatVent','RLatVent','Lthal','Rthal','Lcaud','Rcaud','Lput','Rput','Lpal','Rpal','Lhippo','Rhippo','Lamyg','Ramyg','Laccumb','Raccumb','ICV')
rename_map = setNames(voi_sv_renamed, voi_sv)
names(bd) = ifelse(names(bd) %in% names(rename_map), rename_map[names(bd)], names(bd))

# set header
header_sa = c('SUBJID','L_bankssts_surfavg','L_caudalanteriorcingulate_surfavg','L_caudalmiddlefrontal_surfavg','L_cuneus_surfavg','L_entorhinal_surfavg','L_fusiform_surfavg','L_inferiorparietal_surfavg','L_inferiortemporal_surfavg','L_isthmuscingulate_surfavg','L_lateraloccipital_surfavg','L_lateralorbitofrontal_surfavg','L_lingual_surfavg','L_medialorbitofrontal_surfavg','L_middletemporal_surfavg','L_parahippocampal_surfavg','L_paracentral_surfavg','L_parsopercularis_surfavg','L_parsorbitalis_surfavg','L_parstriangularis_surfavg','L_pericalcarine_surfavg','L_postcentral_surfavg','L_posteriorcingulate_surfavg','L_precentral_surfavg','L_precuneus_surfavg','L_rostralanteriorcingulate_surfavg','L_rostralmiddlefrontal_surfavg','L_superiorfrontal_surfavg','L_superiorparietal_surfavg','L_superiortemporal_surfavg','L_supramarginal_surfavg','L_frontalpole_surfavg','L_temporalpole_surfavg','L_transversetemporal_surfavg','L_insula_surfavg','R_bankssts_surfavg','R_caudalanteriorcingulate_surfavg','R_caudalmiddlefrontal_surfavg','R_cuneus_surfavg','R_entorhinal_surfavg','R_fusiform_surfavg','R_inferiorparietal_surfavg','R_inferiortemporal_surfavg','R_isthmuscingulate_surfavg','R_lateraloccipital_surfavg','R_lateralorbitofrontal_surfavg','R_lingual_surfavg','R_medialorbitofrontal_surfavg','R_middletemporal_surfavg','R_parahippocampal_surfavg','R_paracentral_surfavg','R_parsopercularis_surfavg','R_parsorbitalis_surfavg','R_parstriangularis_surfavg','R_pericalcarine_surfavg','R_postcentral_surfavg','R_posteriorcingulate_surfavg','R_precentral_surfavg','R_precuneus_surfavg','R_rostralanteriorcingulate_surfavg','R_rostralmiddlefrontal_surfavg','R_superiorfrontal_surfavg','R_superiorparietal_surfavg','R_superiortemporal_surfavg','R_supramarginal_surfavg','R_frontalpole_surfavg','R_temporalpole_surfavg','R_transversetemporal_surfavg','R_insula_surfavg','LThickness','RThickness','LSurfArea','RSurfArea','ICV')
header_ct = c('SUBJID','L_bankssts_thickavg','L_caudalanteriorcingulate_thickavg','L_caudalmiddlefrontal_thickavg','L_cuneus_thickavg','L_entorhinal_thickavg','L_fusiform_thickavg','L_inferiorparietal_thickavg','L_inferiortemporal_thickavg','L_isthmuscingulate_thickavg','L_lateraloccipital_thickavg','L_lateralorbitofrontal_thickavg','L_lingual_thickavg','L_medialorbitofrontal_thickavg','L_middletemporal_thickavg','L_parahippocampal_thickavg','L_paracentral_thickavg','L_parsopercularis_thickavg','L_parsorbitalis_thickavg','L_parstriangularis_thickavg','L_pericalcarine_thickavg','L_postcentral_thickavg','L_posteriorcingulate_thickavg','L_precentral_thickavg','L_precuneus_thickavg','L_rostralanteriorcingulate_thickavg','L_rostralmiddlefrontal_thickavg','L_superiorfrontal_thickavg','L_superiorparietal_thickavg','L_superiortemporal_thickavg','L_supramarginal_thickavg','L_frontalpole_thickavg','L_temporalpole_thickavg','L_transversetemporal_thickavg','L_insula_thickavg','R_bankssts_thickavg','R_caudalanteriorcingulate_thickavg','R_caudalmiddlefrontal_thickavg','R_cuneus_thickavg','R_entorhinal_thickavg','R_fusiform_thickavg','R_inferiorparietal_thickavg','R_inferiortemporal_thickavg','R_isthmuscingulate_thickavg','R_lateraloccipital_thickavg','R_lateralorbitofrontal_thickavg','R_lingual_thickavg','R_medialorbitofrontal_thickavg','R_middletemporal_thickavg','R_parahippocampal_thickavg','R_paracentral_thickavg','R_parsopercularis_thickavg','R_parsorbitalis_thickavg','R_parstriangularis_thickavg','R_pericalcarine_thickavg','R_postcentral_thickavg','R_posteriorcingulate_thickavg','R_precentral_thickavg','R_precuneus_thickavg','R_rostralanteriorcingulate_thickavg','R_rostralmiddlefrontal_thickavg','R_superiorfrontal_thickavg','R_superiorparietal_thickavg','R_superiortemporal_thickavg','R_supramarginal_thickavg','R_frontalpole_thickavg','R_temporalpole_thickavg','R_transversetemporal_thickavg','R_insula_thickavg','LThickness','RThickness','LSurfArea','RSurfArea','ICV')
header_sv = c('SUBJID','LLatVent','RLatVent','Lthal','Rthal','Lcaud','Rcaud','Lput','Rput','Lpal','Rpal','Lhippo','Rhippo','Lamyg','Ramyg','Laccumb','Raccumb','ICV')
header_covs = c('SUBJID','DX','SEX','AGE',paste0('PC',1:10),'SITE','ETHNICITY','ANCESTRY','SES','IQ','HAND','LENGTH_OF_ILLNESS','DISEASE_CATEGORY','SUBTYPE')

# add missing columns with NA
requiredCols = c(header_ct,header_sa,header_sv,header_covs)
missingCols = requiredCols[!(requiredCols %in% names(bd))]
message(sprintf('The following columns are missing: %s',paste(missingCols,collapse = ' | ')))
for (col in missingCols) { bd[[col]] = NA }

# keep cases with available age, sex, site, and PC1:10 information (from initial imaging visit)
message(' - keeping cases with available age and site information (from initial imaging visit)')
bd = bd[complete.cases(bd[,c('AGE','SEX','SITE',paste0('PC',1:10))]),]
message(sprintf(' - %d cases remaining', nrow(bd)))

# keep cases with available information from LThickness, RThickness, LSurfArea, RSurfArea, ICV
message(' - keeping cases with available information from LThickness, RThickness, LSurfArea, RSurfArea, ICV')
bd = bd[complete.cases(bd[,c('LThickness','RThickness','LSurfArea','RSurfArea','ICV')]),]
message(sprintf(' - %d cases remaining', nrow(bd)))

# create output data frame
sa = bd[,header_sa]
ct = bd[,header_ct]
sv = bd[,header_sv]
covs = bd[,header_covs]

# get summary
summary(sa) # no L_temporalpole_surfavg & R_temporalpole_surfavg data, 1 NA in R_parahippocampal_surfavg
summary(ct) # no L_temporalpole_surfavg & R_temporalpole_surfavg data, 1 NA in R_parahippocampal_surfavg
summary(sv) # 1 NA in Lthal in Rthal
summary(covs)  # no data in DX, SES, IQ, LENGTH_OF_ILLNESS, DISEASE_CATEGORY, SUBTYPE | 15 NAs in ETHNICITY, 

# export CSV files
write.csv(covs,'code/ENIGMA-brainage-local/Covariates.csv', row.names=FALSE)
write.csv(sa,'code/ENIGMA-brainage-local/CorticalMeasuresENIGMA_SurfAvg.csv', row.names=FALSE)
write.csv(ct,'code/ENIGMA-brainage-local/CorticalMeasuresENIGMA_ThickAvg.csv', row.names=FALSE)
write.csv(sv,'code/ENIGMA-brainage-local/SubcorticalMeasuresENIGMA_VolAvg.csv', row.names=FALSE)

#!/bin/bash
# =====================================================
# === run photon.ai models in singularity container ===
# =====================================================

# needs to be added


#!/usr/bin/env Rscript
# ================================================================================
# === This script will save the phenotype and covariates in "pheno_covars.txt" ===
# ================================================================================

# unzip files produced by photon.ai scripts
cd /slow/projects/enigma_brainage/results
unzip UKBB_backup_20250408.zip
unzip UKBB_pheno_20250408.zip
cd ..

# load required packages
library(dplyr)

# read in the dataset with the brain age phenotype and the covariates 
load("results/output-backup/data_pheno_covariates.Rdata") # this file was created during the phenotypic analysis

# Add the first two columns as FID (Family ID) and IID (Individual ID).
# remove columns with NA only
# select relevant variables (for UKBB: 10 PCs as computed by the UKBB)
data$FID = data$IID = data$SUBJID
data = data[,colSums(is.na(data))<nrow(data)]
data = data %>% select(FID, IID, SEX, AGE, AGE2, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, SITE, ETHNICITY, ANCESTRY, SITE.f2, SITE.f3, SITE.f4, predAge, devAge, AE, ICV) 

# UKBB covariates should (ideally) also include: array batch (Field Code: 22000)
# array batch converted to binary variable (genotyping array)
pj = read.delim('/slow/projects/ukb_brainage/results/mri/r2024.vars.txt')
pj = pj[pj$discovery == 1 | !is.na(pj$pan),]
pj2k = read.delim('/slow/projects/ukb_brainage/data/traits/replicate.2k.txt')
df = inner_join(data, pj[,c('IID','discovery','pan','array')], by = 'IID')
df$heldout = as.numeric(df$IID %in% pj2k$IID)

# remove spaces from ETHNICITY
df$ETHNICITY = gsub(" ", "_", df$ETHNICITY)

# get discovery, replication, and left-out sample
cols = c('FID','IID','SEX','AGE','AGE2','array',paste0('PC',1:10),'SITE','ETHNICITY','ANCESTRY','SITE.f2','SITE.f3','SITE.f4','predAge','devAge','AE','ICV')
discov1 = df[df$discovery==1,cols]
discov1_male = df[df$discovery==1 & df$SEX==0,cols]
discov1_female = df[df$discovery==1 & df$SEX==1 ,cols]
discov2 = df[df$discovery!=1 & df$pan=='EUR' & df$heldout==0,cols]
discov2_male = df[df$discovery!=1 & df$pan=='EUR' & df$heldout==0 & df$SEX==0,cols]
discov2_female = df[df$discovery!=1 & df$pan=='EUR' & df$heldout==0 & df$SEX==1,cols]
discov1_2 = df[df$discovery==1 | (df$discovery!=1 & df$pan=='EUR' & df$heldout==0),cols]
discov1_2_male = df[(df$discovery==1 | (df$discovery!=1 & df$pan=='EUR' & df$heldout==0)) & df$SEX==0,cols]
discov1_2_female = df[(df$discovery==1 | (df$discovery!=1 & df$pan=='EUR' & df$heldout==0)) & df$SEX==1,cols]
heldout = df[df$pan!='EUR' | df$heldout==1,]

# save to text file
write.table(discov1, file = 'results/output-backup/pheno_covars_discov1.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov2, file = 'results/output-backup/pheno_covars_discov2.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov1_2, file = 'results/output-backup/pheno_covars_discov1-2.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov1_male, file = 'results/output-backup/pheno_covars_discov1_male.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov1_female, file = 'results/output-backup/pheno_covars_discov1_female.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov2_male, file = 'results/output-backup/pheno_covars_discov2_male.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov2_female, file = 'results/output-backup/pheno_covars_discov2_female.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov1_2_male, file = 'results/output-backup/pheno_covars_discov1-2_male.txt', sep = " ", col.names = T, row.names = F, quote = F)
write.table(discov1_2_female, file = 'results/output-backup/pheno_covars_discov1-2_female.txt', sep = " ", col.names = T, row.names = F, quote = F)

#!/bin/bash
# ===================================================================================
# === Prepare genetics for ukbb.discovery1, ukbb.discovery, and ukbb.discovery1-2 === 
# ===================================================================================

# create bgen-1.2 with 8 bits (instead of 16)
for sample in "" "_EUR" "_EURjoined"; do # ukbb.discovery1: "" | ukbb.discovery2: "_EUR" | ukbb.discovery1-2: "_EURjoined"
	N=2; (for i in {1..22}; do 
	   ((j=j%N)); ((j++==0)) && wait
	      (mkdir -p data/genetics/chr${i}/imp_mri_qc${sample}/bgen_8bits
	      plink2 \
	      --pfile data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc \
	      --export bgen-1.2 bits=8 \
	      --out data/genetics/chr${i}/imp_mri_qc${sample}/bgen_8bits/chr${i}_mri_qc
	   ) &
	done
	wait)
done

	# create diploid X chromosome data for males by setting the sex of all individuals to female
	# https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html
	i=X
	for sample in "" "_EUR" "_EURjoined"; do # ukbb.discovery1: "" | ukbb.discovery2: "_EUR" | ukbb.discovery1-2: "_EURjoined"
		mkdir -p data/genetics/chr${i}/imp_mri_qc${sample}/bgen_8bits
		# cp data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.psam data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.original.psam
		# awk 'NR==1 { print } NR > 1 { print $1, $2, 2 }' OFS='\t' data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.psam > data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.allfemale.psam
		\cp data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.allfemale.psam data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.psam
		plink2 \
		  --pfile data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc \
		  --export bgen-1.2 bits=8 \
		  --out data/genetics/chr${i}/imp_mri_qc${sample}/bgen_8bits/chr${i}_mri_qc
		\cp data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.original.psam data/genetics/chr${i}/imp_mri_qc${sample}/chr${i}_mri_qc.psam
	done

# prune QC-ed genetic data for BOLT-LMM
# qc done so far: selection of biallelic variants with maf 0.001; no hwe or geno-thresholds applied post imputation
# suggested qc for .bed files (hard-calls): geno = 0.02, hwe = 1e-8, maf = 0.01
for sample in "" "_EUR" "_EURjoined"; do # ukbb.discovery1: "" | ukbb.discovery2: "_EUR" | ukbb.discovery1-2: "_EURjoined"
	for i in {1..22} X; do
		mkdir -p "data/genetics/chr${i}/imp_mri_qc${sample}/bed/pruned"
		plink2 --bfile "data/genetics/chr${i}/imp_mri_qc${sample}/bed/chr${i}_mri_qc" \
			--geno 0.02 \
			--hwe 1e-8 \
			--maf 0.01 \
			--indep-pairwise 200 50 0.2 \
			--threads 100 \
			--out "data/genetics/chr${i}/imp_mri_qc${sample}/bed/pruned/chr${i}_mri_qc"
		plink2 --bfile "data/genetics/chr${i}/imp_mri_qc${sample}/bed/chr${i}_mri_qc" \
			--extract "data/genetics/chr${i}/imp_mri_qc${sample}/bed/pruned/chr${i}_mri_qc.prune.in" \
			--make-bed \
			--out "data/genetics/chr${i}/imp_mri_qc${sample}/bed/pruned/chr${i}_mri_qc_pruned"
		chmod 770 "data/genetics/chr${i}/imp_mri_qc${sample}/bed/pruned/chr${i}_mri_qc"*
	done
done

# ~500,000 SNPs should remain in the pruned file
cat $(for i in {1..22} X; do echo "data/genetics/chr${i}/imp_mri_qc/bed/pruned/chr${i}_mri_qc_pruned.bim"; done) | wc -l # 539,428 remaining (including column header) in ukbb.discovery1
cat $(for i in {1..22} X; do echo "data/genetics/chr${i}/imp_mri_qc_EUR/bed/pruned/chr${i}_mri_qc_pruned.bim"; done) | wc -l # 536,967 remaining (including column header) in ukbb.discovery2
cat $(for i in {1..22} X; do echo "data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/pruned/chr${i}_mri_qc_pruned.bim"; done) | wc -l # 537,793 remaining (including column header) in ukbb.discovery2

# =======================================
# === Run BOLT-LMM in ukbb.discovery1 === 
# =======================================

# set working directory and load conda environment
cd /slow/projects/enigma_brainage
conda activate envs/default

# create temporary list of bgen and sample files
targetDir="results/bolt"
mkdir -p "${targetDir}"/tmp
> "${targetDir}"/tmp/bolt.filelist
for i in {1..22} X; do
	echo data/genetics/chr${i}/imp_mri_qc/bgen_8bits/chr${i}_mri_qc.bgen data/genetics/chr${i}/imp_mri_qc/bgen_8bits/chr${i}_mri_qc.sample >> "${targetDir}"/tmp/bolt.filelist
done

# create temporary symbolic links to pruned bfiles
for i in {1..22}; do
	for type in bed bim fam; do
		ln -sf $(pwd)/data/genetics/chr${i}/imp_mri_qc/bed/pruned/chr${i}_mri_qc_pruned.${type} ${targetDir}/tmp/chr${i}_mri_qc_pruned.${type}
	done
done
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc/bed/pruned/chrX_mri_qc_pruned.bed ${targetDir}/tmp/chr23_mri_qc_pruned.bed
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc/bed/pruned/chrX_mri_qc_pruned.bim ${targetDir}/tmp/chr23_mri_qc_pruned.bim
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc/bed/pruned/chrX_mri_qc_pruned.fam ${targetDir}/tmp/chr23_mri_qc_pruned.fam

# model 1 (devAge) + model 2 (predAge)
BOLT_PATH="/fast/software/BOLT-LMM/BOLT-LMM_v2.4.1/" # path to BOLT-LMM
PHENO_PATH="$(pwd)/results/output-backup" # Path to 'pheno_covars.txt'
for PHENO in predAge devAge; do
	bolt \
		--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
		--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
		--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
		--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
		--LDscoresMatchBp \
		--lmmForceNonInf \
		--noBgenIDcheck \
		--covarUseMissingIndic \
		--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
		--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
		--numThreads=56 \
		--lmm \
		--phenoFile=${PHENO_PATH}/pheno_covars_discov1.txt \
		--phenoCol=${PHENO} \
		--covarFile=${PHENO_PATH}/pheno_covars_discov1.txt \
		--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
		--qCovarCol={PC{1:10},AGE,AGE2,SEX,ICV} \
		--covarMaxLevels 1000 \
		--verboseStats \
		--h2gGuess=0.10 \
		--statsFile="${targetDir}"/ukbb.discov1.chr1-23_non_imputed_${PHENO}.stats \
		--statsFileBgenSnps="${targetDir}"/ukbb.discov1.chr1-23_imputed_${PHENO}.stats \
		--numLeaveOutChunks 2 \
		--bgenMinMAF=0.01 \
		--bgenMinINFO=0.3 \
		2>&1 | tee ${targetDir}/ukbb.discov1.chr1-23_imputed_${PHENO}_$(date +'%Y%m%d').log
	pigz -f "${targetDir}"/ukbb.discov1.chr1-23_non_imputed_${PHENO}.stats
	pigz -f "${targetDir}"/ukbb.discov1.chr1-23_imputed_${PHENO}.stats
done
chmod 770 ${targetDir}/*

# model 3 (noICV)
bolt \
	--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
	--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
	--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
	--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
	--LDscoresMatchBp \
	--lmmForceNonInf \
	--noBgenIDcheck \
	--covarUseMissingIndic \
	--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
	--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
	--numThreads=56 \
	--lmm \
	--phenoFile=${PHENO_PATH}/pheno_covars_discov1.txt \
	--phenoCol=devAge \
	--covarFile=${PHENO_PATH}/pheno_covars_discov1.txt \
	--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
	--qCovarCol={PC{1:10},AGE,AGE2,SEX} \
	--covarMaxLevels 1000 \
	--verboseStats \
	--h2gGuess=0.10 \
	--statsFile="${targetDir}"/ukbb.discov1.chr1-23_non_imputed_devAge_noICV.stats \
	--statsFileBgenSnps="${targetDir}"/ukbb.discov1.chr1-23_imputed_devAge_noICV.stats \
	--numLeaveOutChunks 2 \
	--bgenMinMAF=0.01 \
	--bgenMinINFO=0.3 \
	2>&1 | tee ${targetDir}/ukbb.discov1.chr1-23_imputed_devAge_noICV_$(date +'%Y%m%d').log
chmod 770 ${targetDir}/*
pigz -f "${targetDir}"/ukbb.discov1.chr1-23_non_imputed_devAge_noICV.stats
pigz -f "${targetDir}"/ukbb.discov1.chr1-23_imputed_devAge_noICV.stats

# model 4 + 5 (sex-stratified)
for SEX in male female; do
	bolt \
		--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
		--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
		--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
		--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
		--LDscoresMatchBp \
		--lmmForceNonInf \
		--noBgenIDcheck \
		--covarUseMissingIndic \
		--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
		--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
		--numThreads=56 \
		--lmm \
		--phenoFile=${PHENO_PATH}/pheno_covars_discov1_${SEX}.txt \
		--phenoCol=devAge \
		--covarFile=${PHENO_PATH}/pheno_covars_discov1_${SEX}.txt \
		--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
		--qCovarCol={PC{1:10},AGE,AGE2,ICV} \
		--covarMaxLevels 1000 \
		--verboseStats \
		--h2gGuess=0.10 \
		--statsFile="${targetDir}"/ukbb.discov1.chr1-23_non_imputed_devAge_${SEX}.stats \
		--statsFileBgenSnps="${targetDir}"/ukbb.discov1.chr1-23_imputed_devAge_${SEX}.stats \
		--numLeaveOutChunks 2 \
		--bgenMinMAF=0.01 \
		--bgenMinINFO=0.3 \
		2>&1 | tee ${targetDir}/ukbb.discov1.chr1-23_imputed_devAge_${SEX}_$(date +'%Y%m%d').log
done
chmod 770 ${targetDir}/*
pigz -f "${targetDir}"/ukbb.discov1.chr1-23_non_imputed_devAge_female.stats
pigz -f "${targetDir}"/ukbb.discov1.chr1-23_imputed_devAge_female.stats
pigz -f "${targetDir}"/ukbb.discov1.chr1-23_non_imputed_devAge_male.stats
pigz -f "${targetDir}"/ukbb.discov1.chr1-23_imputed_devAge_male.stats
rm -rf ${targetDir}/tmp

# =======================================
# === Run BOLT-LMM in ukbb.discovery2 === 
# =======================================

# set working directory and load conda environment
cd /slow/projects/enigma_brainage
conda activate envs/default

# create temporary list of bgen and sample files
targetDir="results/bolt"
mkdir -p "${targetDir}"/tmp
> "${targetDir}"/tmp/bolt.filelist
for i in {1..22} X; do
	echo data/genetics/chr${i}/imp_mri_qc_EUR/bgen_8bits/chr${i}_mri_qc.bgen data/genetics/chr${i}/imp_mri_qc_EUR/bgen_8bits/chr${i}_mri_qc.sample >> "${targetDir}"/tmp/bolt.filelist
done

# create temporary symbolic links to pruned bfiles
for i in {1..22}; do
	for type in bed bim fam; do
		ln -sf $(pwd)/data/genetics/chr${i}/imp_mri_qc_EUR/bed/pruned/chr${i}_mri_qc_pruned.${type} ${targetDir}/tmp/chr${i}_mri_qc_pruned.${type}
	done
done
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc_EUR/bed/pruned/chrX_mri_qc_pruned.bed ${targetDir}/tmp/chr23_mri_qc_pruned.bed
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc_EUR/bed/pruned/chrX_mri_qc_pruned.bim ${targetDir}/tmp/chr23_mri_qc_pruned.bim
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc_EUR/bed/pruned/chrX_mri_qc_pruned.fam ${targetDir}/tmp/chr23_mri_qc_pruned.fam

# model 1 (devAge) + model 2 (predAge)
BOLT_PATH="/fast/software/BOLT-LMM/BOLT-LMM_v2.4.1/" # path to BOLT-LMM
PHENO_PATH="$(pwd)/results/output-backup" # Path to 'pheno_covars.txt'
for PHENO in devAge predAge; do
	bolt \
		--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
		--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
		--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
		--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
		--LDscoresMatchBp \
		--lmmForceNonInf \
		--noBgenIDcheck \
		--covarUseMissingIndic \
		--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
		--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
		--numThreads=56 \
		--lmm \
		--phenoFile=${PHENO_PATH}/pheno_covars_discov2.txt \
		--phenoCol=${PHENO} \
		--covarFile=${PHENO_PATH}/pheno_covars_discov2.txt \
		--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
		--qCovarCol={PC{1:10},AGE,AGE2,SEX,ICV} \
		--covarMaxLevels 1000 \
		--verboseStats \
		--h2gGuess=0.10 \
		--statsFile="${targetDir}"/ukbb.discov2.chr1-23_non_imputed_${PHENO}.stats \
		--statsFileBgenSnps="${targetDir}"/ukbb.discov2.chr1-23_imputed_${PHENO}.stats \
		--numLeaveOutChunks 2 \
		--bgenMinMAF=0.01 \
		--bgenMinINFO=0.3 \
		2>&1 | tee ${targetDir}/ukbb.discov2.chr1-23_imputed_${PHENO}_$(date +'%Y%m%d').log
	pigz -f "${targetDir}"/ukbb.discov2.chr1-23_non_imputed_${PHENO}.stats
	pigz -f "${targetDir}"/ukbb.discov2.chr1-23_imputed_${PHENO}.stats
done

# model 3 (noICV)
bolt \
	--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
	--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
	--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
	--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
	--LDscoresMatchBp \
	--lmmForceNonInf \
	--noBgenIDcheck \
	--covarUseMissingIndic \
	--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
	--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
	--numThreads=56 \
	--lmm \
	--phenoFile=${PHENO_PATH}/pheno_covars_discov2.txt \
	--phenoCol=devAge \
	--covarFile=${PHENO_PATH}/pheno_covars_discov2.txt \
	--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
	--qCovarCol={PC{1:10},AGE,AGE2,SEX} \
	--covarMaxLevels 1000 \
	--verboseStats \
	--h2gGuess=0.10 \
	--statsFile="${targetDir}"/ukbb.discov2.chr1-23_non_imputed_devAge_noICV.stats \
	--statsFileBgenSnps="${targetDir}"/ukbb.discov2.chr1-23_imputed_devAge_noICV.stats \
	--numLeaveOutChunks 2 \
	--bgenMinMAF=0.01 \
	--bgenMinINFO=0.3 \
	2>&1 | tee ${targetDir}/ukbb.discov2.chr1-23_imputed_devAge_noICV_$(date +'%Y%m%d').log
pigz -f "${targetDir}"/ukbb.discov2.chr1-23_non_imputed_devAge_noICV.stats
pigz -f "${targetDir}"/ukbb.discov2.chr1-23_imputed_devAge_noICV.stats

# model 4 + 5 (sex-stratified)
for SEX in male female; do
	bolt \
		--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
		--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
		--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
		--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
		--LDscoresMatchBp \
		--lmmForceNonInf \
		--noBgenIDcheck \
		--covarUseMissingIndic \
		--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
		--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
		--numThreads=56 \
		--lmm \
		--phenoFile=${PHENO_PATH}/pheno_covars_discov2_${SEX}.txt \
		--phenoCol=devAge \
		--covarFile=${PHENO_PATH}/pheno_covars_discov2_${SEX}.txt \
		--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
		--qCovarCol={PC{1:10},AGE,AGE2,ICV} \
		--covarMaxLevels 1000 \
		--verboseStats \
		--h2gGuess=0.10 \
		--statsFile="${targetDir}"/ukbb.discov2.chr1-23_non_imputed_devAge_${SEX}.stats \
		--statsFileBgenSnps="${targetDir}"/ukbb.discov2.chr1-23_imputed_devAge_${SEX}.stats \
		--numLeaveOutChunks 2 \
		--bgenMinMAF=0.01 \
		--bgenMinINFO=0.3 \
		2>&1 | tee ${targetDir}/ukbb.discov2.chr1-23_imputed_devAge_${SEX}_$(date +'%Y%m%d').log
	pigz -f "${targetDir}"/ukbb.discov2.chr1-23_non_imputed_devAge_${SEX}.stats
	pigz -f "${targetDir}"/ukbb.discov2.chr1-23_imputed_devAge_${SEX}.stats
done
chmod 770 ${targetDir}/*
rm -rf ${targetDir}/tmp

# =========================================
# === Run BOLT-LMM in ukbb.discovery1+2 === 
# =========================================

# set working directory and load conda environment
cd /slow/projects/enigma_brainage
conda activate envs/default

# create temporary list of bgen and sample files
targetDir="results/bolt"
mkdir -p "${targetDir}"/tmp
> "${targetDir}"/tmp/bolt.filelist
for i in {1..22} X; do
	echo data/genetics/chr${i}/imp_mri_qc_EURjoined/bgen_8bits/chr${i}_mri_qc.bgen data/genetics/chr${i}/imp_mri_qc_EURjoined/bgen_8bits/chr${i}_mri_qc.sample >> "${targetDir}"/tmp/bolt.filelist
done

# create temporary symbolic links to pruned bfiles
for i in {1..22}; do
	for type in bed bim fam; do
		ln -sf $(pwd)/data/genetics/chr${i}/imp_mri_qc_EURjoined/bed/pruned/chr${i}_mri_qc_pruned.${type} ${targetDir}/tmp/chr${i}_mri_qc_pruned.${type}
	done
done
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc_EURjoined/bed/pruned/chrX_mri_qc_pruned.bed ${targetDir}/tmp/chr23_mri_qc_pruned.bed
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc_EURjoined/bed/pruned/chrX_mri_qc_pruned.bim ${targetDir}/tmp/chr23_mri_qc_pruned.bim
ln -sf $(pwd)/data/genetics/chrX/imp_mri_qc_EURjoined/bed/pruned/chrX_mri_qc_pruned.fam ${targetDir}/tmp/chr23_mri_qc_pruned.fam

# model 1 (devAge) + model 2 (predAge)
BOLT_PATH="/fast/software/BOLT-LMM/BOLT-LMM_v2.4.1/" # path to BOLT-LMM
PHENO_PATH="$(pwd)/results/output-backup" # Path to 'pheno_covars.txt'
for PHENO in devAge predAge; do
	bolt \
		--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
		--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
		--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
		--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
		--LDscoresMatchBp \
		--lmmForceNonInf \
		--noBgenIDcheck \
		--covarUseMissingIndic \
		--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
		--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
		--numThreads=56 \
		--lmm \
		--phenoFile=${PHENO_PATH}/pheno_covars_discov1-2.txt \
		--phenoCol=${PHENO} \
		--covarFile=${PHENO_PATH}/pheno_covars_discov1-2.txt \
		--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
		--qCovarCol={PC{1:10},AGE,AGE2,SEX,ICV} \
		--covarMaxLevels 1000 \
		--verboseStats \
		--h2gGuess=0.10 \
		--statsFile="${targetDir}"/ukbb.discov1-2.chr1-23_non_imputed_${PHENO}.stats \
		--statsFileBgenSnps="${targetDir}"/ukbb.discov1-2.chr1-23_imputed_${PHENO}.stats \
		--numLeaveOutChunks 2 \
		--bgenMinMAF=0.01 \
		--bgenMinINFO=0.3 \
		2>&1 | tee ${targetDir}/ukbb.discov1-2.chr1-23_imputed_${PHENO}_$(date +'%Y%m%d').log
	pigz -f "${targetDir}"/ukbb.discov1-2.chr1-23_non_imputed_${PHENO}.stats
	pigz -f "${targetDir}"/ukbb.discov1-2.chr1-23_imputed_${PHENO}.stats
done

# model 3 (noICV)
bolt \
	--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
	--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
	--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
	--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
	--LDscoresMatchBp \
	--lmmForceNonInf \
	--noBgenIDcheck \
	--covarUseMissingIndic \
	--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
	--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
	--numThreads=56 \
	--lmm \
	--phenoFile=${PHENO_PATH}/pheno_covars_discov1-2.txt \
	--phenoCol=devAge \
	--covarFile=${PHENO_PATH}/pheno_covars_discov1-2.txt \
	--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
	--qCovarCol={PC{1:10},AGE,AGE2,SEX} \
	--covarMaxLevels 1000 \
	--verboseStats \
	--h2gGuess=0.10 \
	--statsFile="${targetDir}"/ukbb.discov1-2.chr1-23_non_imputed_devAge_noICV.stats \
	--statsFileBgenSnps="${targetDir}"/ukbb.discov1-2.chr1-23_imputed_devAge_noICV.stats \
	--numLeaveOutChunks 2 \
	--bgenMinMAF=0.01 \
	--bgenMinINFO=0.3 \
	2>&1 | tee ${targetDir}/ukbb.discov1-2.chr1-23_imputed_devAge_noICV_$(date +'%Y%m%d').log
pigz -f "${targetDir}"/ukbb.discov1-2.chr1-23_non_imputed_devAge_noICV.stats
pigz -f "${targetDir}"/ukbb.discov1-2.chr1-23_imputed_devAge_noICV.stats

# model 4 + 5 (sex-stratified)
for SEX in male female; do
	bolt \
		--bed=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bed \
		--bim=${targetDir}/tmp/chr{1:23}_mri_qc_pruned.bim \
		--fam=${targetDir}/tmp/chr1_mri_qc_pruned.fam \
		--bgenSampleFileList=${targetDir}/tmp/bolt.filelist \
		--LDscoresMatchBp \
		--lmmForceNonInf \
		--noBgenIDcheck \
		--covarUseMissingIndic \
		--LDscoresFile=${BOLT_PATH}/tables/LDSCORE.1000G_EUR.tab.gz \
		--geneticMapFile=${BOLT_PATH}/tables/genetic_map_hg19_withX.txt.gz \
		--numThreads=56 \
		--lmm \
		--phenoFile=${PHENO_PATH}/pheno_covars_discov1-2_${SEX}.txt \
		--phenoCol=devAge \
		--covarFile=${PHENO_PATH}/pheno_covars_discov1-2_${SEX}.txt \
		--covarCol={SITE.f2,SITE.f3,SITE.f4,array} \
		--qCovarCol={PC{1:10},AGE,AGE2,ICV} \
		--covarMaxLevels 1000 \
		--verboseStats \
		--h2gGuess=0.10 \
		--statsFile="${targetDir}"/ukbb.discov1-2.chr1-23_non_imputed_devAge_${SEX}.stats \
		--statsFileBgenSnps="${targetDir}"/ukbb.discov1-2.chr1-23_imputed_devAge_${SEX}.stats \
		--numLeaveOutChunks 2 \
		--bgenMinMAF=0.01 \
		--bgenMinINFO=0.3 \
		2>&1 | tee ${targetDir}/ukbb.discov1-2.chr1-23_imputed_devAge_${SEX}_$(date +'%Y%m%d').log
	pigz -f "${targetDir}"/ukbb.discov1-2.chr1-23_non_imputed_devAge_${SEX}.stats
	pigz -f "${targetDir}"/ukbb.discov1-2.chr1-23_imputed_devAge_${SEX}.stats
done
chmod 770 ${targetDir}/*
rm -rf ${targetDir}/tmp

# ==============================================================================================================
# === prepare Jawinski et al. summary statistics for Genomic SEM (excluding 2k individuals for pgs analysis) ===
# ==============================================================================================================

# set working directory
cd /slow/projects/enigma_brainage
conda activate envs/default

# prepare summary statistics for zenodo
targetDir="data/brainage_sumstats/"
mkdir -p "${targetDir}"
for tissue in wm gm gwm; do
	echo "Starting with gap_${tissue}"

	# discovery & replication (excluding 2k individual)
	awk '$2!="MT"' <(pigz -dc "/slow/projects/ukb_brainage/results/gap_${tissue}/gwama/eur/eur.excl.2k/metal.ivweight.qc.gz") > "${targetDir}/brainage2025.full.eur.excl2k.${tissue}"
	chmod 770 "${targetDir}/brainage2025.full.eur.excl2k.${tissue}"
	pigz -f "${targetDir}/brainage2025.full.eur.excl2k.${tissue}"

done

# upload file
files=$(ls ${targetDir}/brainage2025.full.eur.excl2k.*)
hubox FILELINK "${files}"

# ========================================================
# === create phenotypic plots for selected individuals === 
# ========================================================

# needs to be added

# ====================
# === PGS analysis === 
# ====================

# set working directory
cd /slow/projects/enigma_brainage
conda activate envs/default

# Harmonize brainage sumstats
mkdir -p data/brainage_sumstats/harmonized/ data/brainage_sumstats/six-brainage-sumstats/harmonized/

	# brainageFactor (original sumstats contain 1000 genomes frequencies, recover A1 frequencies)
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/brainage_factor_nogenr_excl2k_jawinski_noGC.txt.gz" \
		--outFile "data/brainage_sumstats/harmonized/brainageFactor" \
		--colsIn "SNP,CHR,BP,A1,A2,est,se_c,Pval_Estimate,N" \
		--colsRename "ID,CHR,BP,A1,A2,BETA,SE,P,N" \
		--colsOut "ID,CHR,BP,A1,A2,KGP_A1_FREQ,BETA,SE,P,N" \
		--addRSID TRUE \
		--kgpFiles 'data/1kgp/v5a/ALL.chr\\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz' \
		--keepKGPfreq TRUE

	# ENIGMA
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/six-brainage-sumstats/METAANALYSIS_ENIGMA_combinedUKBB_GCon_only-samplesize_info0.6_2025-05-27.txt.gz" \
		--outFile "data/brainage_sumstats/six-brainage-sumstats/harmonized/enigma" \
		--colsIn "MarkerName,CHR,BP,Allele1,Allele2,Freq1,Zscore,P,N" \
		--colsRename "ID,CHR,BP,A1,A2,A1_FREQ,Z,P,N" \
		--colsOut "ID,CHR,BP,A1,A2,A1_FREQ,BETA,SE,P,N" \
		--compMAF TRUE \
		--compBetaSe TRUE

	# Kaufmann et al.
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/six-brainage-sumstats/Kaufman_sumstats.gz" \
		--outFile "data/brainage_sumstats/six-brainage-sumstats/harmonized/kaufmann" \
		--colsIn "SNP,CHR,BP,A1,A2,BETA,SE,PVAL" \
		--colsRename "ID,CHR,BP,A1,A2,BETA,SE,P" \
		--colsOut "ID,CHR,BP,A1,A2,BETA,SE,P,N" \
		--setN 20170
		awk 'NR==FNR { snp[$1":"$5":"$6":"$7]=$8; next }
			 FNR==1 { print $0, "A1_FREQ"; next }
			 $2":"$3":"$4":"$5 in snp { print $0, snp[$2":"$3":"$4":"$5]; next }
			 $2":"$3":"$5":"$4 in snp { print $0, 1-snp[$2":"$3":"$5":"$4]; next }' OFS='\t' /fast/software/gctb/resources/ukbEUR_Imputed/snp.info <(gzip -dc "data/brainage_sumstats/six-brainage-sumstats/harmonized/kaufmann.gz") \
			 > data/brainage_sumstats/six-brainage-sumstats/harmonized/kaufmann_freqUKB
		pigz data/brainage_sumstats/six-brainage-sumstats/harmonized/kaufmann_freqUKB
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/six-brainage-sumstats/Kaufman_sumstats.gz" \
		--outFile "data/brainage_sumstats/six-brainage-sumstats/harmonized/kaufmann_freq1KGP" \
		--colsIn "SNP,CHR,BP,A1,A2,BETA,SE,PVAL" \
		--colsRename "ID,CHR,BP,A1,A2,BETA,SE,P" \
		--colsOut "ID,CHR,BP,A1,A2,KGP_A1_FREQ,BETA,SE,P,N" \
		--setN 20170 \
		--addRSID TRUE \
		--kgpFiles 'data/1kgp/v5a/ALL.chr\\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz' \
		--keepKGPfreq TRUE

	# Leonardsen et al.
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/six-brainage-sumstats/leonardsen.txt.gz" \
		--outFile "data/brainage_sumstats/six-brainage-sumstats/harmonized/leonardsen" \
		--colsIn "SNP,CHR,BP,A1,A2,FRQ,BETA,SE,P,N" \
		--colsRename "ID,CHR,BP,A1,A2,A1_FREQ,BETA,SE,P,N" \
		--colsOut "ID,CHR,BP,A1,A2,A1_FREQ,BETA,SE,P,N"
	zcat data/brainage_sumstats/six-brainage-sumstats/harmonized/leonardsen.gz | awk 'NR==1 { print; next } $7==0 && $8==0 { $8=0.01 } { print }' > data/brainage_sumstats/six-brainage-sumstats/harmonized/leonardsen
	pigz -f data/brainage_sumstats/six-brainage-sumstats/harmonized/leonardsen

	# Smith et al.
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/six-brainage-sumstats/Smith_V0140_sumstats_A2swap.txt.gz" \
		--outFile "data/brainage_sumstats/six-brainage-sumstats/harmonized/smith" \
		--colsIn "rsid,chr,pos,a1,a2,af,beta,SE,P" \
		--colsRename "ID,CHR,BP,A1,A2,A1_FREQ,BETA,SE,P" \
		--colsOut "ID,CHR,BP,A1,A2,A1_FREQ,BETA,SE,P,N" \
		--setN 10612

	# Wen et al.
	Rscript code/sumstats.harmonize.R \
		--sumstats "data/brainage_sumstats/six-brainage-sumstats/oag_pheno_normalized_residualized.Brain_age_gap.glm.linear.gz" \
		--outFile "data/brainage_sumstats/six-brainage-sumstats/harmonized/wen" \
		--colsIn "ID,CHR,BP,A1,A2,BETA,SE,P,OBS_CT" \
		--colsRename "ID,CHR,BP,A1,A2,BETA,SE,P,N" \
		--colsOut "ID,CHR,BP,A1,A2,BETA,SE,P,N"
		awk 'NR==FNR { snp[$1":"$5":"$6":"$7]=$8; next }
			 FNR==1 { print $0, "A1_FREQ"; next }
			 $2":"$3":"$4":"$5 in snp { print $0, snp[$2":"$3":"$4":"$5]; next }
			 $2":"$3":"$5":"$4 in snp { print $0, 1-snp[$2":"$3":"$5":"$4]; next }' OFS='\t' /fast/software/gctb/resources/ukbEUR_Imputed/snp.info <(gzip -dc "data/brainage_sumstats/six-brainage-sumstats/harmonized/wen.gz") \
			 > data/brainage_sumstats/six-brainage-sumstats/harmonized/wen_freqUKB
		pigz data/brainage_sumstats/six-brainage-sumstats/harmonized/wen_freqUKB

# get SBayesRC weights
mkdir -p results/pgsweights
trait=(brainageFactor enigma jawinski kaufmann leonardsen smith wen)
sumstats=(data/brainage_sumstats/harmonized/brainageFactor.gz data/brainage_sumstats/six-brainage-sumstats/harmonized/enigma.gz data/brainage_sumstats/brainage2025.full.eur.excl2k.gwm.gz data/brainage_sumstats/six-brainage-sumstats/harmonized/kaufmann_freqUKB.gz data/brainage_sumstats/six-brainage-sumstats/harmonized/leonardsen.gz data/brainage_sumstats/six-brainage-sumstats/harmonized/smith.gz data/brainage_sumstats/six-brainage-sumstats/harmonized/wen_freqUKB.gz)
cols=("ID,A1,A2,KGP_A1_FREQ,BETA,SE,P,N" "ID,A1,A2,A1_FREQ,BETA,SE,P,N" "ID,A1,A2,A1_FREQ,BETA,SE,P,N" "ID,A1,A2,A1_FREQ,BETA,SE,P,N" "ID,A1,A2,A1_FREQ,BETA,SE,P,N" "ID,A1,A2,A1_FREQ,BETA,SE,P,N" "ID,A1,A2,A1_FREQ,BETA,SE,P,N") # set input column names that correspond to header columns SNP A1 A2 freq b se p N
ldmFolder="/fast/software/gctb/resources/ukbEUR_Imputed" # ldmFolder="/lustre/psychologie/jawinskp/software/gctb/resources/ukbEUR_Imputed" 
annotFile="/fast/software/gctb/resources/annot_baseline2.2.txt" # annotFile="/lustre/psychologie/jawinskp/software/gctb/resources/annot_baseline2.2.txt"
threads=100
sbayes="sbayesrc"
imputation=1
mcmc=1
for ((i=0; i<${#trait[@]}; i++)); do
	targetDir="results/pgsweights/${trait[i]}"
	mkdir -p ${targetDir}/deprecated
	mv ${targetDir}/* ${targetDir}/deprecated
	./code/sbayes.sh "${targetDir}" "${sumstats[i]}" "${cols[i]}" "${ldmFolder}" "${annotFile}" "${threads}" "${sbayes}" "${imputation}" "${mcmc}" # taskset -c "${core}" d
done

# calculate polygenic scores in held-out samples
pgenFileHandlerBeta='data/genetics/chr${i}/imp_mri_qc_ANCESTRY/chr${i}_mri_qc'
idCol="Name"
a1Col="A1"
effectCol="A1Effect"
maf=0.01
threads=100
mkdir -p results/pgs
for ancestry in AFR CSA EAS EUR; do # do not use AMR and MID due to low sample sizes (n < 100)
	for trait in enigma; do # brainageFactor enigma jawinski kaufmann leonardsen smith wen
		pgenFileHandler=$(echo ${pgenFileHandlerBeta} | sed "s/ANCESTRY/${ancestry}/g")
		outFile="results/pgs/pgs.sbayesrc.${trait}.${ancestry}"
		weightFile="results/pgsweights/${trait}/sbayesrc.snpRes.gz"
		code/pgs.predict.sh "${weightFile}" "${pgenFileHandler}" "${idCol}" "${a1Col}" "${effectCol}" "${maf}" "${threads}" "${outFile}"
	done
done

# add photon predictions to phenotype files
Rscript -e "\
load('results/output-backup/data_pheno_covariates.Rdata'); \
repl = read.delim('/slow/projects/ukb_brainage/data/traits/replicate.2k.txt', header = T, sep = '\t');
repl = dplyr::left_join(repl,data[,c('SUBJID','predAge','devAge')], by = dplyr::join_by(FID == SUBJID));
write.table(repl, 'results/pgs/pgs.phenotypes.eur.txt', row.names = FALSE, quote = FALSE, sep = '\t')"
Rscript -e "\
load('results/output-backup/data_pheno_covariates.Rdata'); \
repl = read.delim('/slow/projects/ukb_brainage/data/traits/replicate.txt', header = T, sep = '\t');
repl = dplyr::left_join(repl,data[,c('SUBJID','predAge','devAge')], by = dplyr::join_by(FID == SUBJID));
write.table(repl, 'results/pgs/pgs.phenotypes.all.txt', row.names = FALSE, quote = FALSE, sep = '\t')"

# correlate polygenic scores with phenotypes in replication sample - non-European ancestry individuals
joinVar="IID"
pgsVar="SCORE1_SUM"
covs="sex,age,age2,ac1,ac2,ac3,TIV,array,PanC1,PanC2,PanC3,PanC4,PanC5,PanC6,PanC7,PanC8,PanC9,PanC10"
phenoVar=(devAge gap_gwm)
phenoLabel=(BAGenigma BAGjawinski)
for ancestry in EUR AFR CSA EAS; do
	for pgs in enigma jawinski kaufmann leonardsen smith wen; do
		for ((i=0; i<${#phenoVar[@]}; i++)); do
			if [[ "$ancestry" == "EUR" ]]; then
			    phenoFile="results/pgs/pgs.phenotypes.eur.txt"
			else
			    phenoFile="results/pgs/pgs.phenotypes.all.txt"
			fi
			pgsFile="results/pgs/pgs.sbayesrc.${pgs}.${ancestry}.score" # results/brainageFactor/pgs/pgs.sbayesrc.hm3.EUR.score
			outFile="results/pgs/pgs.sbayesrc.${pgs}.${ancestry}.${phenoLabel[i]}.assoc" # results/brainageFactor/pgs/pgs.sbayesrc.hm3.EUR.assoc
			code/pgs.corr.R "${pgs}" "${pgsFile}" "${phenoFile}" "${joinVar}" "${pgsVar}" "${phenoVar[i]}" "${covs}" "${outFile}"
		done
	done	
done

# join PRS association results
for phenotype in BAGenigma BAGjawinski; do
	awk -v phenotype="${phenotype}" 'NR==1 { output="ancestry\t"$1; for (i=2; i<=NF; i++) { output=output"\t"phenotype"_"$i }; print output }' OFS='\t' "results/pgs/pgs.sbayesrc.brainageFactor.EUR.BAGenigma.assoc" > "results/pgs/pgs.${phenotype}.assoc.txt"
	for ancestry in EUR AFR CSA EAS; do
		for pgs in brainageFactor enigma jawinski kaufmann leonardsen smith wen; do
			awk -v ancestry="${ancestry}" '
				FNR==1 { next }
				{ print ancestry, $0 }' OFS='\t' "results/pgs/pgs.sbayesrc.${pgs}.${ancestry}.${phenotype}.assoc" \
				>> "results/pgs/pgs.${phenotype}.assoc.txt"
		done
		echo "" >> "results/pgs/pgs.${phenotype}.assoc.txt"
	done
done
paste \
  <(awk -F'\t' '{print $1, $2, "", $4, $3, $5, $8, $7}' OFS="\t" results/pgs/pgs.BAGenigma.assoc.txt) \
  <(awk -F'\t' '{print "", $4, $3, $5, $8, $7}' OFS="\t" results/pgs/pgs.BAGjawinski.assoc.txt) \
  > results/pgs/pgs.assoc.txt

# =======================
# === PheWAS analysis === 
# =======================

# set working directory
cd /slow/projects/enigma_brainage
conda activate envs/default

# Calculate polygenic scores in full UKB sample
fileHandler='data/genetics/chr${i}/imp/pgen/ukb_imp_chr${i}_v3'
fileType="pgen"
idCol="Name"
a1Col="A1"
effectCol="A1Effect"
maf=0.00
threads=50
for pgs in brainageFactor enigma jawinski kaufmann leonardsen smith wen; do
	outFile="results/phewas/sample/pgs.sbayesrc.${pgs}"
	weightFile="results/pgsweights/${pgs}/sbayesrc.snpRes.gz"
	taskset -c 56-111 code/pgs.predict.sh "${weightFile}" "${fileHandler}" "${fileType}" "${idCol}" "${a1Col}" "${effectCol}" "${maf}" "${threads}" "${outFile}"
done

# Prepare phenotype file for individuals without relatedness to the MRI sample
csvFile="data/basket/20240307_4017567/data/ukb678162.csv"
imagingFile="data/basket/dnanexus/fields.imaging.20250407.csv"
occurenceFile="data/basket/dnanexus/fields.first.occurence.output.txt.csv.gz"
kinshipFile="data/genetics/2024/meta/ukb42032_rel_s487950.dat"
panFile="data/basket/20210327_2008261/data/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
bridgeFile="data/basket/20210327_2008261/data/ukb42032bridge31063.txt"
targetDir="results/phewas/sample/"
phesantVarFile="/fast/software/PHESANT/variable-info/outcome-info.tsv"
Rscript code/sampleFiltering.R "${csvFile}" "${imagingFile}" "${occurenceFile}" "${kinshipFile}" "${panFile}" "${bridgeFile}" "${targetDir}"

# Run PHESANT on the non-MRI UK Biobank sample using the SLURM cluster (slurm-login.hpc-service.hu-berlin.de)
# Note: One job file was created and submitted per polygenic score to enable parallel execution across multiple jobs (see code/slurm/)

	#!/bin/bash
	#SBATCH --job-name=phesant
	#SBATCH --nodes=1
	#SBATCH --ntasks=1
	#SBATCH --cpus-per-task=48
	#SBATCH --mem=194300M
	#SBATCH --time=4-00:00:00
	#SBATCH --output=logs/phesant.out
	#SBATCH --error=logs/phesant.err
	#SBATCH --partition=gpu
	#SBATCH --gres=gpu

	# Load environment
	cd /lustre/psychologie/jawinskp/projects/enigma_brainage
	source ~/.bashrc
	conda activate envs/phesant

	# settings
	covsFile="results/phewas/sample/r2025.vars.discovery"
	covs="sex,age,array,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
	phenoFile="results/phewas/sample/ukb678162.phesant.2025.wba.csv"
	phesantDir="/lustre/psychologie/jawinskp/software/PHESANT"
	nparts=25
	partList="$(echo {1..25})"
	standardise="TRUE"
	genetic="FALSE"
	sensitivity="FALSE"

	# Run PHESANT for each trait
	for trait in brainageFactor enigma jawinski kaufmann leonardsen smith wen; do
	    targetDir="results/phewas/${trait}/"
	    rm -rf "${targetDir}"
	    mkdir -p "${targetDir}"
	    traitFile="results/phewas/sample/pgs.sbayesrc.${trait}.score"
	    code/phesant.sh "SCORE1_SUM" "${traitFile}" "${covsFile}" "${covs}" \
	        "${targetDir}" "${phenoFile}" "${phesantDir}" "${nparts}" "${partList}" \
	        "${standardise}" "${genetic}" "${sensitivity}" \
	        2>&1 | tee "${targetDir}/phesant.log"
	done

# create plots and get result summary
ncovs=0 # do not adjust degrees of freedom (not super accurate due to the different models used)
imaging=FALSE # show imaging variables
multipleTesting=both # both FDR and Bonferroni lines shown
ylim=20 # y axis limit
ybreaks=5 # y breaks
width=8.0 # plot width
height=6.5 # plot height
repel_nudge_y=10 # annotations (have been commented out in the current script)
desat=FALSE # desatuation of non-significant results
(for trait in brainageFactor enigma jawinski kaufmann leonardsen smith wen; do (
	
	# standard results
	targetDir="results/phewas/${trait}/"
	phesantResults="results/phewas/${trait}/phesant.output/results-combined.txt"
	code/phesant.output.R "${trait}" "${phesantResults}" "${ncovs}" "${targetDir}" "${imaging}" "${multipleTesting}" "${ylim}" "${ybreaks}" "${width}" "${height}" "${repel_nudge_y}" "${desat}"

	# with merged instances to increase sample size
	targetDir="results/phewas/mergedInstances/${trait}/"
	phesantResults="results/phewas/mergedInstances/${trait}/phesant.output/results-combined.txt"
	code/phesant.output.R "${trait}" "${phesantResults}" "${ncovs}" "${targetDir}" "${imaging}" "${multipleTesting}" "${ylim}" "${ybreaks}" "${width}" "${height}" "${repel_nudge_y}" "${desat}"
	) &
done)
wait

# combine phesant plots and result summaries across traits - standard results
traits="brainageFactor,enigma,jawinski,kaufmann,leonardsen,smith,wen"
phesantSummary=$(IFS=,; for trait in $traits; do echo -n "results/phewas/${trait}/phesant.summary.txt,"; done | sed 's/,$//') # phesantSummary="results/phewas/brainageFactor/phesant.summary.txt,results/phewas/enigma/phesant.summary.txt,results/phewas/jawinski/phesant.summary.txt,results/phewas/kaufmann/phesant.summary.txt,results/phewas/leonardsen/phesant.summary.txt,results/phewas/smith/phesant.summary.txt,results/phewas/wen/phesant.summary.txt"
phesantPlot="" # $(IFS=,; for trait in $traits; do echo -n "results/phewas/${trait}/phesant.png,"; done | sed 's/,$//')
outputFile="results/phewas/phewas"
code/phesant.combine.R "${traits}" "${phesantSummary}" "${phesantPlot}" "${outputFile}"

# combine phesant plots and result summaries across traits - merged instances
traits="brainageFactor,enigma,jawinski,kaufmann,leonardsen,smith,wen"
phesantSummary=$(IFS=,; for trait in $traits; do echo -n "results/phewas/mergedInstances/${trait}/phesant.summary.txt,"; done | sed 's/,$//') # phesantSummary="results/phewas/brainageFactor/phesant.summary.txt,results/phewas/enigma/phesant.summary.txt,results/phewas/jawinski/phesant.summary.txt,results/phewas/kaufmann/phesant.summary.txt,results/phewas/leonardsen/phesant.summary.txt,results/phewas/smith/phesant.summary.txt,results/phewas/wen/phesant.summary.txt"
phesantPlot="" # $(IFS=,; for trait in $traits; do echo -n "results/phewas/${trait}/phesant.png,"; done | sed 's/,$//')
outputFile="results/phewas/mergedInstances/phewas"
code/phesant.combine.R "${traits}" "${phesantSummary}" "${phesantPlot}" "${outputFile}"


