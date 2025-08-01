#!/usr/bin/env Rscript

# ====================================
# === Harmonize summary statistics ===
# ====================================

# load required package
library(argparse)

# Create parser
parser = ArgumentParser(description = 'Harmonize summary statistics')

# Define all expected arguments
parser$add_argument('--sumstats', metavar='[FILE]', required=TRUE, help='Path to summary statistics file.')
parser$add_argument('--outFile', metavar='[FILE]', required=TRUE, help='Output file path.')
parser$add_argument('--skipLines', metavar='[val]', type='integer', default=0, help='Number of lines to skip in the input file.')
parser$add_argument('--sep', metavar='[delimiter]', default='auto', help='Field separator in the input file.')
parser$add_argument('--colsIn', metavar='[column1,column2,...]', required=TRUE, help='Input column names as comma-separated string.')
parser$add_argument('--colsRename', metavar='[column1,column2,...]', required=TRUE, help='Column renaming mapping as comma-separated string.')
parser$add_argument('--colsOut', metavar='[column1,column2,...]', required=TRUE, help='Output column names as comma-separated string.')
parser$add_argument('--halveNeff', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Halve the Neff column values (required renaming of Neff to Neff_half).')
parser$add_argument('--setNca', metavar='[val]', type='double', default=-1, help='Set number of cases.')
parser$add_argument('--setNco', metavar='[val]', type='double', default=-1, help='Set number of controls.')
parser$add_argument('--setN', metavar='[val]', type='double', default=-1, help='Set total sample size.')
parser$add_argument('--compNeff', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Compute Neff from cases and controls.')
parser$add_argument('--compNeff_half', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Compute Neff_half by dividing Neff by 2.')
parser$add_argument('--compN', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Compute total N from cases and controls.')
parser$add_argument('--compMAF', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Compute minor allele frequency.')
parser$add_argument('--invertCol', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Invert values in a specified column.')
parser$add_argument('--rmDuplicatesCol', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Remove duplicates by column.')
parser$add_argument('--or2beta', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Convert odds ratio to beta.')
parser$add_argument('--compBetaSe', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Compute Beta and SE (requires Z and MAF).')
parser$add_argument('--stdize', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Standardize betas to unit phenotypic variance')
parser$add_argument('--addRSID', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Add RSID from reference.')
parser$add_argument('--kgpFiles', metavar='[FILEPATH]', help="Path to 1000 Genomes reference file pattern, e.g. 'data/1kgp/v5a/ALL.chr\\\\$chr.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.meta.vcf.gz'")
parser$add_argument('--keepKGPfreq', metavar='[TRUE/FALSE]', type='logical', default=FALSE, help='Retain allele frequencies from KGP.')

# Parse arguments and assign them into the global environment
args = parser$parse_args()
list2env(as.list(args), envir = globalenv())
logInfo = paste0('\n--- Harmonize summary statistics ---',
               '\nsumstats: ', sumstats,
               '\noutFile: ', outFile,
               '\nskipLines: ', skipLines,
               '\nsep: ', sep,
               '\ncolsIn: ', colsIn,
               '\ncolsRename: ', colsRename,
               '\ncolsOut: ', colsOut,
               '\nhalveNeff: ', halveNeff,
               '\nsetNca: ', setNca,
               '\nsetNco: ', setNco,
               '\nsetN: ', setN,
               '\ncompNeff: ', compNeff,
               '\ncompNeff_half: ', compNeff_half,
               '\ncompN: ', compN,
               '\ncompMAF: ', compMAF,
               '\ninvertCol: ', invertCol,
               '\nrmDuplicatesCol: ', rmDuplicatesCol,
               '\nor2beta: ', or2beta,
               '\ncompBetaSe: ', compBetaSe,  
               '\nstdize: ', stdize,                                                                                                                          
               '\naddRSID: ', addRSID,
               '\nkgpFiles: ', kgpFiles,
               '\nkeepKGPfreq: ', keepKGPfreq, '\n')
message(logInfo)

# transform variables
colsIn = stringr::str_split(colsIn, ',')[[1]]
colsRename = stringr::str_split(colsRename, ',')[[1]]
colsOut = stringr::str_split(colsOut, ',')[[1]]
if (sep == "tab") { sep = "\t" }

# import data
message(paste0(' - importing data.'))
if (stringr::str_sub(sumstats,-3,-1) == '.gz') {
  df = data.frame(data.table::fread(cmd=sprintf("gzip -dc %s", sumstats), select = colsIn, fill = T, tmpdir = getwd(), header=T, skip = skipLines, sep = sep, stringsAsFactors=FALSE))
} else {
  df = data.frame(data.table::fread(file = sumstats, select = colsIn, fill = T, tmpdir = getwd(), header=T, skip = skipLines, sep = sep, stringsAsFactors=FALSE))
}

# rename columns
message(paste0(' - renaming columns.'))
names(df) = colsRename

# remove rows with NA fields
message(' - removing rows with NA fields.')
idx = rowSums(is.na(df)) > 0
logInfo = paste0(logInfo,
  '\nrows removed due to NA: ', sum(idx))
df = df[!idx,]

# remove 'chr' in Chromosome column
if ('CHR' %in% colsRename) {
  message(" - removing 'chr' in chromsome column.")
  df$CHR = stringr::str_remove(df$CHR,'chr')
} else if ('CHR_BP' %in% colsRename) {
  message(" - removing 'chr' in chromsome column.")
  df$CHR_BP = stringr::str_remove(df$CHR_BP,'chr')
}

# use upper case for alleles
message(" - use upper case for alleles.")
if ('A1' %in% colsRename) { df$A1 = toupper(df$A1) }
if ('A2' %in% colsRename) { df$A2 = toupper(df$A2) }

# divide chr bp
if ('CHR_BP' %in% colsRename) {
  message(" - dividing CHR_BP into CHR and BP column.")
  CHR_BP = stringr::str_split_fixed(df$CHR_BP,':',4)
  df$CHR = CHR_BP[,1]
  df$BP = CHR_BP[,2]
  df = df[,c('CHR','BP',colsRename[!(colsRename %in% 'CHR_BP')])]
}

# make sure that P and BP column is numeric
if('P' %in% names(df)) { df$P = as.numeric(df$P) }
if('BP' %in% names(df)) { df$BP = as.numeric(df$BP) }

# extract rsID from rs:a1:a2 identifiers
if('ID' %in% names(df) & length(grep(':',df$ID)) > 1 & length(stringr::str_starts(df$ID,'rs')) > 1) {
    message(" - extracting rsid from ID column (rs:a1:a2).")
    idx = stringr::str_starts(df$ID,'rs')
    rsid = stringr::str_split_fixed(df$ID[idx],':',2)
    df$ID[idx] = rsid[,1]
}

# add N, Nca, Nco | compute Neff and N (requires Nca and Nco)
message(' - adding N / computing Neff / computing N.')
if (halveNeff == TRUE ) { df$Neff_half = df$Neff_half/2 }
if (setNca > 0 ) { df$Nca = setNca }
if (setNco > 0 ) { df$Nco = setNco }
if (setN > 0 ) { df$N = setN }
if (compNeff == TRUE) { df$Neff = 4/(1/df$Nca + 1/df$Nco) }
if (compNeff_half == TRUE) { df$Neff_half = 2/(1/df$Nca + 1/df$Nco) } # same as 4*Nca*Nco/(2*(Nca+Nco)); do not confuse with Neff, which is defined as 4/(1/Nca + 1/Nco)
if (compN == TRUE) { df$N = df$Nca + df$Nco }

# compute MAF (requires A1_FREQ)
if (compMAF == TRUE) { 
  message(' - computing MAF.')
  df$MAF = df$A1_FREQ
  df$MAF[df$A1_FREQ > 0.5] = 1-df$A1_FREQ[df$A1_FREQ > 0.5]
}

# invert effect column
if (invertCol != FALSE) { 
  message(' - inverting effect column.')
  df[,invertCol] = -df[,invertCol]
}

# remove duplicates
if (rmDuplicatesCol != FALSE) { 
  message(' - removing duplicates (first occurence will be kept).')
  idx = duplicated(df[,rmDuplicatesCol])
  logInfo = paste0(logInfo,
    '\nduplicates removed: ', sum(idx))
  df = df[!idx,]
}

# recode odds ratio to beta
if (or2beta == TRUE) { 
  message(' - converting odds ratio to beta coeffcient.')
  df$BETA = log(df$OR)
}

# add rsid
if (addRSID == TRUE) {

  # run shell sript
  message(' - adding 1000 Genomes Project rsID.')
  wd = funr::get_script_path() # wd = "code/genetics"
  data.table::fwrite(df[,c('CHR','BP','A1','A2')], file = sprintf('%s.tmp.txt', outFile), sep = '\t')
  system(sprintf('%s/sumstats.harmonize.rsid.sh %s.tmp.txt %s', wd, outFile, kgpFiles))

  # import results
  message('\n - importing rsIDs of matched variants.')
  kgp = data.table::fread(sprintf('%s.tmp.txt', outFile))
  kgp$CHR = as.character(kgp$CHR)
  df = dplyr::left_join(df,kgp, by = c('CHR','BP','A1','A2'))

  # how many rows will be removed due to NA?
  idx = is.na(df$KGP_REF)
  logInfo = paste0(logInfo,
    '\nrows removed due to non-existent 1000 genomes match: ', sum(idx),
    '\nrows remaining: ', sum(!idx))
  df = df[!idx,]

  # keep allele frequencies?
  if (keepKGPfreq == TRUE) {
    message(' - adding 1000 Genomes A1 allele frequency and MAF column.')
    df$KGP_A1_FREQ = df$KGP_EUR_AF
    idx = df$A1 == df$KGP_REF & !is.na(df$KGP_REF)
    df$KGP_A1_FREQ[idx] = 1-df$KGP_EUR_AF[idx]
    df$KGP_MAF = df$KGP_A1_FREQ
    idx = df$KGP_A1_FREQ > 0.5 & !is.na(df$KGP_A1_FREQ)
    df$KGP_MAF[idx] = 1-df$KGP_MAF[idx]
  }
  
  df = df[,c('CHR','BP','KGP_RSID',names(df)[!(names(df) %in% c('CHR','BP','KGP_RSID','KGP_REF','KGP_ALT','KGP_EUR_AF'))])]
  names(df)[3] = 'ID'
  system(sprintf('rm -f %s.tmp.txt', outFile))
}

# compute BETA and SE from Z and MAF
if (compBetaSe == TRUE) { 
  message(' - calculating BETA and SE from Z and MAF.')

  if ('A1_FREQ' %in% names(df)) {
    MAF = df$A1_FREQ
    MAF[df$A1_FREQ > 0.5] = 1-df$A1_FREQ[df$A1_FREQ > 0.5] 
  } else { 
    MAF = df$KGP_A1_FREQ
    MAF[df$KGP_A1_FREQ > 0.5] = 1-df$KGP_A1_FREQ[df$KGP_A1_FREQ > 0.5] 
  }
  df$BETA = (df$Z)/sqrt(2*MAF*(1-MAF)*(df$N+(df$Z)^2))
  df$SE = 1/sqrt(2*MAF*(1-MAF)*(df$N+(df$Z)^2))
}

# standardize
if (stdize == TRUE && sum(c('BETA','SE','A1_FREQ','N') %in% names(df)) == 4) {
  message(' - standardizing BETA and SE.')
  A1_FREQ = df$A1_FREQ
  A1_FREQ[A1_FREQ == 1] = max(A1_FREQ[A1_FREQ!=1])
  A1_FREQ[A1_FREQ == 0] = min(A1_FREQ[A1_FREQ!=0])
  BETA = (df$BETA/df$SE)/sqrt(2*A1_FREQ*(1-A1_FREQ)*(df$N+(df$BETA/df$SE)^2))
  SE = 1/sqrt(2*A1_FREQ*(1-A1_FREQ)*(df$N+(df$BETA/df$SE)^2))
  df$BETA = BETA
  df$SE = SE
}

# write output
message(sprintf(' - writing %s.gz', outFile))
data.table::fwrite(df[,colsOut], file = sprintf('%s.gz', outFile), sep = '\t', compress = 'gzip')
system(sprintf('chmod -R 770 %s.gz', outFile))

# save log file
message(sprintf(' - writing %s.log', outFile))
sink(sprintf('%s.log', outFile))
sink(stdout(), type = "message")
message(logInfo)
sink()
system(sprintf('chmod -R 770 %s.log', outFile))
message('--- Completed: Harmonize summary statistics ---')

