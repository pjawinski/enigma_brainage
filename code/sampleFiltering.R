#!/usr/bin/env Rscript

# ========================
# === Sample Filtering === 
# ========================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=7) {
  stop(paste0('expected 7 arguments, but ', length(args), ' arguments provided.'), call.=FALSE)
}

# get arguments from command line
csvFile = args[1] # csvFile="data/basket/20240307_4017567/data/ukb678162.csv"
imagingFile = args[2] # imagingFile="data/basket/dnanexus/fields.imaging.20250407.csv"
occurenceFile = args[3] # occurenceFile="data/basket/dnanexus/fields.first.occurence.output.txt.csv.gz"
kinshipFile = args[4] # kinshipFile="data/genetics/2024/meta/ukb42032_rel_s487950.dat"
panFile= args[5] # panFile="data/basket/20210327_2008261/data/Files for retman/all_pops_non_eur_pruned_within_pop_pc_covs.tsv"
bridgeFile = args[6] # bridgeFile="data/basket/20210327_2008261/data/ukb42032bridge31063.txt"
phesantVarFile = args[7] # phesantVarFile="/fast/software/PHESANT/variable-info/outcome-info.tsv"
targetDir = args[8] # targetDir="results/sample/"

logInfo = paste0('\n--- Sample Filtering: Settings ---',
               '\ncsvFile: ', csvFile,
               '\nimagingFile: ', imagingFile,
               '\noccurenceFile: ', occurenceFile,
               '\nkinshipFile: ', kinshipFile,
               '\npanFile: ', panFile,
               '\nbridgeFile: ', bridgeFile,
               '\nphesantVarFile: ', phesantVarFile,
               '\ntargetDir: ', targetDir,'\n')
message(logInfo)

# load required packages
for (pkg in c('dplyr','igraph','data.table','future')) { eval(bquote(suppressWarnings(suppressPackageStartupMessages(require(.(pkg)))))) }

# -------------------------------------------
# --- Filter data released until Feb 2024 ---
# -------------------------------------------
message('\nSelecting individuals.')

# load basket data
message(sprintf(' - loading basket file %s', csvFile))
bd = data.frame(data.table::fread(file = csvFile, header = T, sep = ','))
names(bd) = gsub("^X", "f.", names(bd))
names(bd)[1] = 'f.eid'
message(sprintf(' - %d rows and %d columns.', dim(bd)[1],dim(bd)[2]))

# apply filter criteria
# - reported and genetic sex mismatch
# - sex aneuploidy
# - outliers in heterozygosity and missingness rates
# - kinship information available
message(' - applying filter criteria.')
sex = bd$f.31.0.0==bd$f.22001.0.0; sex[is.na(sex)] = FALSE
noAneuploidy = is.na(bd$f.22019.0.0)
noHetMiss = is.na(bd$f.22027.0.0)
kinship = rep(TRUE, nrow(bd)); kinship[is.na(bd$f.22021) | bd$f.22021 == -1 ] = FALSE
bd = bd[sex & noAneuploidy & noHetMiss & kinship,]
message(sprintf(' - %d cases remaining.', dim(bd)[1]))

# load imaging dataset
message(' - loading imaging data extracted from ukb-rap (release 31/03/2025).')
imaging = data.frame(data.table::fread(imagingFile, header = T, sep = ','))
names(imaging) = gsub("^x", "f.", names(imaging))
names(imaging) = gsub("_", ".", names(imaging))
names(imaging)[1] = "f.eid"
message(sprintf(' - %d/%d imaging variables in full dataset', sum(names(imaging) %in% names(bd)), length(names(imaging))))

  # remove old imaging fields from data frame and add updated fields
  idx = which(names(bd) %in% names(imaging)[-1])
  bd = bd[,-idx]
  bd = left_join(bd,imaging,"f.eid")

# load occurence dataset
message(' - loading first occurence data extracted from ukb-rap (release 31/03/2025).')
occurence = data.frame(data.table::fread(occurenceFile, header = T, sep = ','))
names(occurence) = gsub("^X", "f.", names(occurence))
names(occurence) = gsub("_", ".", names(occurence))
names(occurence)[1] = "f.eid"
message(sprintf(' - %d/%d occurence variables in full dataset', sum(names(occurence) %in% names(bd)), length(names(occurence))))

  # remove available fields from data frame and add updated fields
  idx = which(names(bd) %in% names(occurence)[-1])
  bd = bd[,-idx]
  bd = left_join(bd,occurence,"f.eid")

# define function to select unrelated individuals
select.unrelated <- function(kinshipTable, iid, seed, keep = NULL) {

  # set random number seed
  set.seed(seed)

  # shrink kinshipTable to pairs of individuals included in sample of interest
  kinshipTable <- kinshipTable[(kinshipTable$ID1 %in% iid) & (kinshipTable$ID2 %in% iid),]

  # identify dyadic pairs
  message(' - identify and remove dyadic pair kinships.')
  cat <- c(kinshipTable$ID1,kinshipTable$ID2)
  nodeSubjects <- cat[duplicated(cat)]
  dyadic = kinshipTable[!(kinshipTable$ID1 %in% nodeSubjects) & !(kinshipTable$ID2 %in% nodeSubjects), ]

  # remove one out of two subjects in dyadic pair
  exclude = sample(1:2, nrow(dyadic), replace=T)
  dyadicExclude = c()
  pb = txtProgressBar(min = 1, max = nrow(dyadic), style = 3)
  for (i in 1:nrow(dyadic)) {
    setTxtProgressBar(pb, i)
    if (!is.null(keep) & sum(dyadic[i,c('ID1','ID2')] %in% keep) > 0) { dyadicExclude = c(dyadicExclude, dyadic[i,which(!(dyadic[i,c('ID1','ID2')] %in% keep))]) }
    else { dyadicExclude = c(dyadicExclude, dyadic[i,exclude[i]]) }
  }

  # identify subjects with multiple relationships
  multi = kinshipTable[kinshipTable$ID1 %in% nodeSubjects | kinshipTable$ID2 %in% nodeSubjects, ]
  multi$ID1 = factor(multi$ID1)
  multi$ID2 = factor(multi$ID2)
  adj = get.adjacency(graph.edgelist(as.matrix(multi[,c("ID1","ID2")]), directed = FALSE))

  # Identify trios
  message('\n - identify and remove trio kinships.')
  trios = c(NULL) # identify trios
  degrees = c(NULL)
  pb = txtProgressBar(min = 1, max = nrow(adj), style = 3)
  for (i in 1:nrow(adj)) {
  setTxtProgressBar(pb, i)
    trio = c(row.names(adj)[i])
    degree = c(sum(adj[,i]))
    
    for (j in which(adj[,i]==1)) {
      if (!(row.names(adj)[j] %in% trio) && !(row.names(adj)[j] %in% trios)) {
        trio = c(trio, row.names(adj)[j])
        degree = c(degree, sum(adj[,j]))
        if (length(trio) > 3) { break }
        
        for (k in which(adj[,j]==1)) {
          if (!(row.names(adj)[k] %in% trio) && !(row.names(adj)[k] %in% trios)) {
            trio = c(trio, row.names(adj)[k])
            degree = c(degree, sum(adj[,k]))
            if (length(trio) > 3) { break }

            for (l in which(adj[,k]==1)) {
              if (!(row.names(adj)[l] %in% trio) && !(row.names(adj)[l] %in% trios)) {
                trio = c(trio, row.names(adj)[l])
                degree = c(degree, sum(adj[,l]))
                if (length(trio) > 3) { break }
              }
            }
          }
        }
      }
    }
    if (length(trio) == 3) {
      trios = rbind(trios, trio)
      degrees = rbind(degrees, degree)
    }
  }

  # remove "hub" subject from type 1 trios (A related with B, B related with C) 
  triosType1 = matrix(trios[rowSums(degrees) == 4,], ncol = 3)
  degreesType1 = matrix(degrees[rowSums(degrees) == 4,], ncol = 3)
  triosType1Exclude = c()
  for (i in 1:nrow(triosType1)) {
    if (!is.null(keep) & sum(triosType1[i,] %in% keep) > 0 & identical(as.numeric(degreesType1[i,which(triosType1[i,] %in% keep)]),2)) { 
        triosType1Exclude = c(triosType1Exclude, triosType1[i,which(!(triosType1[i,] %in% keep))]) }
    else {
        triosType1Exclude = c(triosType1Exclude, triosType1[i, degreesType1[i,]==2]) }
  }

  # randomly remove 2 subjects from type 2 trios (a, b and c related to each other)
  triosType2 = matrix(trios[rowSums(degrees) == 6,], ncol = 3)
  select = sample(1:3, nrow(triosType2), replace=T)
  triosType2Exclude = c()
  for (i in 1:nrow(triosType2)) {
    if (!is.null(keep) & sum(triosType2[i,] %in% keep) > 0) { triosType2Exclude = c(triosType2Exclude, triosType2[i,which(!(triosType2[i] %in% keep))]) }
    else { triosType2Exclude = c(triosType2Exclude, triosType2[i, -select[i]]) }
  }

  # get list of remaining subjects
  kinshipExclude = c(dyadicExclude, triosType1Exclude, triosType2Exclude)
  iid = iid[!(iid %in% kinshipExclude)]

  # identify remaining kinships
  message('\n - identify and remove remaining kinships.')
  kinshipTable <- kinshipTable[(kinshipTable$ID1 %in% iid) & (kinshipTable$ID2 %in% iid),]
  cat = c(kinshipTable$ID1,kinshipTable$ID2)
  multi = kinshipTable[kinshipTable$ID1 %in% cat | kinshipTable$ID2 %in% cat, ]
  multi$ID1 = factor(multi$ID1)
  multi$ID2 = factor(multi$ID2)
  adj = get.adjacency(graph.edgelist(as.matrix(multi[,c("ID1","ID2")]), directed = FALSE))

  # iteratively remove subject with highest degree until no kinship remains
    # a) remove subjects related to 'keep' subjects
    if (!is.null(keep) & sum(row.names(adj) %in% keep) > 0) {
      idx = c()
      for (i in which(row.names(adj) %in% keep)) {
        idx = c(idx, as.numeric(which(adj[i,] == 1)))
      }
      idx = unique(idx)
      kinshipExclude = c(kinshipExclude,row.names(adj)[idx])
      adj = adj[-idx,-idx]
    }
    # b) remove remaining relationships
    edges = sum(adj)
    pb = txtProgressBar(min = 1, max = edges, style = 3)
    while (sum(adj) > 0) {
      setTxtProgressBar(pb, edges-sum(adj))
      idx = sample.int(nrow(adj),nrow(adj)) # randomly shuffle participants in case of ties
      adj = adj[idx,idx]
      subExclude = row.names(adj)[Matrix::rowSums(adj) %>% order(decreasing = T)][1]
      kinshipExclude = c(kinshipExclude,subExclude)
      idx = which(row.names(adj) %in% subExclude)
      adj = adj[-idx,-idx]
    }

  # return list of individuals to keep
  message('\n')
  iid = iid[!(iid %in% kinshipExclude)]
  return(iid)
}

# apply function and get list of unrelated individuals
message(sprintf(' - loading kinship file %s', kinshipFile))
kinshipTable = read.table(kinshipFile, head=TRUE)

  # remove all imaging participants and their relatives
  imaging_iid = bd$f.eid[!is.na(bd$f.53.2.0) | !is.na(bd$f.53.3.0)]
  idx = imaging_iid %in% kinshipTable$ID1 | imaging_iid %in% kinshipTable$ID2
  exclude = unique(c(imaging_iid,kinshipTable$ID1[idx],kinshipTable$ID2[idx]))
  bd = bd[!(bd$f.eid %in% exclude),]

  # remove remaining relatedness
  unrelated = select.unrelated(kinshipTable,bd$f.eid,86609)
  bd = bd[bd$f.eid %in% unrelated,]

# Sanity check: Test for two subjects in a row in relatedness table
if (sum(kinshipTable$ID1 %in% bd$f.eid & kinshipTable$ID2 %in% bd$f.eid) == 0) {
  message(sprintf(' - relatedness successfully removed\n - %d cases remaining.', dim(bd)[1]))
} else {
  message(' - relatedness has not been removed, exiting script.')
  quit()
}

# =============================
# === add pan ancestry data ===
# =============================
message('\nAdding pan ancestry data.')

# load dataset
message(sprintf(' - loading %s',panFile))
pan = read.delim(panFile, sep = '\t', header = TRUE)
message(sprintf(' - loading %s',bridgeFile))
bridge = read.delim(bridgeFile, sep = ' ', header = FALSE)

# prepare pan-ancestry data
# remove duplicates with related flag 'true'
# add pan to dataset
dupids = pan$s[duplicated(pan$s)]
pan = pan[!(pan$s %in% dupids & pan$related == 'true'),]
names(bridge) = c('f.eid', 's')
pan = inner_join(pan, bridge, 's')
bd = left_join(bd, pan, 'f.eid')

# ==============================================================
# === Create variable file and ancestry-stratified iid files ===
# ==============================================================
message('\nCreating variable file and ancestry-stratified iid files.')

# get caucasian ancestry and pan-ancestry population data
message(' - preparing output variables.')
caucasian = as.numeric(!is.na(bd$f.22006))
pan = bd$pop

# discovery or replication?
bd$discovery = NA
bd$discovery[caucasian==1] = 1
bd$discovery[caucasian==0 & !is.na(pan)] = 0
bd = bd[!is.na(bd$discovery),]
message(sprintf(' - discovery n = %d | replication n = %d',sum(bd$discovery==1),sum(bd$discovery==0)))

# calculate exact age
birth_year_month = paste(as.numeric(bd$f.34.0.0),as.numeric(bd$f.52.0.0),"15", sep="-")
age = as.numeric((as.Date(bd$f.53.0.0,"%Y-%m-%d") - as.Date(birth_year_month,"%Y-%m-%d"))/365.24219052)

# prepare assessment center
ac.f = factor(bd$f.54.0.0)
ac.dummies = model.matrix(~ac.f)
ac = matrix(NA, nrow = length(ac.f), ncol = ncol(ac.dummies)-1)
ac[!is.na(ac.f),] = ac.dummies[,c(2:ncol(ac.dummies))]
ac = data.frame(ac)
names(ac) = paste0('ac',1:ncol(ac))
ac.dummy = ac

# genotyping array dichotom
array = bd$f.22000.0.0
array[array>0] = 1
array[array<0] = 0

# Covs data.frame
df = data.frame(
          FID = bd$f.eid,
          IID = bd$f.eid,
          wba = bd$discovery,
          pan = bd$pop,
          sex = bd$f.31.0.0,
          age = age,
          age2 = age^2,
          ageXsex = age*bd$f.31.0.0,
          age2Xsex = age^2*bd$f.31.0.0,
          headScale = bd$f.25000.2.0,
          x = bd$f.25756.2.0,
          y = bd$f.25757.2.0, 
          z = bd$f.25758.2.0,
          rfMRI = bd$f.25741.2.0, 
          tfMRI = bd$f.25742.2.0,
          array = array,
          ac.dummy) %>%
      cbind(data.frame(
          PC1 = bd$f.22009.0.1, PC2 = bd$f.22009.0.2, PC3 = bd$f.22009.0.3, PC4 = bd$f.22009.0.4, PC5 = bd$f.22009.0.5,
          PC6 = bd$f.22009.0.6, PC7 = bd$f.22009.0.7, PC8 = bd$f.22009.0.8, PC9 = bd$f.22009.0.9, PC10 = bd$f.22009.0.10,
          PC11 = bd$f.22009.0.11, PC12 = bd$f.22009.0.12, PC13 = bd$f.22009.0.13, PC14 = bd$f.22009.0.14, PC15 = bd$f.22009.0.15,
          PC16 = bd$f.22009.0.16, PC17 = bd$f.22009.0.17, PC18 = bd$f.22009.0.18, PC19 = bd$f.22009.0.19, PC20 = bd$f.22009.0.20,
          PC21 = bd$f.22009.0.21, PC22 = bd$f.22009.0.22, PC23 = bd$f.22009.0.23, PC24 = bd$f.22009.0.24, PC25 = bd$f.22009.0.25,
          PC26 = bd$f.22009.0.26, PC27 = bd$f.22009.0.27, PC28 = bd$f.22009.0.28, PC29 = bd$f.22009.0.29, PC30 = bd$f.22009.0.30,
          PC31 = bd$f.22009.0.31, PC32 = bd$f.22009.0.32, PC33 = bd$f.22009.0.33, PC34 = bd$f.22009.0.34, PC35 = bd$f.22009.0.35,
          PC36 = bd$f.22009.0.36, PC37 = bd$f.22009.0.37, PC38 = bd$f.22009.0.38, PC39 = bd$f.22009.0.39, PC40 = bd$f.22009.0.40,
          PanC1 = bd$PC1, PanC2 = bd$PC2, PanC3 = bd$PC3, PanC4 = bd$PC4, PanC5 = bd$PC5, 
          PanC6 = bd$PC6, PanC7 = bd$PC7, PanC8 = bd$PC8, PanC9 = bd$PC9, PanC10 = bd$PC10,
          PanC11 = bd$PC11, PanC12 = bd$PC12, PanC13 = bd$PC13, PanC14 = bd$PC14, PanC15 = bd$PC15,
          PanC16 = bd$PC16, PanC17 = bd$PC17, PanC18 = bd$PC18, PanC19 = bd$PC19, PanC20 = bd$PC20))

# write subject IIDs and vars
message(' - writing output files.')
system(sprintf('mkdir -p %s',targetDir))
data.table::fwrite(df, file = sprintf('%s/r2025.vars.gz',targetDir), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, compress = 'gzip')
data.table::fwrite(df[df$wba == 1,] , file = sprintf('%s/r2025.vars.discovery.gz',targetDir), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, compress = 'gzip')
data.table::fwrite(df[df$wba == 0,] , file = sprintf('%s/r2025.vars.replication.gz',targetDir), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE, compress = 'gzip')
write.table(data.frame(FID = df$FID[df$wba==1], IID = df$IID[df$wba==1]), file = sprintf('%s/iid.discovery.txt',targetDir), quote = F, sep = '\t', row.names = F)
write.table(data.frame(FID = df$FID[df$wba==0], IID = df$IID[df$wba==0]), file = sprintf('%s/iid.replication.txt',targetDir), quote = F, sep = '\t', row.names = F)

# write subject IIDs for each ancestry (replication sample only)
repl = df[df$wba==0,]
for (anc in names(table(repl$pan))) {
  tmp = data.frame(FID = repl$FID[repl$pan==anc & !is.na(repl$pan)], IID = repl$IID[repl$pan==anc & !is.na(repl$pan)])
  write.table(tmp, file = sprintf('%s/iid.replication.%s.txt',targetDir,anc), sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}
system(sprintf('chmod 750 %s/iid*',targetDir))
system(sprintf('chmod 750 %s/r2025*',targetDir))

# write phesant file for all individuals
setnames(bd,old = names(bd),new = gsub("^f_", "x", gsub("\\.", "_", names(bd))))
names(bd)[1] = "eid"
bd = bd[,1:max(grep("^x[0-9]+_[0-9]+_[0-9]+$", names(bd)))]
data.table::fwrite(bd, file = sprintf("%s/ukb678162.phesant.2025.csv",targetDir), sep = ",", quote = TRUE, na = "", row.names = FALSE, col.names = TRUE, compress = 'none')

  # write phesant file for white-British (discovery sample)
  wba = bd[df$wba==1,]
  data.table::fwrite(wba, file = sprintf("%s/ukb678162.phesant.2025.wba.csv",targetDir), sep = ",", quote = TRUE, na = "", row.names = FALSE, col.names = TRUE, compress = 'none')

# merge instances
message(" - merging multiple instances.")
varInfo = read.delim(phesantVarFile, header = T, sep = '\t')
fieldId = varInfo$FieldID[varInfo$CAT_SINGLE_TO_CAT_MULT!='YES-INSTANCES']

# Example:
# bd has columns like bd$x20195_0_0, bd$x20195_1_0, bd$x20195_2_0, etc.
# fieldId is a vector like c(20195, 30000, ...)
merge_instances <- function(bd, fieldId) {
  total_fields <- length(fieldId)
  pb <- txtProgressBar(min = 0, max = total_fields, style = 3)
  
  counter <- 0
  
  for (fid in fieldId) {
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
    
    # Get all columns for this field
    cols <- grep(paste0("^x", fid, "_[0-9]+_[0-9]+$"), names(bd), value = TRUE)
    if (length(cols) == 0) next
    
    # Extract all instances and arrays
    parsed <- data.frame(
      col = cols,
      instance = as.numeric(sub(paste0("^x", fid, "_([0-9]+)_.*$"), "\\1", cols)),
      array = sub(".*_([0-9]+)$", "\\1", cols),
      stringsAsFactors = FALSE
    )
    
    instances <- sort(unique(parsed$instance))
    arrays <- unique(parsed$array)
    
    # Determine target instance (prefer 0 if exists, otherwise lowest)
    target_instance <- if (0 %in% instances) 0 else min(instances)
    
    for (array in arrays) {
      subset_cols <- parsed[parsed$array == array, ]
      subset_cols <- subset_cols[order(subset_cols$instance), ]  # sort by instance
      
      target_col <- paste0("x", fid, "_", target_instance, "_", array)
      
      # If target column doesn't exist, create it with the same type as the first instance column
      if (!(target_col %in% names(bd))) {
        first_col <- subset_cols$col[1]
        bd[[target_col]] <- as(NA, class(bd[[first_col]]))
      }
      
      # Vectorized merging
      for (col in subset_cols$col) {
        if (col == target_col) next
        
        source_vals <- bd[[col]]
        target_vals <- bd[[target_col]]
        
        # Fill only where target is missing AND source is not missing
        update_idx <- (is.na(target_vals) | target_vals == "") & 
                      !(is.na(source_vals) | source_vals == "")
        
        # Ensure no NA values in update_idx
        update_idx[is.na(update_idx)] <- FALSE
        
        # Update only the missing positions
        if (any(update_idx)) {
          bd[[target_col]][update_idx] <- source_vals[update_idx]
        }
      }
    }
  }
  
  close(pb)
  return(bd)
}

# apply function and save phesant file for all individuals
bd_merged = merge_instances(bd, fieldId)
data.table::fwrite(bd_merged, file = sprintf("%s/ukb678162.phesant.2025.mergedInstances.csv",targetDir), sep = ",", quote = TRUE, na = "", row.names = FALSE, col.names = TRUE, compress = 'none')

  # write phesant file for white-British (discovery sample)
  wba_merged = bd_merged[df$wba==1,]
  data.table::fwrite(wba_merged, file = sprintf("%s/ukb678162.phesant.2025.wba.mergedInstances.csv",targetDir), sep = ",", quote = TRUE, na = "", row.names = FALSE, col.names = TRUE, compress = 'none')

message('\n--- Completed: Sample Filtering ---')

