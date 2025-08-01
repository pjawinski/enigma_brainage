#!/usr/bin/env Rscript

# ==========================================================
# === create phesant summary text file and combined plot ===
# ==========================================================

# get arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=4) {
  stop(paste0('expected 4 arguments, but ', length(args), ' argument(s) provided.'), call.=FALSE)
}

# set arguments
traits=args[1] # traits="rc_auc,rc_phi,rc_range"
phesantSummary=args[2] # phesantSummary="results/rc_auc/phesant/phewas_summary.txt,results/rc_phi/phesant/phewas_summary.txt,results/rc_range/phesant/phewas_summary.txt"
phesantPlot=args[3] # phesantPlot="results/rc_auc/phesant/phewas.png,results/rc_phi/phesant/phewas.png,results/rc_range/phesant/phewas.png"
outputFile=args[4] # outputFile="results/combined/suppl.phewas"

message(paste0('\n--- Create phesant summary text file and combined plot ---',
               '\ntraits: ', traits,
               '\nphesantSummary: ', phesantSummary,
               '\nphesantPlot: ', phesantPlot,
               '\noutputFile: ', outputFile,'\n'))

# attach packages to current R session
for (pkg in c('data.table', 'dplyr', 'ggplot2', 'ggpubr', 'plotly', 'patchwork','magick','stringr')) {
  eval(bquote(suppressPackageStartupMessages(require(.(pkg)))))
}

# transform variables
phesantSummary = str_split(phesantSummary, ',')[[1]]
phesantPlot = str_split(phesantPlot, ',')[[1]]
traits = str_split(traits, ',')[[1]]

# open and merge datasets
for (i in 1:length(traits)) {
  
  # open datasets and keep relevant variables
  tmp = read.delim(phesantSummary[i], sep = '\t', header = T)
  if (i ==1) {
    tmp = tmp[,c('description', 'varName', 'custom_category', 'Path', 'varType', 'resType', 'n', 'ntotal', 'beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'rho')]
  } else { 
    tmp = tmp[,c('varName', 'beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'rho')]
  }
  tmp$rhoAbs = abs(tmp$rho)
  
  # create unique variable names
  names(tmp)[which(names(tmp) %in% c('beta', 'lower', 'upper', 'se', 'pvalue', 'fdr', 'z', 'rho', 'rhoAbs'))] = 
    paste0(traits[i],c('_beta', '_lower', '_upper', '_se', '_pvalue', '_fdr', '_z', '_rho', '_rhoAbs'))
  
  # merge datasets
  if (i ==1) { df = tmp } else { df = inner_join(df,tmp, by = 'varName') }
}
 
# get ranking by rho
rho_cols <- grep("_rhoAbs$", names(df), value = TRUE)
df$trait_ranking <- apply(df[, rho_cols], 1, function(x) {
  traits_sorted <- rho_cols[order(-x)] |> sub("_rhoAbs$", "", x = _)
  paste(traits_sorted, collapse = " > ")
})

# get top p-value and top |rho|
df$top_pvalue =  df[,grep('pvalue',names(df))] %>% apply(1, FUN = min)
df$top_fdr = df[,grep('fdr',names(df))] %>% apply(1, FUN = min)
df$top_rhoAbs = df[,grep('rhoAbs',names(df))] %>% apply(1, FUN = max)

# create output
message('Writing txt file.')
df$`NA` = ""
cols = c('varName', 'description', 'ntotal', 'top_rhoAbs','top_pvalue','top_fdr', 'trait_ranking')
for (i in 1:length(traits)) {
  cols = c(cols,'NA',paste0(traits[i],c('_beta', '_se', '_rho', '_pvalue', '_fdr')))
}
cols = c(cols, c('NA', 'n', 'varType', 'resType', 'custom_category', 'Path'))
output = df[,cols]
output = output[order(output$top_fdr,output$top_pvalue),]
write.table(output, paste0(outputFile,'.txt'), sep = '\t', quote = F, row.names = F)

# combine phesant plots
if (length(phesantPlot) > 0 && phesantPlot != "") {
  for (i in 1:length(traits)) {
    tmp = image_read(phesantPlot[i])
    tmp = ggplot() + background_image(tmp) + coord_fixed(ratio = image_info(tmp)$height/image_info(tmp)$width)
    if (i ==1) { plot = tmp } else { plot = plot / tmp }
  }

  # draw plot with annotations
  plot = plot + plot_annotation(tag_levels ='a') & 
    theme(plot.tag.position = c(0, 1), plot.tag = element_text(size = 22))

  # save plot
  message('Writing png file.')
  png(width = 8.7, height = 5.0*length(traits), units = "in", res = 600, filename = paste0(outputFile,'.png'))
  plot
  invisible(dev.off())
}

# get summary statistics
message('Calculating summary statistics.')

  # Get the number of phenotypes for Bonferroni correction
  n_phenotypes = nrow(df)
  bonferroni_threshold = 0.05 / n_phenotypes
  message(paste0("Bonferroni threshold: ", bonferroni_threshold))

  # Helper function to format summary tables
  table_to_text <- function(title, df_table) {
    lines <- c(title,
               capture.output(print(df_table, row.names = FALSE)),
               "") # Empty line after each table
    return(lines)
  }

  # Helper function to compute average rank
  get_avg_rank <- function(df, trait) {
    ranks <- sapply(df$trait_ranking, function(ranking) {
      traits_ordered <- unlist(strsplit(ranking, " > "))
      match(trait, traits_ordered)  # returns position of trait
    })
    mean(ranks, na.rm = TRUE)
  }

  # overall summary statistics
  overall_summary = data.frame(Measure = c("Total", "meanRank", "Nominal", "FDR", "Bonferroni"))
  for (trait in traits) {
    total_assoc = sum(!is.na(df[[paste0(trait, "_pvalue")]]))
    nominal_sig = sum(df[[paste0(trait, "_pvalue")]] < 0.05, na.rm = TRUE)
    fdr_sig = sum(df[[paste0(trait, "_fdr")]] < 0.05, na.rm = TRUE)
    bonf_sig = sum(df[[paste0(trait, "_pvalue")]] < bonferroni_threshold, na.rm = TRUE)
    avg_rank = get_avg_rank(df, trait)
    overall_summary[[trait]] = c(total_assoc, sprintf("%.2f", avg_rank), nominal_sig, fdr_sig, bonf_sig)
  }

  # Print overall summary table
  log_text = c() 
  log_text = c(log_text, table_to_text("\nOverall summary", overall_summary))

  # category-wise summary statistics
  categories = unique(df$custom_category)
  summary_stats = data.frame()
  for (cat in categories) {
    df_cat = df[df$custom_category == cat, ]
    cat_table = data.frame(Measure = c("Total", "meanRank", "Nominal", "FDR", "Bonferroni"))

    for (trait in traits) {
      total_assoc = sum(!is.na(df_cat[[paste0(trait, "_pvalue")]]))
      nominal_sig = sum(df_cat[[paste0(trait, "_pvalue")]] < 0.05, na.rm = TRUE)
      fdr_sig = sum(df_cat[[paste0(trait, "_fdr")]] < 0.05, na.rm = TRUE)
      bonf_sig = sum(df_cat[[paste0(trait, "_pvalue")]] < bonferroni_threshold, na.rm = TRUE)
      avg_rank = get_avg_rank(df_cat, trait)
      cat_table[[trait]] = c(total_assoc, sprintf("%.2f", avg_rank), nominal_sig, fdr_sig, bonf_sig)
    }

    # Print category label and table
    log_text = c(log_text, table_to_text(paste0("Category: ", cat), cat_table))
  }

  # Convert all log lines into one single string
  final_log = paste(log_text, collapse = "\n")
  message(final_log)
  writeLines(final_log, paste0(outputFile, ".stats.txt"))

system(paste0('chmod 770 ', outputFile, '*'))
message('-- Completed: Creating suppl. table of PheWAS results ---')
