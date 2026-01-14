#BiocManager::install("sangerseqR")
library(sangerseqR)
library(Biostrings)
library(dplyr)
library(tidyr)
library(ggplot2)

# fuction: calculate the methylation level and merge data 
run_all_direct_BSP <- function(samples, C_rev, BASE_PATH, output_prefix) {
  merged_df <- NULL
  #step1: read ab1 and calculate methylation level
  for(sample_name in samples_order){
    ab1_file <- paste0(sample_name, "_DMR1_DMR1_R.ab1")
    
    ab1 <- read.abif(file.path(BASE_PATH, ab1_file))
    sseq <- sangerseq(ab1)
    trace <- sseq@traceMatrix
    A_signal <- trace[,1]  # A
    G_signal <- trace[,3]  # G

    df <- C_rev
    meth_col <- paste0(sample_name, "_meth_ratio")
    df[[meth_col]] <- NA
    
    # calculate methylation ratio
    for(j in seq_len(nrow(df))){
      pos <- df$reverse_pos[j]
      if(pos <= length(A_signal)){
        A <- A_signal[pos]
        G <- G_signal[pos]
        if((A + G) > 0){
          df[[meth_col]][j] <- G / (G + A)}}}
    
    # save each sample
    out_file <- paste0(sample_name, "_", output_prefix, "_meth.tsv")
    write.table(df, file.path(BASE_PATH, out_file), sep="\t", row.names=FALSE, quote=FALSE)
    message("Saved ", out_file)
    }
  #step2: merge all the sample
  for(sample_name in samples_order){
    file <- file.path(BASE_PATH, paste0(sample_name, "_", output_prefix, "_meth.tsv"))
    df <- read.table(file, header=TRUE, sep="\t", stringsAsFactors = FALSE)
    
    meth_col <- grep("_meth_ratio$", names(df), value = TRUE)
    df <- df %>% select(genomic_pos, context, all_of(meth_col))
    if(is.null(merged_df)){
      merged_df <- df
    } else {
      # full join 
      merged_df <- merge(merged_df, df, by=c("genomic_pos","context"), all=TRUE)}}
  
  # step3:arrange by genomic_pos 
  merged_df <- merged_df[order(merged_df$genomic_pos), ]
  
  # step4:save
  out_file <- file.path(BASE_PATH, paste0(output_prefix, "_all_samples_meth_ratio.tsv"))
  write.table(merged_df,out_file,sep = "\t",row.names = FALSE,quote = FALSE)
  message("Saved merged data to ", out_file)
  return(merged_df)
}

#fuction2: plot
run_BSP_analysis <- function(merged_df, BASE_PATH, context,Tr_cols, Bi_cols , Ti_cols,
                             sample_order = NULL, output_prefix) {
  merged_df <- merged_df %>% mutate(context = as.character(context))
  
  # filter by context
  plot_df <- merged_df %>%
    filter(context == !!context) %>%
    pivot_longer(
      cols = ends_with("_meth_ratio"),
      names_to = "sample",
      values_to = "meth_ratio") %>%
    mutate(sample = sub("_meth_ratio$", "", sample),
           group = case_when(grepl("^Tr", sample) ~ "Tr",
                             grepl("^Bi", sample) ~ "Bi",
                             grepl("^Ti", sample) ~ "Ti"),
           meth_bin = cut(meth_ratio * 100, breaks = c(0, 20, 40, 60, 80, 100),labels = c("0-20", "20-40", "40-60", "60-80", "80-100"),include.lowest = TRUE))
  
  if(!is.null(sample_order)){plot_df$sample <- factor(plot_df$sample, levels = sample_order)}
  
  #lolipop_plot
  p_position <- ggplot(plot_df, aes(x = genomic_pos, y = sample)) +
    geom_line(aes(group = sample), color = "grey50", linewidth = 0.8) +
    geom_point(aes(fill = meth_bin), shape = 21, size = 4, color = "black",
               position = position_jitter(width = 0.5, height = 0)) +
    scale_fill_manual(
      values = c("0-20"= "white","20-40"="grey85","40-60" = "grey65",
                 "60-80"= "grey35","80-100" = "black"),
      name = "Methylation level (%)") +
    theme_classic(base_size = 14) +
    xlab(paste0("Genomic position (", context, ")")) +
    ylab(NULL)
  
  ggsave(file.path(BASE_PATH, paste0(output_prefix, "_", context, "_lolipop_plot.png")),
         p_position, width = 10, height = 6)
  
  context_df <- merged_df %>% filter(.data$context == !!context)
  Tr_avg <- context_df %>% summarise(across(all_of(Tr_cols), mean, na.rm = TRUE))
  Bi_avg <- context_df %>% summarise(across(all_of(Bi_cols), mean, na.rm = TRUE))
  Ti_avg <- context_df %>% summarise(across(all_of(Ti_cols), mean, na.rm = TRUE))
  
  Tr_avg_pct <- as.numeric(Tr_avg) * 100
  Bi_avg_pct <- as.numeric(Bi_avg) * 100
  Ti_avg_pct <- as.numeric(Ti_avg) * 100
  
  # calculate the average methlation & delata methylation & t-test
  Tr_minus_Bi <- mean(Tr_avg_pct) - mean(Bi_avg_pct)
  Tr_minus_Ti <- mean(Tr_avg_pct) - mean(Ti_avg_pct)
  
  Tr_vs_Bi_pvalue <- t.test(Tr_avg, Bi_avg, var.equal = FALSE)$p.value
  Tr_vs_Ti_pvalue <- t.test(Tr_avg, Ti_avg, var.equal = FALSE)$p.value
  
  # summary table
  Tr_df <- setNames(as.list(Tr_avg_pct), paste0("Tr", seq_along(Tr_avg_pct), "_avg"))
  Bi_df <- setNames(as.list(Bi_avg_pct), paste0("Bi", seq_along(Bi_avg_pct), "_avg"))
  Ti_df <- setNames(as.list(Ti_avg_pct), paste0("Ti", seq_along(Ti_avg_pct), "_avg"))
  
  # summary row
  summary_df <- data.frame(
    Tr_mean = mean(Tr_avg_pct),Bi_mean = mean(Bi_avg_pct),Ti_mean = mean(Ti_avg_pct),
    Tr_minus_Bi = Tr_minus_Bi,Tr_minus_Ti = Tr_minus_Ti,
    Tr_vs_Bi_pvalue = Tr_vs_Bi_pvalue,Tr_vs_Ti_pvalue = Tr_vs_Ti_pvalue)
  
  # merge these result
  result <- cbind(Tr_df, Bi_df, Ti_df, summary_df)
  
  out_summary_file <- file.path(BASE_PATH, paste0(output_prefix, "_", context, "_methylation_summary.tsv"))
  write.table(result, out_summary_file, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Saved summary to ", out_summary_file)
  
  # boxplot
  df_box <- data.frame(
    sample = c(names(Tr_avg), names(Bi_avg), names(Ti_avg)),
    group = factor(rep(c("Tr", "Bi", "Ti"),times = c(length(Tr_avg), length(Bi_avg), length(Ti_avg))),levels = c("Tr", "Bi", "Ti")),
    meth_ratio = c(Tr_avg_pct, Bi_avg_pct, Ti_avg_pct))
  if(!is.null(sample_order)) df_box$sample <- factor(df_box$sample, levels = sample_order)
  group_colors <- c("Tr"="#3F995A", "Bi"="#6295C3", "Ti"="#636196")
  p_box <- ggplot(df_box, aes(x = group, y = meth_ratio, color = group)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.7, color = "black") +
    geom_jitter(width = 0.1, size = 2, color = "black") +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 15) +
    ylab("Average methylation level (%)") +  
    xlab("Group") +
    theme(
      axis.line = element_line(size = 1.1),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(face = "bold"))
  
  ggsave(file.path(BASE_PATH, paste0(output_prefix, "_", context, "_boxplot.png")),
         p_box, width = 6, height = 5)
  
  # barplot
  p_bar <- ggplot(df_box, aes(x = group, y = meth_ratio)) +
    stat_summary(fun = mean, geom = "bar", aes(fill = group), alpha = 0.7, color = "black") +
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
    geom_jitter(width = 0.1, size = 2, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 15) +
    ylab("Average methylation level (%)") +
    xlab("Group") +
    theme(
      axis.line = element_line(size = 1.1),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(face = "bold"))
  
  ggsave(file.path(BASE_PATH, paste0(output_prefix, "_", context, "_barplot_mean.png")),
         p_bar, width = 6, height = 5)
  
  message("Plots saved to ", BASE_PATH)
  
  return(list(summary = result, position_plot = p_position,
              boxplot = p_box, barplot = p_bar))
}
