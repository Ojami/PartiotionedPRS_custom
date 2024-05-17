prep_z_matrix <- function(z_mat, N_mat,
                          corr_cutoff=0.85,
                          rm_traits=NULL,
                          pval_cutoff=NULL) {
  
  # Given a matrix of z-scores (N_variants x M_traits) and vector of median
  # sample sizes per trait:
  # 1) perform final pre-processing steps before bNMF clustering:
  # trait filtering by p-value, trait pruning based on correlation,
  # and z-score scaling based on sample size
  # 2) expand N x M matrix into N x 2M non-negative matrix
  
  
  # if no user input for pval_cutoff, use bonferroni definition
  if (is.null(pval_cutoff)){
    pval_cutoff <- 0.05 / nrow(z_mat) 
  }
  print(paste0(sum(is.na(z_mat)), " missing values before pval pruning."))
  
  df_traits <- data.frame(trait = colnames(z_mat))
  
  df_traits_filtered <- data.frame(trait=as.character(),
                                   result=as.character(),
                                   note=as.character())
  
  # remove manually entered traits
  if (!is.null(rm_traits)) {
    print(sprintf("Removing %i manually entered traits!!!", length(rm_traits)))
    z_mat <- z_mat[, !colnames(z_mat) %in% rm_traits]
    
    df_remove_manual <- data.frame(trait=rm_traits,
                                   result="removed (manual)",
                                   note=NA)
    df_traits_filtered <- rbind(df_traits_filtered, df_remove_manual)
  }
  
  # Filter traits by p-value (min. p-value < 0.05/N_variants)
  print("Filtering traits w/ no pvalues below cutoff...")
  minP_vec <- apply(z_mat, 2, function(x) min(2 * pnorm(abs(x), lower.tail=F), na.rm=T))
  traits_removed <- colnames(z_mat)[minP_vec >= pval_cutoff]
  df_lowP_vec <- data.frame(minPval=minP_vec[minP_vec >= pval_cutoff]) %>%
    rownames_to_column('trait') %>%
    mutate(result="removed (p-value)") %>% 
    mutate(note=paste("min pval=",format(minPval,scientific=T))) %>%
    dplyr::select(trait, result, note)
  
  df_traits_filtered <- rbind(df_traits_filtered, df_lowP_vec)
  cat(paste(sprintf("Removing traits with no variant having p < %.3e",pval_cutoff),
            paste(traits_removed,
                  collapse="\n"),sep=":\n"))
  z_mat <- z_mat[, minP_vec < pval_cutoff]
  
  # Prune traits by correlation (remove traits with Pearson |r| > 0.85)
  
  cat(sprintf("\n\nPrune traits by correlation (remove traits with Pearson |r| > %.2f)\n",
              corr_cutoff))
  trait_cor_mat <- cor(z_mat, use="pairwise.complete.obs")  # Trait-trait correlation matrix
  write.table(trait_cor_mat,"./trait_cor_mat.txt", sep="\t")
  
  # sort by max(z) instead of min(pval)
  remaining_traits <- names(sort(apply(z_mat, 2, max, na.rm=T),decreasing = T))
  print(paste("Initial number of traits:",length(remaining_traits)))
  
  keep_traits <- c()
  df_corr_removed <- c()
  
  while (length(remaining_traits) > 0) {
    # append to list of traits to keep
    keep_traits <- c(keep_traits, remaining_traits[1])
    if (length(remaining_traits)==1) {
      break
    }
    trait_cor_mat <- trait_cor_mat[remaining_traits, remaining_traits]
    print(dim(trait_cor_mat))
    to_remove <- rownames(trait_cor_mat)[abs(trait_cor_mat[, remaining_traits[1]]) >= corr_cutoff]
    
    # track results in dataframe
    if (length(to_remove)>1) {
      df_corr_removed_tmp <- data.frame(trait=to_remove[to_remove!=remaining_traits[1]],
                                        result="removed (correlation)",
                                        note=paste("correlated w/", remaining_traits[1]))
      df_traits_filtered <- rbind(df_traits_filtered, df_corr_removed_tmp)
    }
    
    if (length(to_remove) > 1) {
      cat(paste(sprintf("Correlated trait being removed for %s:",remaining_traits[1]),
                paste(to_remove[!to_remove %in% remaining_traits[1]],collapse="\n"),"\n",sep="\n"))
    }
    
    # remove all correlated traits from remaining list
    remaining_traits <- setdiff(
      remaining_traits, 
      to_remove
    )
  }
  
  df_traits_filtered <- df_traits_filtered %>%
    right_join(df_traits, by="trait") %>%
    mutate(result = ifelse(is.na(result), "trait kept", result))
  
  pruned_traits <- df_traits_filtered %>%
    filter(result!="trait kept") %>%
    pull(trait)
  cat(paste("Traits removed in pruning process:", 
            paste(pruned_traits, collapse="\n"),sep = "\n"))
  
  cat(paste0("\nNumber of remaining traits: ",length(keep_traits),"\n"))
  
  # cat(paste("Remaing traits:", 
  #           paste(keep_traits, collapse="\n"),sep = "\n")) 
  print(paste0(sum(is.na(z_mat)), " missing values before corr pruning."))
  z_mat <- z_mat[, keep_traits]
  
  # Adjust z-scores by sample size for each variant-trait combo
  # i.e. (z = z / sqrt(medN) * mean(sqrt(medN_all_traits)))
  cat("\n\n")
  print("Performing sample size adjustment...")
  medN_vec <- apply(N_mat[, colnames(z_mat)], 2, median, na.rm=T)
  z_mat <- z_mat / sqrt(N_mat[, colnames(z_mat)]) * mean(sqrt(medN_vec))
  
  # Replace missing values with zero
  print("Replacing remaining missing values with zero...")
  print(paste0(sum(is.na(z_mat)), " missing values were replaced."))
  z_mat[is.na(z_mat)] <- 0
  
  # Expand into N x 2M non-negative matrix
  print("Expanding z-score matrix into non-negative matrix (N-variants x 2M-traits)...")
  z_mat_pos <- z_mat
  z_mat_pos[z_mat_pos < 0] <- 0
  colnames(z_mat_pos) <- paste0(colnames(z_mat), "_pos")
  z_mat_neg <- -z_mat
  z_mat_neg[z_mat_neg < 0] <- 0
  colnames(z_mat_neg) <- paste0(colnames(z_mat), "_neg")
  final_z_mat <- cbind(z_mat_pos, z_mat_neg)
  
  output <- list(final_z_mat=final_z_mat,
                 df_traits=df_traits_filtered)
}