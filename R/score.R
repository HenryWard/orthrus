######
# SCORING CODE
######

# Inner function to scale values between 0 and 1
scale_values <- function(x) {
  val <- (x-min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
}

#' Scores conditions against a single control.
#' 
#' Scores guides for any number of condition screens against a control screen
#' (e.g. for directly comparing drug response to DMSO response). Do NOT use 
#' this function if you would like to compare guides against a null model 
#' derived from single-gene effects. In that case, call \code{score_combn_vs_single}
#' instead. After running this function, pass the resulting dataframe to 
#' \code{call_condition_hits} to call significant effects.
#' 
#' @param guides A list of guides returned from \code{split_guides}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param control_screen_name Name of a control screen to test condition screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param separate_orientation If true, then guide values are scored separately across each 
#'   orientation (default FALSE).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").
#' @param filter_genes List of genes to filter from scoring (default NULL).
#' @param return_residuals If FALSE, returns NA instead of residuals dataframe (default TRUE).
#'   This is recommend if scoring large datasets and memory is a limitation.  
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A list containing two dataframes. The first entry, named "scored_data" in the list,
#'   contains scored data with separate columns given by the specified control and condition
#'   names. The second entry, named "residuals" in the list, is a dataframe containing control,
#'   condition and loess-normalized residuals for all guides.
#' @export
score_conditions_vs_control <- function(guides, screens, control_screen_name, condition_screen_names, 
                                        separate_orientation = FALSE, min_guides = 3, test = "moderated-t", 
                                        loess = TRUE, fdr_method = "BY", filter_genes = NULL,
                                        return_residuals = TRUE, verbose = FALSE) {
  
  # Filters specified genes from dataset before scoring
  if (!is.null(filter_genes)) {
    for (gene in filter_genes) {
      ind <- unlist(lapply(guides, function(x) x[["gene1"]] == gene | x[["gene2"]] == gene))
      if (sum(ind) > 0) {
        cat(paste("Removing", sum(ind), "gene pairs containing gene", gene, "\n"))
        guides <- guides[!ind]
      }
    }
  }
  
  # Runs separately on each orientation if specified
  results1 <- NULL
  results2 <- NULL
  if (!separate_orientation) {
    results <- score_conditions_vs_control_inner(guides, screens, control_screen_name, condition_screen_names, 
                                                 min_guides = min_guides, test = test, loess = loess,
                                                 return_residuals = return_residuals, verbose = verbose)
    return(results)
  } else {
    results1 <- score_conditions_vs_control_inner(guides, screens, control_screen_name, condition_screen_names, 
                                                  min_guides = min_guides, test = test, loess = loess,
                                                  return_residuals = return_residuals, verbose = verbose,
                                                  screen_prefix = "orient1_")
    results2 <- score_conditions_vs_control_inner(guides, screens, control_screen_name, condition_screen_names, 
                                                  min_guides = min_guides, test = test, loess = loess,
                                                  return_residuals = return_residuals, verbose = verbose,
                                                  screen_prefix = "orient2_")
    return(list(results1, results2))
  }
}

# Inner function for the above
score_conditions_vs_control_inner <- function(guides, screens, control_screen_name, condition_screen_names, 
                                              min_guides = 3, test = "moderated-t", 
                                              loess = TRUE, screen_prefix = "", fdr_method = "BY",
                                              return_residuals = TRUE, verbose = FALSE) {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_screen_name
  control_cols <- paste0(screen_prefix, screens[[control_name]][["replicates"]])
  condition_names <- c()
  condition_cols <- list()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
    condition_cols[[condition]] <- paste0(screen_prefix, screens[[condition]][["replicates"]])
  }
  
  # Makes output dataframe
  scores <- data.frame(gene1 = rep(NA, length(guides)), gene2 = NA)
  
  # Makes residual dataframes if necessary
  max_guides <- -1
  condition_residuals <- list()
  if (test == "moderated-t") {
    
    # Gets max number of guides first
    for (guide in guides) {
      guide_names <- names(guide)[!(names(guide) %in% c("gene1", "gene2", "guide_type"))]
      guide_num <- max(unlist(lapply(guide_names, function(x) length(guide[[x]]))))
      max_guides <- max(max_guides, guide_num)
    }
    
    # Makes residual dataframes with columns equal to the max number of guides
    for (name in condition_names) {
      residual_df <- data.frame(matrix(nrow = length(guides), ncol = max_guides))
      colnames(residual_df) <- paste0("guide_residual_", 1:max_guides)
      condition_residuals[[name]] <- residual_df
    }
  }
  
  # Makes loess residual dataframe if specified
  loess_residuals <- NULL
  if (return_residuals & test == "moderated-t") {
    loess_residuals <- data.frame(n = rep(0, max_guides*length(guides)))
    loess_residuals[[paste0("mean_", control_name)]] <- rep(0, nrow(loess_residuals))
    for (name in condition_names) {
      loess_residuals[[paste0("mean_", name)]] <- rep(0, nrow(loess_residuals))
    }
  }
  
  # Appends additional columns for each condition
  new_cols <- c(paste0("n_", control_name), 
                paste0("mean_", control_name),
                paste0("variance_", control_name))
  for (name in condition_names) {
    new_cols <- c(new_cols, c(
      paste0("n_", name), 
      paste0("mean_", name),
      paste0("variance_", name),
      paste0("differential_", name, "_vs_", control_name),
      paste0("pval_", name, "_vs_", control_name),
      paste0("fdr_", name, "_vs_", control_name),
      paste0("significant_", name, "_vs_", control_name)
    ))
  }
  scores[new_cols] <- NA
  
  # Scores guides for each condition
  counter <- 1
  for (i in 1:length(guides)) {
    
    # Gets gene names and control guide values across replicates
    guide_vals <- guides[[i]]
    scores$gene1[i] <- guide_vals$gene1
    scores$gene2[i] <- guide_vals$gene2
    rep_mean_control <- rowMeans(data.frame(guide_vals[control_cols]))
    control_keep_ind <- !is.nan(rep_mean_control) & !is.na(rep_mean_control) & !is.null(rep_mean_control)
    
    # Skips if too few guides
    if (sum(control_keep_ind) < min_guides) {
      next
    }
    
    # Makes loess-normalized residual dataframe if necessary
    ind <- counter:(counter + length(rep_mean_control) - 1)
    if (loess) {
      loess_residuals$n[ind] <- i
      loess_residuals[[paste0("mean_", control_name)]][ind] <- rep_mean_control
    }
    
    # Takes the mean across replicates for all conditions
    for (name in condition_names) {
      rep_mean_condition <- rowMeans(data.frame(guide_vals[condition_cols[[name]]]))
      condition_keep_ind <- !is.nan(rep_mean_condition) & !is.na(rep_mean_condition) & 
        !is.null(rep_mean_condition) & control_keep_ind
      
      # Subsets to same guides and skips if too few guides remaining
      rep_mean_control <- rep_mean_control[condition_keep_ind]
      rep_mean_condition <- rep_mean_condition[condition_keep_ind]
      if (sum(condition_keep_ind) < min_guides) {
        next
      }
      
      # Computes stats
      diff <- rep_mean_condition - rep_mean_control
      scores[[paste0("n_", control_name)]][i] <- length(rep_mean_control)
      scores[[paste0("n_", name)]][i] <- length(rep_mean_condition)
      scores[[paste0("mean_", control_name)]][i] <- mean(rep_mean_control)
      scores[[paste0("mean_", name)]][i] <- mean(rep_mean_condition)
      scores[[paste0("variance_", control_name)]][i] <- stats::var(rep_mean_control)
      scores[[paste0("variance_", name)]][i] <- stats::var(rep_mean_condition)
      scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff)
      
      # Appends mean LFCs for loess-normalization if specified
      if (loess) {
        n_condition <- length(loess_residuals[[paste0("mean_", name)]][ind])
        loess_residuals[[paste0("mean_", name)]][ind] <- 
          c(rep_mean_condition, rep(NA, n_condition - length(rep_mean_condition)))
      }
      
      # Performs the specified type of testing or stores residuals for later testing
      if (test == "rank-sum") {
        scores[[paste0("pval_", name, "_vs_", control_name)]][i] <- 
          suppressWarnings(stats::wilcox.test(rep_mean_condition, rep_mean_control))$p.value
      } else if (test == "moderated-t") {
        if (length(diff) < max_guides) { diff <- c(diff, rep(NA, max_guides - length(diff))) } 
        condition_residuals[[name]][i,1:max_guides] <- diff 
      }
    }
    counter <- counter + length(rep_mean_control)
  }
   
  # Computes loess-normalized residuals if specified
  if (return_residuals & test == "moderated-t") {
    loess_residuals <- loess_residuals[1:counter,]
    control_values <- loess_residuals[[paste0("mean_", control_name)]]
    for (name in condition_names) {
      condition_values <- loess_residuals[[paste0("mean_", name)]]
      temp <- loess_MA(control_values, condition_values)
      loess_residuals[[paste0("loess_residual_", name)]] <- temp[["residual"]]
      loess_residuals[[paste0("loess_predicted_", name)]] <- temp[["predicted"]]
    }
    
    # Replaces residuals with loess-normalized residuals
    for (i in 1:length(guides)) {
      for (name in condition_names) {
        resid <- loess_residuals[[paste0("loess_residual_", name)]][loess_residuals$n == i]
        predicted <- loess_residuals[[paste0("loess_predicted_", name)]][loess_residuals$n == i]
        if (length(resid) < max_guides) { 
          resid <- c(resid, rep(NA, max_guides - length(resid))) 
        } 
        condition_residuals[[name]][i,1:max_guides] <- resid
        mean_resid <- mean(resid, na.rm = TRUE)
        if (!is.nan(mean_resid)) {
          scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(resid, na.rm = TRUE)
        }
      }
    }
  } else if (loess) {
    cat("Warning: loess-normalization is only enabled for the test=\"moderated-t\" option\n")
  }
  
  # Scores condition response with moderated t-test
  if (test == "moderated-t") {
    for (name in condition_names) {
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]]))
      p_val <- ebayes_fit$p.value[,1]
      scores[[paste0("pval_", name, "_vs_", control_name)]] <- p_val
    }   
  }
  
  # Computes FDRs
  for (name in condition_names) {
    scores[[paste0("fdr_", name, "_vs_", control_name)]] <- 
      stats::p.adjust(scores[[paste0("pval_", name, "_vs_", control_name)]], method = fdr_method)
  }
  
  # Removes extra zero row from residuals
  loess_residuals <- loess_residuals[1:(nrow(loess_residuals) - 1),]
  
  # Explicitly returns scored data
  output <- list()
  output[["scored_data"]] <- scores
  if (return_residuals) {
    output[["residuals"]] <- loess_residuals
  } else {
    output[["residuals"]] <- NA
  }
  return(output)
}

#' Scores conditions against a single control.
#' 
#' Scores guides for any number of condition screens against a multiplicative null 
#' model derived from single-gene effects (e.g. for comparing paralog double-knockouts
#' to single-knockouts for each paralogous gene). Do NOT use this function if you 
#' would like to directly compare conditions against a control. In that case, call 
#' \code{score_conditions_vs_control} instead. After running this function, pass the 
#' resulting dataframe to \code{call_combn_hits} to call significant effects.
#' 
#' @param combn_guides A list of exonic-exonic guides returned from \code{split_guides}.
#' @param single_guides A list of exonic-intergenic guides returned from \code{split_guides}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param screen_names A list of screen names to score against a derived null model from
#'   single-gene effects.
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE). 
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default "BY").
#' @param filter_genes List of genes to filter from scoring (default NULL).
#' @param ignore_orientation If TRUE, aggregates guides across both orientations, returning only
#'   one p-value and FDR column with orientation2 p-values set to NA (default FALSE).
#' @param return_residuals If FALSE, doesn't return residuals dataframe (default TRUE).
#'   This is recommend if scoring large datasets and memory is a limitation.  
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A list containing two dataframes. The first entry, named "scored_data" in the list,
#'   contains scored data with separate columns given by the specified control and condition
#'   names. The second entry, named "residuals" in the list, is a dataframe containing control,
#'   condition and loess-normalized residuals for all guides.
#' @export
score_combn_vs_single <- function(combn_guides, single_guides, screens, screen_names, 
                                  min_guides = 3, test = "moderated-t",
                                  loess = TRUE, fdr_method = "BY", filter_genes = NULL,
                                  ignore_orientation = FALSE, return_residuals = TRUE, 
                                  verbose = FALSE) {
  
  # Filters specified genes from dataset before scoring
  if (!is.null(filter_genes)) {
    for (gene in filter_genes) {
      ind <- unlist(lapply(combn_guides, function(x) x[["gene1"]] == gene | x[["gene2"]] == gene))
      if (sum(ind) > 0) {
        print(paste("Removing", sum(ind), "gene pairs containing", gene, "from combn guides\n"))
        combn_guides <- combn_guides[!ind]
      }
    }
  }
  
  # Gets condition names and columns for any number of conditions
  condition_names <- c()
  condition_cols <- list()
  for (condition in screen_names) {
    condition_names <- c(condition_names, condition)
    rep_names <- screens[[condition]][["replicates"]]
    orient1 <- paste0("orient1_", rep_names)
    orient2 <- paste0("orient2_", rep_names)
    condition_cols[[condition]][["orient1"]] <- orient1
    condition_cols[[condition]][["orient2"]] <- orient2
  }
  
  # Makes output dataframe
  scores <- data.frame(gene1 = rep(NA, length(combn_guides)), gene2 = NA)
  
  # Appends additional columns for each condition
  new_cols <- c()
  for (name in condition_names) {
    new_cols <- c(new_cols, c(
      paste0("n_combn_", name), 
      paste0("n_single_", name),
      paste0("mean_combn_", name),
      paste0("mean_single_", name),
      paste0("var_combn_", name),
      paste0("var_single_", name),
      paste0("orientation_agree_", name),
      paste0("differential_combn_vs_single_", name),
      paste0("pval1_combn_vs_single_", name),
      paste0("pval2_combn_vs_single_", name),
      paste0("fdr1_combn_vs_single_", name),
      paste0("fdr2_combn_vs_single_", name),
      paste0("significant_combn_vs_single_", name)
    ))
  }
  scores[new_cols] <- NA
  
  # Makes residual dataframes if necessary
  max_guides <- -1
  condition_residuals <- list()
  if (test == "moderated-t") {
    
    # Gets max number of guides first
    for (guide in combn_guides) {
      guide_names <- names(guide)[!(names(guide) %in% c("gene1", "gene2", "guide_type"))]
      guide_num <- max(unlist(lapply(guide_names, function(x) length(guide[[x]]))))
      max_guides <- max(max_guides, guide_num)
    }
    
    # Makes residual dataframes with columns equal to the max number of guides.
    # Doubles the size if orientation is ignored
    if (!ignore_orientation) {
      for (name in condition_names) {
        residual_df <- data.frame(matrix(nrow = length(combn_guides), ncol = max_guides))
        colnames(residual_df) <- paste0("guide_residual_", 1:max_guides)
        condition_residuals[[name]][[1]] <- residual_df
        condition_residuals[[name]][[2]] <- residual_df
      } 
    } else {
      for (name in condition_names) {
        residual_df <- data.frame(matrix(nrow = length(combn_guides), ncol = max_guides*2))
        colnames(residual_df) <- paste0("guide_residual_", 1:(max_guides*2))
        condition_residuals[[name]][[1]] <- residual_df
        condition_residuals[[name]][[2]] <- residual_df
      } 
    }
  }
  
  # Makes loess residual dataframes if specified, one for each orientation
  loess_residuals <- list()
  if (return_residuals & test == "moderated-t") {
    loess_residuals[[1]] <- data.frame(n = rep(0, max_guides*length(combn_guides)))
    loess_residuals[[2]] <- data.frame(n = rep(0, max_guides*length(combn_guides)))
    for (name in condition_names) {
      loess_residuals[[1]][[paste0("mean_single_", name)]] <- rep(0, nrow(loess_residuals[[1]]))
      loess_residuals[[1]][[paste0("mean_combn_", name)]] <- rep(0, nrow(loess_residuals[[1]]))
      loess_residuals[[2]][[paste0("mean_single_", name)]] <- rep(0, nrow(loess_residuals[[2]]))
      loess_residuals[[2]][[paste0("mean_combn_", name)]] <- rep(0, nrow(loess_residuals[[2]]))
    }
  }
  
  # Scores guides for each condition
  counter1 <- 1
  counter2 <- 1
  for (i in 1:length(combn_guides)) {
    
    # Gets gene names and guide values
    combn_vals <- combn_guides[[i]]
    gene1 <- combn_vals$gene1
    gene2 <- combn_vals$gene2
    scores$gene1[i] <- gene1
    scores$gene2[i] <- gene2
    
    # Finds matching single-targeting guides if any exist
    single_gene1 <- NULL
    single_gene2 <- NULL
    gene1_ind <- unlist(lapply(single_guides, function(x) x[["gene1"]] == gene1))
    gene2_ind <- unlist(lapply(single_guides, function(x) x[["gene1"]] == gene2))
    if (sum(gene1_ind) > 0 & sum(gene2_ind) > 0) {
      single_gene1 <- single_guides[gene1_ind][[1]]
      single_gene2 <- single_guides[gene2_ind][[1]]
    } else {
      if (verbose) {
        cat(paste(gene1, "and", gene2, "skipped because of too few guides\n"))
      }
      next
    }
    
    # Subsets single-targeting guides to matching IDs if given
    if (!("orient1_id1" %in% names(combn_vals) & "orient1_id2" %in% names(combn_vals))) {
      single_gene1 <- single_gene1[single_gene1[["orient1_id1"]] %in% combn_vals[["orient1_id1"]] |
                                     single_gene1[["orient2_id2"]] %in% combn_vals[["orient2_id2"]]]
      single_gene2 <- single_gene2[single_gene2[["orient1_id1"]] %in% combn_vals[["orient2_id1"]] |
                                     single_gene2[["orient2_id2"]] %in% combn_vals[["orient1_id2"]]]
      combn_vals <- combn_vals[combn_vals[["orient1_id1"]] %in% single_gene1[["orient1_id1"]] &
                                 combn_vals[["orient1_id2"]] %in% single_gene2[["orient2_id2"]]]
      combn_vals <- combn_vals[combn_vals[["orient2_id1"]] %in% single_gene2[["orient1_id1"]] &
                                 combn_vals[["orient2_id2"]] %in% single_gene1[["orient2_id2"]]]
    }
    
    # Scores combn-targeting guides vs. single-targeting null model
    increment1 <- 0
    increment2 <- 0
    for (name in condition_names) {
      
      # Gets column names for each orientation and the current condition
      orient1 <- condition_cols[[name]][["orient1"]]
      orient2 <- condition_cols[[name]][["orient2"]]
      
      # Takes the combn-targeting mean across replicates
      combn1 <- rowMeans(data.frame(combn_vals[orient1]))
      combn2 <- rowMeans(data.frame(combn_vals[orient2]))
      
      # Gets guide values for single-targeting guides
      single_gene1_orient1 <- rowMeans(data.frame(single_gene1[orient1]))
      single_gene1_orient2 <- rowMeans(data.frame(single_gene1[orient2]))
      single_gene2_orient1 <- rowMeans(data.frame(single_gene2[orient1]))
      single_gene2_orient2 <- rowMeans(data.frame(single_gene2[orient2]))
      
      # Gets null model by summing different orientations of single-targeting guides
      null1 <- c()
      null2 <- c()
      if (!("orient1_id1" %in% names(combn_vals) & "orient1_id2" %in% names(combn_vals))) {
        null1 <- unlist(apply(expand.grid(single_gene1_orient1, single_gene2_orient2), 1, sum))
        null2 <- unlist(apply(expand.grid(single_gene1_orient2, single_gene2_orient1), 1, sum))
      } else {
        for (j in 1:length(combn1)) {
          id1 <- combn_vals[["orient1_id1"]][j]
          id2 <- combn_vals[["orient1_id2"]][j]
          val1 <- single_gene1_orient1[single_gene1[["orient1_id1"]] == id1]
          val2 <- single_gene2_orient2[single_gene2[["orient2_id2"]] == id2]
          if (length(val1) > 1) {
            val1 <- mean(val1) 
            if (verbose) {
              cat(paste("Warning: taking mean of single-gene effects for", gene1, "/", gene2, "in ", name, "\n"))
            }
          }
          if (length(val2) > 1) { 
            val2 <- mean(val2) 
            if (verbose) {
              cat(paste("Warning: taking mean of single-gene effects for", gene1, "/", gene2, "in ", name, "\n"))
            }
          }
          if (length(val1) > 0 & length(val2) > 0) {
            null1 <- c(null1, val1 + val2)
          } else {
            null1 <- c(null1, NA)
          }
        }
        for (j in 1:length(combn2)) {
          id1 <- combn_vals[["orient2_id1"]][j]
          id2 <- combn_vals[["orient2_id2"]][j]
          val1 <- single_gene1_orient2[single_gene1[["orient2_id2"]] == id2]
          val2 <- single_gene2_orient1[single_gene2[["orient1_id1"]] == id1]
          if (length(val1) > 1) { 
            val1 <- mean(val1) 
            if (verbose) {
              cat(paste("Warning: taking mean of single-gene effects for", gene1, "/", gene2, "in ", name, "\n"))
            }
          }
          if (length(val2) > 1) { 
            val2 <- mean(val2) 
            if (verbose) {
              cat(paste("Warning: taking mean of single-gene effects for", gene1, "/", gene2, "in ", name, "\n"))
            }
          }
          if (length(val1) > 0 & length(val2) > 0) {
            null2 <- c(null2, val1 + val2)
          } else {
            null2 <- c(null2, NA)
          }
        }
      }
      
      # Subsets combinatorial guides to existing expected values
      keep_ind1 <- !is.na(null1) & !is.nan(null1) & !is.null(null1) &
        !is.na(combn1) & !is.nan(combn1) & !is.null(combn1)
      keep_ind2 <- !is.na(null2) & !is.nan(null2) & !is.null(null2) &
        !is.na(combn2) & !is.nan(combn2) & !is.null(combn2)
      combn1 <- combn1[keep_ind1]
      combn2 <- combn2[keep_ind2]
      null1 <- null1[keep_ind1]
      null2 <- null2[keep_ind2]
      
      # Joins guides if orientation is ignored
      all_combn <- NULL
      all_null <- NULL
      if (ignore_orientation) {
        all_combn <- c(combn1, combn2)
        all_null <- c(null1, null2) 
      }
      
      # Skips if too few guides
      if (length(null1) < min_guides | length(null2) < min_guides |
          length(combn1) < min_guides | length(combn2) < min_guides) {
        if (verbose) {
          cat(paste(gene1, "and", gene2, "in", name, "skipped because of too few guides\n"))
        }
        next
      } 
      
      # Appends values for loess-normalization to residual dataframe if necessary
      ind1 <- counter1:(counter1 + length(combn1) - 1)
      ind2 <- counter2:(counter2 + length(combn2) - 1)
      if (return_residuals & test == "moderated-t") {
        loess_residuals[[1]][["n"]][ind1] <- i
        loess_residuals[[2]][["n"]][ind2] <- i
        loess_residuals[[1]][[paste0("mean_single_", name)]][ind1] <- null1
        loess_residuals[[1]][[paste0("mean_combn_", name)]][ind1] <- combn1
        loess_residuals[[2]][[paste0("mean_single_", name)]][ind2] <- null2
        loess_residuals[[2]][[paste0("mean_combn_", name)]][ind2] <- combn2
      }
      
      # Gets residuals
      scores[[paste0("orientation_agree_", name)]][i] <- 
        (sign(mean(combn1) - mean(null1)) == sign(mean(combn2) - mean(null2)))
      scores[[paste0("differential_combn_vs_single_", name)]][i] <- mean(c(combn1, combn2)) - mean(c(null1, null2))
      
      # Performs the specified type of testing or stores residuals for later testing
      if (test == "rank-sum") {
        if(!ignore_orientation) {
          scores[[paste0("pval1_combn_vs_single_", name)]][i] <- 
            suppressWarnings(stats::wilcox.test(combn1, null1))$p.value
          scores[[paste0("pval2_combn_vs_single_", name)]][i] <- 
            suppressWarnings(stats::wilcox.test(combn2, null2))$p.value  
        } else {
          scores[[paste0("pval1_combn_vs_single_", name)]][i] <-
            suppressWarnings(stats::wilcox.test(all_combn, all_null))$p.value
        }
      } else if (test == "moderated-t") {
        residuals1 <- combn1 - null1
        residuals2 <- combn2 - null2
        if(length(residuals1) < max_guides) { 
          residuals1 <- c(residuals1, rep(NA, max_guides - length(residuals1))) 
        }
        if(length(residuals2) < max_guides) { 
          residuals2 <- c(residuals2, rep(NA, max_guides - length(residuals2))) 
        }
        condition_residuals[[name]][[1]][i,1:max_guides] <- residuals1 
        condition_residuals[[name]][[2]][i,1:max_guides] <- residuals2 
      }
      
      # Stores number of guides in dataframe depending on whether orientation is ignored
      if (!ignore_orientation) {
        scores[[paste0("n_combn_", name)]][i] <- length(combn1)
        scores[[paste0("n_single_", name)]][i] <- length(null1) 
      } else {
        scores[[paste0("n_combn_", name)]][i] <- length(c(combn1, combn2))
        scores[[paste0("n_single_", name)]][i] <- length(c(null1, null2)) 
      }

      # Stores values in dataframe
      scores[[paste0("mean_combn_", name)]][i] <- mean(c(combn1, combn2))
      scores[[paste0("mean_single_", name)]][i] <- mean(c(null1, null2))
      scores[[paste0("var_combn_", name)]][i] <- mean(c(stats::var(combn1), stats::var(combn2)))
      scores[[paste0("var_single_", name)]][i] <- mean(c(stats::var(null1), stats::var(null2)))
      
      # Increments counter after looping through all conditions
      increment1 <- length(combn1)
      increment2 <- length(combn2)
    }
    counter1 <- counter1 + increment1
    counter2 <- counter2 + increment2
  }
  
  # Computes loess-normalized residuals if specified
  joined_residuals <- NULL
  if (loess & test == "moderated-t") {
    
    # Normalizes orientations separately if orientation taken into account
    loess_residuals[[1]] <- loess_residuals[[1]][1:counter1,]
    loess_residuals[[2]] <- loess_residuals[[2]][1:counter2,]
    if (!ignore_orientation) {
      for (name in condition_names) {
        for (i in 1:2) {
          null_values <- loess_residuals[[i]][[paste0("mean_single_", name)]]
          condition_values <- loess_residuals[[i]][[paste0("mean_combn_", name)]]
          temp <- loess_MA(null_values, condition_values)
          loess_residuals[[i]][[paste0("loess_residual_", name)]] <- temp[["residual"]]
          loess_residuals[[i]][[paste0("loess_predicted_", name)]] <- temp[["predicted"]] 
        }
      } 
      
    # Otherwise we join residuals before loess-normalization
    } else {
      joined_residuals <- rbind(loess_residuals[[1]], loess_residuals[[2]])
      joined_residuals <- joined_residuals[order(joined_residuals[["n"]]),]
      for (name in condition_names) {
        null_values <- joined_residuals[[paste0("mean_single_", name)]]
        condition_values <- joined_residuals[[paste0("mean_combn_", name)]]
        temp <- loess_MA(null_values, condition_values)
        joined_residuals[[paste0("loess_residual_", name)]] <- temp[["residual"]]
        joined_residuals[[paste0("loess_predicted_", name)]] <- temp[["predicted"]] 
      } 
    }
    
    # Replaces residuals with loess-normalized residuals. Only runs on 
    # joined residuals if orientation is ignored
    for (i in 1:length(combn_guides)) {
      if (!ignore_orientation) {
        for (name in condition_names) {
          ind1 <- loess_residuals[[1]]$n == i
          ind2 <- loess_residuals[[2]]$n == i
          resid1 <- loess_residuals[[1]][[paste0("loess_residual_", name)]][ind1]
          predicted1 <- loess_residuals[[1]][[paste0("loess_predicted_", name)]][ind1]
          resid2 <- loess_residuals[[2]][[paste0("loess_residual_", name)]][ind2]
          predicted2 <- loess_residuals[[2]][[paste0("loess_predicted_", name)]][ind2]
          if (length(resid1) < max_guides) { 
            resid1 <- c(resid1, rep(NA, max_guides - length(resid1))) 
          } 
          if (length(resid2) < max_guides) { 
            resid2 <- c(resid2, rep(NA, max_guides - length(resid2))) 
          } 
          condition_residuals[[name]][[1]][i,1:max_guides] <- resid1
          condition_residuals[[name]][[2]][i,1:max_guides] <- resid2
          mean_resid <- mean(c(resid1, resid2), na.rm = TRUE)
          if (!is.nan(mean_resid)) {
            scores[[paste0("differential_combn_vs_single_", name)]][i] <- mean_resid
          }
        }
      } else {
        ind <- joined_residuals$n == i
        resid <- joined_residuals[[paste0("loess_residual_", name)]][ind]
        predicted <- joined_residuals[[paste0("loess_predicted_", name)]][ind]
        if (length(resid) < max_guides*2) { 
          resid <- c(resid, rep(NA, max_guides*2 - length(resid))) 
        } 
        condition_residuals[[name]][[1]][i,1:(max_guides*2)] <- resid
        mean_resid <- mean(resid, na.rm = TRUE)
        if (!is.nan(mean_resid)) { 
          scores[[paste0("differential_combn_vs_single_", name)]][i] <- mean_resid
        }
      }
    }
  } else if (loess) {
    cat("Warning: loess-normalization is only enabled for the test=\"moderated-t\" option\n")
  }
  
  # Scores condition response with moderated t-test if specified
  if (test == "moderated-t") {
    if (!ignore_orientation) {
      for (name in condition_names) {
        ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]][[1]]))
        p_val1 <- ebayes_fit$p.value[,1]
        ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]][[2]]))
        p_val2 <- ebayes_fit$p.value[,1]
        scores[[paste0("pval1_combn_vs_single_", name)]] <- p_val1
        scores[[paste0("pval2_combn_vs_single_", name)]] <- p_val2
      } 
    } else {
      
      # Appends both conditions together if loess-normalization not used
      if(!loess) {
        for (name in condition_names) {
          colnames(condition_residuals[[name]][[2]]) <- paste0(colnames(condition_residuals[[name]][[2]]), "_orientation2")
          condition_residuals[[name]][[1]][,(max_guides + 1):ncol(condition_residuals[[name]][[1]])] <-
            condition_residuals[[name]][[2]][,1:max_guides]
        }
      }
      for (name in condition_names) {
        resid <- condition_residuals[[name]][[1]]
        ebayes_fit <- limma::eBayes(limma::lmFit(resid))
        p_val <- ebayes_fit$p.value[,1]
        scores[[paste0("pval1_combn_vs_single_", name)]] <- p_val
        scores[[paste0("pval2_combn_vs_single_", name)]] <- NA
      } 
    }
  }
  
  # Gets FDRs
  for (name in condition_names) {
    if (!ignore_orientation) {
      scores[[paste0("fdr1_combn_vs_single_", name)]] <- 
        stats::p.adjust(scores[[paste0("pval1_combn_vs_single_", name)]], method = fdr_method)
      scores[[paste0("fdr2_combn_vs_single_", name)]] <- 
        stats::p.adjust(scores[[paste0("pval2_combn_vs_single_", name)]], method = fdr_method)
    } else {
      scores[[paste0("fdr1_combn_vs_single_", name)]] <- 
        stats::p.adjust(scores[[paste0("pval1_combn_vs_single_", name)]], method = fdr_method)
      scores[[paste0("fdr2_combn_vs_single_", name)]] <- NA
    }
  }
  
  # Removes extra zero row from residuals
  loess_residuals[[1]] <- loess_residuals[[1]][1:(nrow(loess_residuals[[1]]) - 1),]
  loess_residuals[[2]] <- loess_residuals[[2]][1:(nrow(loess_residuals[[2]]) - 1),]
  
  # Explicitly returns scored data
  output <- list()
  output[["scored_data"]] <- scores
  if (return_residuals) {
    output[["residuals"]] <- loess_residuals
  } else {
    output[["residuals"]] <- NA
  }
  return(output)
}

#' Call significant responses for scored data.
#' 
#' Run this to call significant responses for data returned from 
#' \code{score_conditions_vs_control}. Do NOT run this on data returned from
#' \code{score_combn_vs_single}.
#' 
#' @param scores Dataframe returned from \code{score_conditions_vs_control}.
#' @param control_screen_name Name of a control screen to test condition screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_condition_hits <- function(scores, control_screen_name, condition_screen_names,
                                fdr_threshold = 0.1, differential_threshold = 0.5,
                                neg_type = "Negative", pos_type = "Positive") {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_screen_name
  condition_names <- c()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
  }
  
  # Calls significant differences for each condition against the control
  for (name in condition_names) {
    scores[[paste0("significant_", name, "_vs_", control_name)]] <- 
      scores[[paste0("fdr_", name, "_vs_", control_name)]] < fdr_threshold
  }
  
  # Makes thresholded calls for significant negative and positive effects
  for (name in condition_names) {
    response_col <- paste0("effect_type_", name)
    scores[[response_col]] <- "None"
    diffs <- scores[[paste0("differential_", name, "_vs_", control_name)]]
    sig <- scores[[paste0("significant_", name, "_vs_", control_name)]]
    scores[[response_col]][sig & diffs < 0 & abs(diffs) > differential_threshold] <- neg_type
    scores[[response_col]][sig & diffs > 0 & abs(diffs) > differential_threshold] <- pos_type
  }
  
  # Explicitly returns scored data
  return(scores)
}

#' Call significant responses for scored paired data.
#' 
#' Run this to call significant responses for data returned from 
#' \code{score_combn_vs_single}. Do NOT run this on data returned from
#' \code{score_conditions_vs_control}.
#' 
#' @param scores Dataframe returned from \code{score_combn_vs_single}.
#' @param screen_names A list of screen names scored in \code{score_combn_vs_single}
#'   for which to call significant effects.
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_combn_hits <- function(scores, screen_names,
                            fdr_threshold = 0.1, differential_threshold = 0.5,
                            neg_type = "Negative", pos_type = "Positive") {
  
  # Gets condition names and columns for any number of conditions
  condition_names <- c()
  condition_cols <- list()
  for (condition in screen_names) {
    condition_names <- c(condition_names, condition)
  }
  
  # Checks if orientation was ignored
  ignore_orientation <- TRUE
  for (name in condition_names) {
    pval2_col <- paste0("pval2_combn_vs_single_", name)
    if (sum(is.na(scores[[pval2_col]])) != nrow(scores)) {
      ignore_orientation <- FALSE
    } 
  }
  
  # Calls significant differences for dual-targeting vs. single-targeting null model
  if (!ignore_orientation) {
    for (name in condition_names) {
      scores[[paste0("significant_combn_vs_single_", name)]] <- 
        scores[[paste0("fdr1_combn_vs_single_", name)]] < fdr_threshold &
        scores[[paste0("fdr2_combn_vs_single_", name)]] < fdr_threshold &
        scores[[paste0("orientation_agree_", name)]]
    } 
  } else {
    for (name in condition_names) {
      scores[[paste0("significant_combn_vs_single_", name)]] <- 
        scores[[paste0("fdr1_combn_vs_single_", name)]] < fdr_threshold
    } 
  }
  
  # Makes thresholded calls for significant negative and positive effects
  for (name in condition_names) {
    response_col <- paste0("effect_type_", name)
    scores[[response_col]] <- "None"
    diffs <- scores[[paste0("differential_combn_vs_single_", name)]]
    sig <- scores[[paste0("significant_combn_vs_single_", name)]]
    scores[[response_col]][sig & diffs < 0 & abs(diffs) > differential_threshold] <- neg_type
    scores[[response_col]][sig & diffs > 0 & abs(diffs) > differential_threshold] <- pos_type
  }
  
  # Explicitly returns scored data
  return(scores)
}

# Fits a loess curve to predict y given x
loess_MA <- function(x, y, sp = 0.4, dg = 2, binSize = 100) {
  #this concept is based on pythagoras and cancels out sqrt, square and factor 2
  #it also ignores the factor sqrt(2) as factor between y and x vs distance of x,y from diagonal x = y
  if(all(x == y, na.rm = T)) { #if e.g. wt scored against itself
    gi <- rep(NA, length(x))
  }
  else {
    m <- y - x
    a <- y + x
    A <- (a - stats::median(a, na.rm = T)) / stats::mad(a, na.rm = T) #scale to generate bins along m
    B <- seq(floor(min(A, na.rm = T)), ceiling(max(A, na.rm = T)), .1) #define bins
    b <- c() #indices for model training
    for(i in 1:(length(B) - 1)) {
      temp_b <- which(A > B[i] & A < B[i+1])
      if(length(temp_b) > binSize) { #sample if more events in bin than max bin size
        set.seed(1)
        temp_b <- sample(temp_b, binSize)
      }
      b <- c(b, temp_b)
    }
    I <- is.finite(m[b]) & is.finite(a[b]) #only use finite values
    model <- stats::loess(m[b][I] ~ a[b][I], span = sp, degree = dg) #train model on m ~ a (approx. y ~ x)
    expected <- stats::predict(model, a) #predict expected m ~ a
    gi <- m - expected
  }
  result <- list()
  result[["residual"]] <- gi
  result[["predicted"]] <- expected
  return(result)
}

#' Scores multiple conditions against multiple controls
#' 
#' Takes in an input .tsv file with two columns for "Screen" and "Control" and scores
#' all screens listed in "Screen" against their corresponding screens listed in
#' "Control." Outputs all files and plots in the specified folder. 
#' 
#' @param guides A list of guides returned from \code{split_guides}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param batch_table Either a dataframe or a path to .tsv file mapping screens to their controls 
#'   for scoring, with two columns for "Screen" and "Control." Screens to score against derived 
#'   null-models with the combn scoring mode must have their respective control labeled as "combn." 
#' @param output_folder Folder to output scored data and plots to. 
#' @param separate_orientation If true, then guide values are scored separately across each 
#'   orientation (default FALSE).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param filter_genes List of genes to filter from scoring (default NULL).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
score_conditions_batch <- function(guides, screens, batch_table, output_folder, 
                                   separate_orientation = FALSE, min_guides = 3, test = "moderated-t", 
                                   loess = TRUE, filter_genes = NULL, fdr_method = "BY",
                                   fdr_threshold = 0.1, differential_threshold = 0.5, 
                                   neg_type = "Negative", pos_type = "Positive",
                                   plot_type = "png") {
  
  # Checks batch file and loads it
  batch <- load_batch_table(batch_table, screens)
  
  # Makes output folders if nonexistent
  lfc_folder <- file.path(output_folder, "lfc")
  plot_folder <- file.path(output_folder, "plots")
  if (!dir.exists(output_folder)) { dir.create(output_folder) }
  if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
  if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
  
  # Scores guides for each batch
  all_scores <- NULL
  for (i in 1:nrow(batch)) {
    condition <- batch[i,1]
    control <- batch[i,2]
    if (control != "combn") {
      screen_lfc_folder <- file.path(lfc_folder, paste0(condition, "_", control))
      temp <- score_conditions_vs_control(guides, screens, control, condition, test = test, 
                                          min_guides = min_guides, loess = loess, 
                                          filter_genes = NULL, fdr_method = fdr_method)
      scores <- temp[["scored_data"]]
      residuals <- temp[["residuals"]]
      scores <- call_condition_hits(scores, control, condition,
                                    neg_type = neg_type, pos_type = pos_type,
                                    fdr_threshold = fdr_threshold, 
                                    differential_threshold = differential_threshold)
      plot_condition_residuals(scores, residuals, control, condition, screen_lfc_folder, 
                               neg_type = neg_type, pos_type = pos_type,
                               plot_type = plot_type)
      plot_condition_response(scores, control, condition, plot_folder,
                              neg_type = neg_type, pos_type = pos_type,
                              plot_type = plot_type)
      if (is.null(all_scores)) {
        all_scores <- scores
      } else {
        all_scores <- cbind(all_scores, scores[,3:ncol(scores)])
      }
    }
  }
  if (!is.null(all_scores)) {
    utils::write.table(all_scores, file.path(output_folder, "condition_gene_calls.tsv"), sep = "\t",
                       row.names = FALSE, col.names = TRUE, quote = FALSE) 
  }
}

#' Scores multiple combinatorial screens against derived null models
#' 
#' Takes in an input .tsv file with two columns for "Screen" and "Control" and scores
#' all screens listed in "Screen" against their corresponding screens listed in
#' "Control." Outputs all files and plots in the specified folder. 
#' 
#' @param combn_guides A list of exonic-exonic guides returned from \code{split_guides}.
#' @param single_guides A list of exonic-intergenic guides returned from \code{split_guides}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param batch_table Either a dataframe or a path to .tsv file mapping screens to their controls 
#'   for scoring, with two columns for "Screen" and "Control." Screens to score against derived 
#'   null-models with the combn scoring mode must have their respective control labeled as "combn." 
#' @param output_folder Folder to output scored data and plots to. 
#' @param separate_orientation If true, then guide values are scored separately across each 
#'   orientation (default FALSE).
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param filter_genes List of genes to filter from scoring (default NULL).
#' @param ignore_orientation If TRUE, aggregates guides across both orientations, returning only
#'   one p-value and FDR column with orientation2 p-values set to NA (default FALSE).
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
score_combn_batch <- function(combn_guides, single_guides, screens, batch_table, output_folder, 
                              separate_orientation = FALSE, min_guides = 3, test = "moderated-t", 
                              loess = TRUE, filter_genes = NULL, ignore_orientation = FALSE, 
                              fdr_method = "BY", fdr_threshold = 0.1, differential_threshold = 0.5, 
                              neg_type = "Negative", pos_type = "Positive", plot_type = "png") {
  
  # Checks batch file and loads it
  batch <- load_batch_table(batch_table, screens)
  
  # Makes output folders if nonexistent
  lfc_folder <- file.path(output_folder, "lfc")
  plot_folder <- file.path(output_folder, "plots")
  if (!dir.exists(output_folder)) { dir.create(output_folder) }
  if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
  if (!dir.exists(plot_folder)) { dir.create(plot_folder) }
  
  # Scores guides for each batch
  all_scores <- NULL
  for (i in 1:nrow(batch)) {
    condition <- batch[i,1]
    control <- batch[i,2]
    if (control == "combn") {
      screen_lfc_folder <- file.path(lfc_folder, paste0(condition, "_combn"))
      temp <- score_combn_vs_single(combn_guides, single_guides, screens, condition, test = test, 
                                    min_guides = min_guides, loess = loess, filter_genes = NULL,
                                    ignore_orientation = ignore_orientation, fdr_method = fdr_method)
      scores <- temp[["scored_data"]]
      residuals <- temp[["residuals"]]
      scores <- call_combn_hits(scores, condition,
                                neg_type = neg_type, pos_type = pos_type,
                                fdr_threshold = fdr_threshold, 
                                differential_threshold = differential_threshold)
      plot_combn_residuals(scores, residuals, condition, screen_lfc_folder, 
                           neg_type = neg_type, pos_type = pos_type,
                           plot_type = plot_type)
      plot_combn_response(scores, condition, plot_folder,
                          neg_type = neg_type, pos_type = pos_type,
                          plot_type = plot_type)
      if (is.null(all_scores)) {
        all_scores <- scores
      } else {
        all_scores <- cbind(all_scores, scores[,3:ncol(scores)])
      }
    }
  }
  if (!is.null(all_scores)) {
    utils::write.table(all_scores, file.path(output_folder, "combn_gene_calls.tsv"), sep = "\t",
                       row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
}

#' All-in-one wrapper for Orthrus
#' 
#' Takes in three .tsv files for per-replicate readcounts, a mapping of screens to replicates, 
#' and a list of screens to score, and scores all given screens before outputting results to a 
#' given folder. The format for each .tsv file is described below.
#' 
#' reads_file: A tab-separated file with readcounts for each guide and each replicate. In addition
#'   to numerical readcounts columns, genes must be labeled in columns named "gene1" and "gene2."
#'   Optionally, the dataset may also contain columns with guide-specific IDs, such as guide
#'   sequences, which the user can specify the names of in the parameters id_col1 and id_col2. 
#'   These are required to compute matching guide-specific residuals when test is specified
#'   as "moderated-t."
#' sample_file: 
#' batch_file:
#' 
#' @param reads_file Path to .tsv file where each row corresponds to a single guide construct that 
#'   targets two genomic locations, which contains the columns gene1, gene2 and optionally 
#'   the columns specified by id_col1 and id_col2.
#' @param sample_file Path to .tsv file with three columns: Screen, Replicates, and NormalizeTo. 
#'   Screen is the name of the screen, replicates are each technical replicate separated by
#'   semicolons, and NormalizeTo is either the screen to normalize against or NA if unnecessary 
#'   (e.g. for T0 screens).
#' @param batch_file Path to .tsv file mapping screens to their controls for scoring, with two 
#'   columns for "Screen" and "Control." Screens to score against dervied null-models with the
#'   combn scoring mode must have their respective control labeled as "combn." 
#' @param output_folder Folder to output scored data and plots to.
#' @param id_col1 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in the column gene1 (optional, default NULL).  
#' @param id_col2 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in the column gene2 (optional, default NULL).  
#' @param filter_names List of screen names to filter based on read counts by. 
#' @param filter_genes List of genes to filter from scoring (default NULL).
#' @param min_reads Minimum number of reads to keep (default 30, anything
#'   below this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (default 10000, anything
#'   above this value will be filtered out).
#' @param display_numbers Whether or not to include PCC values in heatmap (default TRUE).
#' @param negative_controls List of negative control genes to append to default list of
#'   non-essential genes (default NULL).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE).
#' @param ignore_orientation If TRUE, aggregates guides across both orientations, returning only
#'   one p-value and FDR column with orientation2 p-values set to NA (default FALSE).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @param fdr_method Type of FDR to compute. One of "BH", "BY" or "bonferroni" (default
#'   "BY")
#' @param fdr_threshold Threshold below which to call gene effects as significant 
#'   (default 0.1).
#' @param differential_threshold Absolute value threshold on differential effects, 
#'   below which gene effects are not called as significant (default 0.5).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
orthrus_wrapper <- function(reads_file, sample_file, batch_file, output_folder,
                            id_col1 = NULL, id_col2 = NULL, 
                            filter_names = NULL, filter_genes = NULL,
                            min_reads = 30, max_reads = 10000,
                            display_numbers = FALSE, negative_controls = NULL, 
                            test = "moderated-t", loess = TRUE,
                            ignore_orientation = FALSE,
                            neg_type = "Negative", pos_type = "Positive", 
                            fdr_method = "BY", fdr_threshold = 0.1, 
                            differential_threshold = 0.5, plot_type = "png") {
  
  # Makes output folders if nonexistent
  qc_folder <- file.path(output_folder, "qc")
  lfc_folder <- file.path(qc_folder, "lfc")
  if (!dir.exists(output_folder)) { dir.create(output_folder, recursive = TRUE) }
  if (!dir.exists(qc_folder)) { dir.create(qc_folder) }
  if (!dir.exists(lfc_folder)) { dir.create(lfc_folder) }
  
  # Reads in data
  df <- utils::read.csv(reads_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  screens <- add_screens_from_table(sample_file)
  
  # Processes screens and generates QC plots and metrics
  plot_reads_qc(df, screens, qc_folder, display_numbers = display_numbers, plot_type = plot_type)
  df <- normalize_screens(df, screens, filter_names = filter_names, min_reads = min_reads,
                          max_reads = max_reads)
  plot_lfc_qc(df, screens, qc_folder, display_numbers = display_numbers, plot_type = plot_type, 
              negative_controls = negative_controls)
  guides <- split_guides(df, screens, id_col1, id_col2)
  dual <- guides[["single_gene_dual_targeted"]]
  single <- guides[["exonic_intergenic"]]
  paralogs <- guides[["exonic_exonic"]]
  
  # Scores dual-targeting guides
  score_conditions_batch(dual, screens, batch_file, output_folder, 
                         test = test, loess = loess, filter_genes = filter_genes,
                         neg_type = neg_type, pos_type = pos_type,
                         fdr_threshold = fdr_threshold, 
                         differential_threshold = differential_threshold,
                         plot_type = plot_type)
  
  # Scores combinatorial-targeting guides
  score_combn_batch(paralogs, single, screens, batch_file, output_folder, 
                    test = test, loess = loess, filter_genes = filter_genes,
                    neg_type = neg_type, pos_type = pos_type,
                    fdr_threshold = fdr_threshold, 
                    differential_threshold = differential_threshold,
                    plot_type = plot_type,
                    ignore_orientation = ignore_orientation)
  
}