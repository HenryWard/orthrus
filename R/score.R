######
# SCORING CODE
######

#' Scores conditions against a single control.
#' 
#' Scores guides for any number of condition screens against a control screen
#' (e.g. for directly comparing drug response to DMSO response). Do NOT use 
#' this function if you would like to compare guides against a null model 
#' derived from single-gene effects. In that case, call \code{score_combn_vs_single}
#' instead. After running this function, pass the resulting dataframe to 
#' \code{call_significant_response} to call significant effects.
#' 
#' @param guides A list of guides returned from \code{retrieve_guides_by_label} or 
#'   from an entry of \code{split_guides_by_type}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param control_screen_name Name of a control screen to test condition screens against.
#' @param condition_screen_names A list of condition screen names to score against the 
#'   control screen.
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE). 
#' @return A dataframe of scored data with separate columns given by the specified control
#'   and condition names.
#' @export
score_conditions_vs_control <- function(guides, screens, control_screen_name, condition_screen_names, 
                                        min_guides = 3, test = "moderated-t", loess = TRUE) {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_screen_name
  control_cols <- screens[[control_name]][["replicates"]]
  condition_names <- c()
  condition_cols <- list()
  for (condition in condition_screen_names) {
    condition_names <- c(condition_names, condition)
    condition_cols[[condition]] <- screens[[condition]][["replicates"]]
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
  if (loess & test == "moderated-t") {
    loess_residuals <- data.frame(n = rep(0, max_guides*length(guides)))
    loess_residuals[[paste0("mean_", control_name)]] <- rep(0, nrow(loess_residuals))
    for (name in condition_names) {
      loess_residuals[[paste0("mean_", name)]] <- rep(0, nrow(loess_residuals))
    }
  }
  
  # Appends additional columns for each condition
  new_cols <- c(paste0("n_", control_name), 
                paste0("mean_", control_name))
  for (name in condition_names) {
    new_cols <- c(new_cols, c(
      paste0("n_", name), 
      paste0("mean_", name),
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
    
    # Makes loess-normalized residual dataframe if necessary
    ind <- counter:(counter + length(rep_mean_control) - 1)
    if (loess) {
      loess_residuals$n[ind] <- i
      loess_residuals[[paste0("mean_", control_name)]][ind] <- rep_mean_control
    }
    
    # Takes the mean across replicates for all conditions
    for (name in condition_names) {
      rep_mean_condition <- rowMeans(data.frame(guide_vals[condition_cols[[name]]]))
      diff <- rep_mean_condition - rep_mean_control
      scores[[paste0("n_", control_name)]][i] <- length(rep_mean_control)
      scores[[paste0("n_", name)]][i] <- length(rep_mean_condition)
      scores[[paste0("mean_", control_name)]][i] <- mean(rep_mean_control)
      scores[[paste0("mean_", name)]][i] <- mean(rep_mean_condition)
      scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff)
      
      # Appends mean LFCs for loess-normalization if specified
      if (loess) {
        loess_residuals[[paste0("mean_", name)]][ind] <- rep_mean_condition
      }
      
      # Performs the specified type of testing or stores residuals for later testing
      if (test == "rank-sum") {
        scores[[paste0("pval_", name, "_vs_", control_name)]][i] <- 
          suppressWarnings(wilcox.test(rep_mean_condition, rep_mean_control))$p.value
      } else if (test == "moderated-t") {
        if (length(diff) < max_guides) { diff <- c(diff, rep(NA, max_guides - length(diff))) } 
        condition_residuals[[name]][i,1:max_guides] <- diff 
      }
    }
    counter <- counter + length(rep_mean_control)
  }
   
  # Computes loess-normalized residuals if specified
  if (loess & test == "moderated-t") {
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
        scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(resid, na.rm = TRUE)
        # scores[[paste0("loess_predicted_", name)]][i] <- mean(predicted)
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
      p.adjust(scores[[paste0("pval_", name, "_vs_", control_name)]], method = "BH")
  }
  
  # Removes genes with too few observations
  scores <- scores[scores[[paste0("n_", control_name)]] >= min_guides,]
  
  # Explicitly returns scored data
  return(scores)
}

#' Scores conditions against a single control.
#' 
#' Scores guides for any number of condition screens against a multiplicative null 
#' model derived from single-gene effects (e.g. for comparing paralog double-knockouts
#' to single-knockouts for each paralogous gene). Do NOT use this function if you 
#' would like to directly compare conditions against a control. In that case, call 
#' \code{score_conditions_vs_control} instead. After running this function, pass the 
#' resulting dataframe to \code{call_significant_response_dual} to call significant 
#' effects.
#' 
#' @param combn_guides A list of exonic-exonic guides returned from \code{retrieve_guides_by_label} 
#'   or from the exonic-exonic entry of \code{split_guides_by_type}.
#' @param single_guides A list of exonic-intergenic guides returned from \code{retrieve_guides_by_label} 
#'   or from the exonic-intergenic entry of \code{split_guides_by_type}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param screen_names A list of screen names to score against a derived null model from
#'   single-gene effects.
#' @param test Type of hypothesis testing to run. Must be one of "rank-sum" for Wilcoxon
#'   rank-sum testing or "moderated-t" for moderated t-testing (default "moderated-t").
#' @param loess If true, loess-normalizes residuals before running hypothesis testing.
#'   Only works when test = "moderated-t" (default TRUE). 
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @param verbose If true, prints verbose output (default FALSE). 
#' @return A dataframe of scored data with separate columns given by the specified control
#'   and condition names.
#' @export
score_combn_vs_single <- function(combn_guides, single_guides, screens, screen_names, 
                                  min_guides = 3, test = "moderated-t",
                                  loess = TRUE, verbose = FALSE) {
  
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
    
    # Makes residual dataframes with columns equal to the max number of guides
    for (name in condition_names) {
      residual_df <- data.frame(matrix(nrow = length(combn_guides), ncol = max_guides))
      colnames(residual_df) <- paste0("guide_residual_", 1:max_guides)
      condition_residuals[[name]][[1]] <- residual_df
      condition_residuals[[name]][[2]] <- residual_df
    }
  }
  
  # Makes loess residual dataframes if specified, one for each orientation
  loess_residuals <- list()
  if (loess & test == "moderated-t") {
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
      combn1 <- combn1[!is.na(null1)]
      combn2 <- combn2[!is.na(null2)]
      null1 <- null1[!is.na(null1)]
      null2 <- null2[!is.na(null2)]
      
      # Skips if too few guides
      if (length(null1) < 2 | length(null2) < 2 |
          length(combn1) < 2 | length(combn2) < 2) {
        if (verbose) {
          cat(paste(gene1, "and", gene2, "in", name, "skipped because of too few guides\n"))
        }
        next
      }
      
      # Appends values for loess-normalization to residual dataframe if necessary
      ind1 <- counter1:(counter1 + length(combn1) - 1)
      ind2 <- counter2:(counter2 + length(combn2) - 1)
      if (loess & test == "moderated-t") {
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
        scores[[paste0("pval1_combn_vs_single_", name)]][i] <- suppressWarnings(wilcox.test(combn1, null1))$p.value
        scores[[paste0("pval2_combn_vs_single_", name)]][i] <- suppressWarnings(wilcox.test(combn2, null2))$p.value 
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

      # Stores values in dataframe
      scores[[paste0("n_combn_", name)]][i] <- length(combn1)
      scores[[paste0("n_single_", name)]][i] <- length(null1)
      scores[[paste0("mean_combn_", name)]][i] <- mean(c(combn1, combn2))
      scores[[paste0("mean_single_", name)]][i] <- mean(c(null1, null2))
      
      # Increments counter after looping through all conditions
      increment1 <- length(combn1)
      increment2 <- length(combn2)
    }
    counter1 <- counter1 + increment1
    counter2 <- counter2 + increment2
  }
  
  # Computes loess-normalized residuals if specified
  if (loess & test == "moderated-t") {
    loess_residuals[[1]] <- loess_residuals[[1]][1:counter1,]
    loess_residuals[[2]] <- loess_residuals[[2]][1:counter2,]
    for (name in condition_names) {
      for (i in 1:2) {
        null_values <- loess_residuals[[i]][[paste0("mean_single_", name)]]
        condition_values <- loess_residuals[[i]][[paste0("mean_combn_", name)]]
        temp <- loess_MA(null_values, condition_values)
        loess_residuals[[i]][[paste0("loess_residual_", name)]] <- temp[["residual"]]
        loess_residuals[[i]][[paste0("loess_predicted_", name)]] <- temp[["predicted"]] 
      }
    }
    
    # Replaces residuals with loess-normalized residuals
    for (i in 1:length(combn_guides)) {
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
        scores[[paste0("differential_combn_vs_single_", name)]][i] <- mean(c(resid1, resid2), na.rm = TRUE)
        #scores[[paste0("loess_predicted_", name)]][i] <- mean(c(predicted1, predicted2))
      }
    }
  } else if (loess) {
    cat("Warning: loess-normalization is only enabled for the test=\"moderated-t\" option\n")
  }
  
  # Scores condition response with moderated t-test if specified
  if (test == "moderated-t") {
    for (name in condition_names) {
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]][[1]]))
      p_val1 <- ebayes_fit$p.value[,1]
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]][[2]]))
      p_val2 <- ebayes_fit$p.value[,1]
      scores[[paste0("pval1_combn_vs_single_", name)]] <- p_val1
      scores[[paste0("pval2_combn_vs_single_", name)]] <- p_val2
    }
  }
  
  # Gets FDRs
  for (name in condition_names) {
    scores[[paste0("fdr1_combn_vs_single_", name)]] <- 
      p.adjust(scores[[paste0("pval1_combn_vs_single_", name)]], method = "BH")
    scores[[paste0("fdr2_combn_vs_single_", name)]] <- 
      p.adjust(scores[[paste0("pval2_combn_vs_single_", name)]], method = "BH")
  }
  
  # Explicitly returns scoresd data
  return(scores)
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
#'   below which gene effects are not called as significant (default 0).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_significant_response <- function(scores, control_screen_name, condition_screen_names,
                                      fdr_threshold = 0.1, differential_threshold = 0,
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
#'   below which gene effects are not called as significant (default 0).
#' @param neg_type Label for significant effects with a negative differential effect
#'   (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   (default "Positive").
#' @return Dataframe of scored data with differential effects called as significant
#'   for the specified conditions. 
#' @export
call_significant_response_combn <- function(scores, screen_names,
                                            fdr_threshold = 0.1, differential_threshold = 0,
                                            neg_type = "Negative", pos_type = "Positive") {
  
  # Gets condition names and columns for any number of conditions
  condition_names <- c()
  condition_cols <- list()
  for (condition in screen_names) {
    condition_names <- c(condition_names, condition)
  }
  
  # Calls significant differences for dual-targeting vs. single-targeting null model
  for (name in condition_names) {
    scores[[paste0("significant_combn_vs_single_", name)]] <- 
      scores[[paste0("fdr1_combn_vs_single_", name)]] < fdr_threshold &
      scores[[paste0("fdr2_combn_vs_single_", name)]] < fdr_threshold &
      scores[[paste0("orientation_agree_", name)]]
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
    A <- (a - median(a, na.rm = T)) / mad(a, na.rm = T) #scale to generate bins along m
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
    model <- loess(m[b][I] ~ a[b][I], span = sp, degree = dg) #train model on m ~ a (approx. y ~ x)
    expected <- predict(model, a) #predict expected m ~ a
    gi <- m - expected
  }
  result <- list()
  result[["residual"]] <- gi
  result[["predicted"]] <- expected
  return(result)
}
