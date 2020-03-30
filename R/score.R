######
# SCORING CODE
######

#' Scores conditions against a single control.
#' 
#' Scores guides for any number of condition screens against a control screen
#' (e.g. for directly comparing drug response to DMSO response). Do NOT use 
#' this function if you would like to compare guides against a null model 
#' derived from single-gene effects. In that case, call \code{score_dual_vs_single}
#' instead. After running this function, pass the resulting dataframe to 
#' \code{call_significant_response} to call significant effects.
#' 
#' @param guides A list of guides returned from \code{retrieve_guides_by_label} or 
#'   from an entry of \code{split_guides_by_type}.
#' @param control_cols A list where the first entry is a desired name for a condition, 
#'   e.g. DMSO, and all subsequent values are column names of the technical replicates 
#'   for that condition. All technical replicates specified must have been passed to
#'   \code{retrieve_guides_by_gene} or \code{retrieve_guides_by_gene} in the guide_cols
#'   argument. For example, c("DMSO", "DMSO_T15A", "DMSO_T15B", "DMSO_T15C"). 
#' @param ... Any number of conditions to score against the specified control, passed 
#'   in as separate lists formatted the same as control_cols. 
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @return A dataframe of scored data with separate columns given by the specified control
#'   and condition names.
#' @export
score_conditions_vs_control <- function(guides, control_cols, ..., 
                                        min_guides = 3) {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_cols[1]
  control_cols <- control_cols[2:length(control_cols)]
  condition_names <- c()
  condition_cols <- list()
  for (condition in list(...)) {
    name <- condition[1]
    condition_names <- c(condition_names, name)
    condition_cols[[name]] <- condition[2:length(condition)]
  }
  
  # Makes output dataframes
  scores <- data.frame(gene1 = rep(NA, length(guides)), gene2 = NA)
  condition_residuals <- list()
  for (name in condition_names) {
    condition_residuals[[name]] <- data.frame(residual1 = rep(NA, length(guides)), 
                                              residual2 = NA, residual3 = NA,
                                              residual4 = NA, residual5 = NA, 
                                              residual6 = NA)
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
  for (i in 1:length(guides)) {
    
    # Gets gene names and control guide values across replicates
    guide_vals <- guides[[i]]
    scores$gene1[i] <- guide_vals$gene1
    scores$gene2[i] <- guide_vals$gene2
    rep_mean_control <- rowMeans(data.frame(guide_vals[control_cols]))
    
    # Takes the mean across replicates for all conditions
    for (name in condition_names) {
      rep_mean_condition <- rowMeans(data.frame(guide_vals[condition_cols[[name]]]))
      diff <- rep_mean_condition - rep_mean_control
      scores[[paste0("n_", control_name)]][i] <- length(rep_mean_control)
      scores[[paste0("n_", name)]][i] <- length(rep_mean_condition)
      scores[[paste0("mean_", control_name)]][i] <- mean(rep_mean_control)
      scores[[paste0("mean_", name)]][i] <- mean(rep_mean_condition)
      scores[[paste0("differential_", name, "_vs_", control_name)]][i] <- mean(diff)
      
      # Stores residuals for moderated t-test
      if (length(diff) < 6) { diff <- c(diff, rep(NA, 6 - length(diff))) } 
      condition_residuals[[name]][i,1:6] <- diff 
    }
  }
  
  # Scores condition response with moderated t-test
  for (name in condition_names) {
    ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]]))
    p_val <- ebayes_fit$p.value[,1]
    scores[[paste0("pval_", name, "_vs_", control_name)]] <- p_val
    scores[[paste0("fdr_", name, "_vs_", control_name)]] <- 
      p.adjust(scores[[paste0("pval_", name, "_vs_", control_name)]], method = "BH")
  }
  
  # Removes genes with too few observations
  scores <- scores[scores[[paste0("n_", control_name)]] >= min_guides,]
  
  # Explicitly returns scoresd data
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
#' @param dual_guides A list of exonic-exonic guides returned from \code{retrieve_guides_by_label} 
#'   or from the exonic-exonic entry of \code{split_guides_by_type}.
#' @param single_guides A list of exonic-intergenic guides returned from \code{retrieve_guides_by_label} 
#'   or from the exonic-intergenic entry of \code{split_guides_by_type}.
#' @param ... Any number of lists to score against a multiplicative null model derived 
#'   from single-gene effects. Each list must be formatted where the first entry is a 
#'   desired name for a condition, e.g. DMSO, and all subsequent values are column names 
#'   of the technical replicates for that condition. All technical replicates specified 
#'   must have been passed to \code{retrieve_guides_by_gene} or 
#'   \code{retrieve_guides_by_label} in the guide_cols argument. For example, 
#'   \code{c("DMSO", "DMSO_T15A", "DMSO_T15B", "DMSO_T15C")}. 
#' @param test The desired type of hypothesis testing to run. Must be one of "rank-sum" or
#'   "moderated-t".
#' @param min_guides The minimum number of guides per gene pair required to score data 
#'   (default 3).
#' @return A dataframe of scored data with separate columns given by the specified control
#'   and condition names.
#' @export
score_dual_vs_single <- function(dual_guides, single_guides, ..., 
                                 test = "rank-sum", min_guides = 3) {
  
  # Gets condition names and columns for any number of conditions
  condition_names <- c()
  condition_cols <- list()
  for (condition in list(...)) {
    name <- condition[1]
    condition_names <- c(condition_names, name)
    orient1 <- paste0("orient1_", condition[2:length(condition)])
    orient2 <- paste0("orient2_", condition[2:length(condition)])
    condition_cols[[name]][["orient1"]] <- orient1
    condition_cols[[name]][["orient2"]] <- orient2
  }
  
  # Makes output dataframe
  scores <- data.frame(gene1 = rep(NA, length(dual_guides)), gene2 = NA)
  
  # Appends additional columns for each condition
  new_cols <- c()
  for (name in condition_names) {
    new_cols <- c(new_cols, c(
      paste0("n_dual_", name), 
      paste0("n_single_", name),
      paste0("mean_dual_", name),
      paste0("mean_single_", name),
      paste0("orientation_agree_", name),
      paste0("differential_dual_vs_single_", name),
      paste0("pval1_dual_vs_single_", name),
      paste0("pval2_dual_vs_single_", name),
      paste0("fdr1_dual_vs_single_", name),
      paste0("fdr2_dual_vs_single_", name),
      paste0("significant_dual_vs_single_", name)
    ))
  }
  scores[new_cols] <- NA
  
  # Makes residual dataframes if necessary
  max_guides <- -1
  condition_residuals <- list()
  if (test == "moderated-t") {
    
    # Gets max number of guides first
    for (guide in dual_guides) {
      guide_names <- names(guide)[!(names(guide) %in% c("gene1", "gene2", "guide_type"))]
      guide_num <- max(unlist(lapply(guide_names, function(x) length(guide[[x]]))))
      max_guides <- max(max_guides, guide_num)
    }
    
    # Makes residual dataframes with columns equal to the max number of guides
    for (name in condition_names) {
      residual_df <- data.frame(matrix(nrow = length(dual_guides), ncol = max_guides))
      colnames(residual_df) <- paste0("guide_residual_", 1:max_guides)
      condition_residuals[[name]][[1]] <- residual_df
      condition_residuals[[name]][[2]] <- residual_df
    }
  }
  
  # Scores guides for each condition
  for (i in 1:length(dual_guides)) {
    
    # Gets gene names and guide values
    dual_vals <- dual_guides[[i]]
    gene1 <- dual_vals$gene1
    gene2 <- dual_vals$gene2
    scores$gene1[i] <- gene1
    scores$gene2[i] <- gene2
    
    # Finds matching single-targeting guides
    single_gene1 <- single_guides[unlist(lapply(single_guides, function(x) x[["gene1"]] == gene1))][[1]]
    single_gene2 <- single_guides[unlist(lapply(single_guides, function(x) x[["gene1"]] == gene2))][[1]]
    
    # scoress dual-targeting guides vs. single-targeting null model
    for (name in condition_names) {
      
      # Gets column names for each orientation and the current condition
      orient1 <- condition_cols[[name]][["orient1"]]
      orient2 <- condition_cols[[name]][["orient2"]]
      
      # Gets guide values for single-targeting guides
      single_gene1_orient1 <- rowMeans(data.frame(single_gene1[orient1]))
      single_gene1_orient2 <-  rowMeans(data.frame(single_gene1[orient2]))
      single_gene2_orient1 <- rowMeans(data.frame(single_gene2[orient1]))
      single_gene2_orient2 <- rowMeans(data.frame(single_gene2[orient2]))
      
      # Gets null model by summing different orientations of single-targeting guides
      null1 <- unlist(apply(expand.grid(single_gene1_orient1, single_gene2_orient2), 1, sum))
      null2 <- unlist(apply(expand.grid(single_gene1_orient2, single_gene2_orient1), 1, sum))
      
      # Takes the dual-targeting mean across replicates
      dual1 <- rowMeans(data.frame(dual_vals[orient1]))
      dual2 <- rowMeans(data.frame(dual_vals[orient2]))
      
      # Gets residuals
      scores[[paste0("orientation_agree_", name)]] <- 
        (sign(mean(dual1) - mean(null1)) == sign(mean(dual2) - mean(null2)))
      scores[[paste0("differential_dual_vs_single_", name)]][i] <- mean(c(dual1, dual2)) - mean(c(null1, null2))
      
      # Performs the specified type of testing or stores residuals for later testing
      if (test == "rank-sum") {
        scores[[paste0("pval1_dual_vs_single_", name)]][i] <- suppressWarnings(wilcox.test(dual1, null1))$p.value
        scores[[paste0("pval2_dual_vs_single_", name)]][i] <- suppressWarnings(wilcox.test(dual2, null2))$p.value 
      } else if (test == "moderated-t") {
        if (length(dual1) != length(null1)) {
          print(paste("Length of dual1: ", length(dual1)))
          print(paste("Length of null1: ", length(null1)))
          print(paste("Gene name 1: ", gene1))
          print(paste("Gene name 2: ", gene2))
        }
        residuals1 <- dual1 - null1
        residuals2 <- dual2 - null2
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
      scores[[paste0("n_dual_", name)]][i] <- length(dual1)
      scores[[paste0("n_single_", name)]][i] <- length(null1)
      scores[[paste0("mean_dual_", name)]][i] <- mean(c(dual1, dual2))
      scores[[paste0("mean_single_", name)]][i] <- mean(c(null1, null2))
    }
  }
  
  # Scores condition response with moderated t-test if specified
  if (test == "moderated-t") {
    for (name in condition_names) {
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]][[1]]))
      p_val1 <- ebayes_fit$p.value[,1]
      ebayes_fit <- limma::eBayes(limma::lmFit(condition_residuals[[name]][[2]]))
      p_val2 <- ebayes_fit$p.value[,1]
      scores[[paste0("pval1_dual_vs_single_", name)]][i] <- p_val1
      scores[[paste0("pval2_dual_vs_single_", name)]][i] <- p_val2
    }
  }
  
  # Gets FDRs
  for (name in condition_names) {
    scores[[paste0("fdr1_dual_vs_single_", name)]] <- 
      p.adjust(scores[[paste0("pval1_dual_vs_single_", name)]], method = "BH")
    scores[[paste0("fdr2_dual_vs_single_", name)]] <- 
      p.adjust(scores[[paste0("pval2_dual_vs_single_", name)]], method = "BH")
  }
  
  # Removes genes with too few observations
  # to_remove <- rep(FALSE, nrow(scores))
  # for (name in condition_names) {
  #   to_remove <- to_remove | scores[[paste0("n_dual", name)]] < min_guides
  #   to_remove <- to_remove | scores[[paste0("n_single", name)]] < min_guides
  # }
  # scores <- scores[!to_remove,]
  
  # Explicitly returns scoresd data
  return(scores)
}

#' Call significant responses for scored data.
#' 
#' Run this to call significant responses for data returned from 
#' \code{score_conditions_vs_control}. Do NOT run this on data returned from
#' \code{score_dual_vs_single}.
#' 
#' @param scores Dataframe returned from \code{score_conditions_vs_control}.
#' @param control_cols A list where the first entry is a desired name for a condition, 
#'   e.g. DMSO, and all subsequent values are column names of the technical replicates 
#'   for that condition. See documentation for \code{score_conditions_vs_control}. 
#' @param ... Any number of conditions to score against the specified control, passed 
#'   in as separate lists formatted the same as control_cols.
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
call_significant_response <- function(scores, control_cols, ...,
                                      fdr_threshold = 0.1, differential_threshold = 0,
                                      neg_type = "Negative", pos_type = "Positive") {
  
  # Gets condition names and columns for any number of conditions
  control_name <- control_cols[1]
  control_cols <- control_cols[2:length(control_cols)]
  condition_names <- c()
  condition_cols <- list()
  for (condition in list(...)) {
    name <- condition[1]
    condition_names <- c(condition_names, name)
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
#' \code{score_dual_vs_single}. Do NOT run this on data returned from
#' \code{score_conditions_vs_control}.
#' 
#' @param scores Dataframe returned from \code{score_dual_vs_single}.
#' @param control_cols A list where the first entry is a desired name for a condition, 
#'   e.g. DMSO, and all subsequent values are column names of the technical replicates 
#'   for that condition. See documentation for \code{score_dual_vs_single}. 
#' @param ... Any number of lists formatted where the first entry is a desired name for a 
#'   condition, e.g. DMSO, and all subsequent values are column names of the technical 
#'   replicates for that condition. All technical replicates specified must have been passed 
#'   to \code{retrieve_guides_by_gene} or \code{retrieve_guides_by_label} 
#'   in the guide_cols argument. For example, 
#'   \code{c("DMSO", "DMSO_T15A", "DMSO_T15B", "DMSO_T15C")}. 
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
call_significant_response_dual <- function(scores, ...,
                                           fdr_threshold = 0.1, differential_threshold = 0,
                                           neg_type = "Negative", pos_type = "Positive") {
  
  # Gets condition names and columns for any number of conditions
  condition_names <- c()
  condition_cols <- list()
  for (condition in list(...)) {
    name <- condition[1]
    condition_names <- c(condition_names, name)
  }
  
  # Calls significant differences for dual-targeting vs. single-targeting null model
  for (name in condition_names) {
    scores[[paste0("significant_dual_vs_single_", name)]] <- 
      scores[[paste0("fdr1_dual_vs_single_", name)]] < fdr_threshold &
      scores[[paste0("fdr2_dual_vs_single_", name)]] < fdr_threshold &
      scores[[paste0("orientation_agree_", name)]]
  }
  
  # Makes thresholded calls for significant negative and positive effects
  for (name in condition_names) {
    response_col <- paste0("effect_type_", name)
    scores[[response_col]] <- "None"
    diffs <- scores[[paste0("differential_dual_vs_single_", name)]]
    sig <- scores[[paste0("significant_dual_vs_single_", name)]]
    scores[[response_col]][sig & diffs < 0 & abs(diffs) > differential_threshold] <- neg_type
    scores[[response_col]][sig & diffs > 0 & abs(diffs) > differential_threshold] <- pos_type
  }
  
  # Explicitly returns scored data
  return(scores)
}

# Fits a loess curve to predict y given x
gi_MA_loess_Func <- function(x, y, sp = 0.4, dg = 2, binSize = 100) {
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
    model <- loess(m[b][I] ~ a[b][I], span=sp, degree = dg) #train model on m ~ a (approx. y ~ x)
    expected <- predict(model, a) #predict expected m ~ a
    gi <- m - expected
  }
  return(gi)
}
