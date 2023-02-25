######
# SCORING HELPER FUNCTIONS
######

#' Gets null model for each orientation
#' 
#' When guide IDs are not specified, all possible combinations of gene1/gene2 
#' and gene2/gene1 are computed as null values. When guide IDs are specified, 
#' single gene LFCs with matching IDs are instead used to compute gene1/gene2
#' and gene2/gene1 null values. In this case, if there are multiple single gene
#' LFCs matching to a single combinatorial ID, they are averaged before 
#' computing null values.
#' 
#' @param orient1 Column name of first orientation
#' @param orient2 Colum name of second orientation
#' @param single_gene1 One element of exonic-intergenic guides returned from 
#'   \code{split_guides}.
#' @param single_gene2 One element of exonic-intergenic guides returned from 
#'   \code{split_guides}.
#' @param combn_vals One element of exonic-exonic guides returned from 
#'   \code{split_guides}.
#' @param combn1 Replicate-averaged LFCs for the first orientation of combn_vals.
#' @param combn2 Replicate-averaged LFCs for the second orientation of combn_vals.
#' @param collapse_single_targeting If TRUE, takes the mean of single-targeting 
#'   controls when there are multiple controls that match a given gene pair.
#' @param ignore_orientation If TRUE, aggregates guides across both orientations, 
#'   returning only one p-value and FDR column with orientation2 p-values set to NA.
#' @param subset_to_matching_guide_id If TRUE, filters single and combinatorial guides to 
#'   matching guide IDs in both (default TRUE). 
#'
#' @return A list containing two vectors of orientation-specific null models.
compute_null_model <- function(orient1,
                               orient2,
                               single_gene1, 
                               single_gene2, 
                               combn_vals,
                               combn1,
                               combn2,
                               collapse_single_targeting,
                               ignore_orientation,
                               subset_to_matching_guide_id) {
  
  # Gets guide values for single-targeting guides
  single_gene1_orient1 <- rowMeans(data.frame(single_gene1[orient1]))
  single_gene1_orient2 <- rowMeans(data.frame(single_gene1[orient2]))
  single_gene2_orient1 <- rowMeans(data.frame(single_gene2[orient1]))
  single_gene2_orient2 <- rowMeans(data.frame(single_gene2[orient2]))
  
  # Gets null model by summing different orientations of single-targeting guides
  null1 <- c()
  null2 <- c()
  if (!subset_to_matching_guide_id) {
    null1 <- unlist(apply(expand.grid(single_gene1_orient1, single_gene2_orient2), 1, sum))
    null2 <- unlist(apply(expand.grid(single_gene1_orient2, single_gene2_orient1), 1, sum))
  } else {
    for (j in 1:length(combn1)) {
      id1 <- combn_vals[["orient1_id1"]][j]
      id2 <- combn_vals[["orient1_id2"]][j]
      val1 <- single_gene1_orient1[single_gene1[["orient1_id1"]] == id1]
      val2 <- single_gene2_orient2[single_gene2[["orient2_id2"]] == id2]
      if ((collapse_single_targeting == TRUE) & (ignore_orientation == TRUE)) {
        if (length(val1) > 1) {
          val1 <- mean(val1) 
        }
        if (length(val2) > 1) { 
          val2 <- mean(val2) 
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
      if ((collapse_single_targeting == TRUE) & (ignore_orientation == TRUE)) {
        if (length(val1) > 1) { 
          val1 <- mean(val1) 
        }
        if (length(val2) > 1) { 
          val2 <- mean(val2) 
        }
      }
      if (length(val1) > 0 & length(val2) > 0) {
        null2 <- c(null2, val1 + val2)
      } else {
        null2 <- c(null2, NA)
      }
    }
  }
  
  # Explicitly returns orientation-specific null models
  return(list(null1, null2))
}

#' Subsets combn and single-targeting guides to matching IDs
#' 
#' @param combn_vals One element of exonic-exonic guides returned from 
#'   \code{split_guides}.
#' @param single_gene1 One element of exonic-intergenic guides returned from 
#'   \code{split_guides}.
#' @param single_gene2 One element of exonic-intergenic guides returned from 
#'   \code{split_guides}.
#'   
#' @return A list containing three lists of matched combn_vals, single_gene1, and single_2
subset_to_matching_ids <- function(combn_vals, 
                                   single_gene1,
                                   single_gene2) {
  if ("orient1_id1" %in% names(combn_vals) & "orient1_id2" %in% names(combn_vals)) {
    single_gene1 <- single_gene1[single_gene1[["orient1_id1"]] %in% combn_vals[["orient1_id1"]] |
                                   single_gene1[["orient2_id2"]] %in% combn_vals[["orient2_id2"]]]
    single_gene2 <- single_gene2[single_gene2[["orient1_id1"]] %in% combn_vals[["orient2_id1"]] |
                                   single_gene2[["orient2_id2"]] %in% combn_vals[["orient1_id2"]]]
    combn_vals <- combn_vals[combn_vals[["orient1_id1"]] %in% single_gene1[["orient1_id1"]] &
                               combn_vals[["orient1_id2"]] %in% single_gene2[["orient2_id2"]]]
    combn_vals <- combn_vals[combn_vals[["orient2_id1"]] %in% single_gene2[["orient1_id1"]] &
                               combn_vals[["orient2_id2"]] %in% single_gene1[["orient2_id2"]]]
  }
  return(list(combn_vals, single_gene1, single_gene2))
}

#' Reproducibly sets residuals to a given maximum number
#' 
#' @param resid_grid A numeric list of per-guide condition and control values
#'   created with expand.grid, where column 1 is conditions and 2 is controls
#' @param max_resid Maximum number of residuals to keep
#'   
#' @return If length(resid_grid) > max_resid, subsamples resid_grid to max_resid. 
#'   Otherwise, adds NA values to resid_grid so that length(resid_grid) == 
#'   max_resid
fix_residual_length <- function(resid_grid, max_resid) {
  if (nrow(resid_grid) > max_resid) {
    prev_seed <- .Random.seed
    set.seed(123456)
    resid_grid <- resid_grid[sample(nrow(resid_grid), max_resid),]
    .Random.seed <- prev_seed
  } else if (nrow(resid_grid) < max_resid) {
    additional_rows <- max_resid - nrow(resid_grid)
    to_append <- data.frame(matrix(NA, nrow = additional_rows, ncol = 2))
    colnames(to_append) <- colnames(resid_grid)
    resid_grid <- rbind(resid_grid, to_append) 
  }
  return(resid_grid)
}