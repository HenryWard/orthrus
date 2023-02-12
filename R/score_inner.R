######
# SCORING HELPERS
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
#'
#' @return A list containing two vectors of orientation-specific null models.
compute_null_model <- function(orient1,
                               orient2,
                               single_gene1, 
                               single_gene2, 
                               combn_vals,
                               combn1,
                               combn2) {
  
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