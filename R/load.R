######
# PRE-PROCESSING CODE
######

#' Get unique gene pairs
#' 
#' This function returns all unique gene pairs from two columns containing gene names in a dataframe.
#' This is necessary for scoring combinatorial data where each row corresponds to a single guide 
#' construct that targets two genes. Intergenic regions must be listed as "---" values, and guides
#' which only target one gene must have the non-targeted gene listed as "None". 
#' 
#' @param df A dataframe where each row corresponds to a single guide construct that targets 
#'   two genes, which contains the columns named in the variables gene_col1 and gene_col2. 
#' @param gene_col1 A column containing the first of two gene names. Intergenic regions must be 
#'   denoted as "---" and non-targeted genes must be denoted as "None". 
#' @param gene_col2 See above.
#' @return A dataframe with two gene name columns where each row contains one unique gene pair. 
#' @export
unique_gene_pairs <- function(df, gene_col1, gene_col2) {
  
  # Gets unique gene pairs
  index <- unlist(apply(df, 1, function(x) sorted_merge(x, gene_col1, gene_col2)))
  index <- index[!duplicated(index)]
  
  # Converts to dataframe
  pair_df <- data.frame(gene1 = rep(NA, length(index)),
                        gene2 = rep(NA, length(index)))
  for (i in 1:nrow(pair_df)) {
    temp <- strsplit(index[i], ";")[[1]]
    pair_df$gene1[i] <- temp[1]
    pair_df$gene2[i] <- temp[2]
  }
  return(pair_df)
}

#' Retrieve guide values for gene pairs
#' 
#' Gets all guides in each orientation for all gene pairs and the given columns, assuming that 
#' guide orientation is determined by the gene name columns. Must be run after \code{unique_gene_pairs}
#' to take the dataframe returned from that function as input. 
#' 
#' @param df A dataframe where each row corresponds to a single guide construct that targets 
#'   two genes, which contains the columns named in the variables gene_col1, gene_col2 and 
#'   guide_cols. 
#' @param gene_pairs Dataframe returned by \code{unique_gene_pairs}.
#' @param gene_col1 A column containing the first of two gene names. Intergenic regions must be 
#'   denoted as "---" and non-targeted genes must be denoted as "None". 
#' @param gene_col2 See above. 
#' @param id_col1 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in gene_col1 (optional, default NULL).  
#' @param id_col2 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in gene_col2 (optional, default NULL).  
#' @param guide_cols A list of all columns containing log-normalized fold-changes.
#' @return A list indexed by unique gene pairs that contains the gene names in the keys 
#'   "gene1" and "gene2". Each column specified in guide_cols has its guide values 
#'   stored in two keys, one for each orientation. For instance, the column "DMSO_T15"
#'   will have its values stored in the keys "orient1_DMSO_T15" and "orient2_DMSO_T15".
#' @export
retrieve_guides_by_gene <- function(df, gene_pairs, gene_col1, gene_col2, guide_cols,
                                    id_col1 = NULL, id_col2 = NULL) {
  
  # Constructs output list of guide-level values for the given columns
  guides <- list()
  for (i in 1:nrow(gene_pairs)) {
    
    # Gets guide values for both orientations
    gene1 <- gene_pairs$gene1[i]
    gene2 <- gene_pairs$gene2[i]
    index1 <- df[[gene_col1]] == gene1 & df[[gene_col2]] == gene2
    index2 <- df[[gene_col1]] == gene2 & df[[gene_col2]] == gene1
    guide_list <- list()
    guide_list[["gene1"]] <- gene1
    guide_list[["gene2"]] <- gene2
    if (sum(index1) >= 0 & sum(index2) > 0) {
      for (col in guide_cols) {
        guide_name1 <- paste0("orient1_", col)
        guide_name2 <- paste0("orient2_", col)
        guide_list[[guide_name1]] <- df[index1, col]
        guide_list[[guide_name2]] <- df[index2, col]
      }
      if (!is.null(id_col1)) {
        guide_list[["orient1_id1"]] <- df[index1, id_col1]
        guide_list[["orient1_id2"]] <- df[index1, id_col2]
      }
      if (!is.null(id_col2)) {
        guide_list[["orient2_id1"]] <- df[index2, id_col1]
        guide_list[["orient2_id2"]] <- df[index2, id_col2]
      }
      guide_list[which(guide_list == numeric())] <- NA
    } else {
      warning(paste("No unique guides found for ", gene1, "and", gene2))
    }
    guides[[i]] <- guide_list
  }
  return(guides)
}

#' Splits guides by type
#'
#' Gets all guides in each orientation for all unique gene pairs and the given columns.
#' Assumes that orientation is determined by label columns which specify whether each
#' guide targets an exonic region or an intergenic region. Exonic-targeting guides must
#' be labeled as "exonic" and intergenic-targeting guides must be labeled as "intergenic".
#' The column specified by label_col1 corresponds to the gene names specified by the 
#' column gene_col1, and similarly for label_col2 and gene_col2. 
#' 
#' 
#' @param df A dataframe where each row corresponds to a single guide construct that targets 
#'   two genes, which contains the columns named in the variables gene_col1, gene_col2, 
#'   label_col1, label_col2 and guide_cols.
#' @param gene_pairs Dataframe returned by \code{unique_gene_pairs}.
#' @param gene_col1 A column containing the first of two gene names. Intergenic regions must be 
#'   denoted as "---" and non-targeted genes must be denoted as "None". 
#' @param gene_col2 See above. 
#' @param label_col1 A column specify whether the guide in gene_col1 targeted an exonic or
#'   intergenic region, which must be labeled as "exonic" or "intergenic". 
#' @param label_col2 See above, but for gene_col2. 
#' @param guide_cols A list of all columns containing log-normalized fold-changes.
#' @param id_col1 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in gene_col1 (optional, default NULL).  
#' @param id_col2 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in gene_col2 (optional, default NULL). 
#' @return A list indexed by unique gene pairs that contains the gene names in the keys 
#'   "gene1" and "gene2". Each column specified in guide_cols has its guide values 
#'   stored in two keys, one for each orientation. For instance, the column "DMSO_T15"
#'   will have its values stored in the keys "orient1_DMSO_T15" and "orient2_DMSO_T15".
#'   The guide type is additionally contained in the key "guide_type". 
#' @export
retrieve_guides_by_label <- function(df, gene_pairs, gene_col1, gene_col2, 
                                     label_col1, label_col2, guide_cols,
                                     id_col1 = NULL, id_col2 = NULL) {
  
  # Constructs output list of guide-level values for the given columns
  guides <- list()
  for (i in 1:nrow(gene_pairs)) {
    
    # Gets guide values for both orientations and all indices
    gene1 <- gene_pairs$gene1[i]
    gene2 <- gene_pairs$gene2[i]
    all_index <- df[[gene_col1]] == gene1 & df[[gene_col2]] == gene2 |
      df[[gene_col1]] == gene2 & df[[gene_col2]] == gene1
    
    # Creates empty guide list to fill depending on guide type
    guide_list <- list()
    guide_list[["gene1"]] <- gene1
    guide_list[["gene2"]] <- gene2
    
    # Handles guides targeting intergenic-intergenic regions
    if (gene1 == "---" & gene2 == "---") {
      guide_list[["guide_type"]] <- "intergenic_intergenic"
      for (col in guide_cols) {
        guide_name <- col
        guide_list[[guide_name]] <- df[all_index, col]
      }
      if (!is.null(id_col1)) {
        guide_list[["id1"]] <- df[all_index, id_col1]
      }
      if (!is.null(id_col2)) {
        guide_list[["id2"]] <- df[all_index, id_col2]
      }
    }
    
    # Handles guides targeting a single gene twice
    else if (gene2 == "None") {
      guide_list[["guide_type"]] <- "single_gene_dual_targeted"
      for (col in guide_cols) {
        guide_name <- col
        guide_list[[guide_name]] <- df[all_index, col]
      }
      if (!is.null(id_col1)) {
        guide_list[["id1"]] <- df[all_index, id_col1]
      }
      if (!is.null(id_col2)) {
        guide_list[["id2"]] <- df[all_index, id_col2]
      }
    }
    
    # Handles exonic-intergenic guides
    else if (gene2 == "---") {
      guide_list[["guide_type"]] <- "exonic_intergenic"
      orient1_index <- all_index & df[[label_col1]] == "exonic"
      orient2_index <- all_index & df[[label_col2]] == "exonic"
      for (col in guide_cols) {
        guide_name1 <- paste0("orient1_", col)
        guide_name2 <- paste0("orient2_", col)
        guide_list[[guide_name1]] <- df[orient1_index, col]
        guide_list[[guide_name2]] <- df[orient2_index, col]
      }
      if (!is.null(id_col1)) {
        guide_list[["orient1_id1"]] <- df[orient1_index, id_col1]
        guide_list[["orient1_id2"]] <- df[orient1_index, id_col2]
      }
      if (!is.null(id_col2)) {
        guide_list[["orient2_id1"]] <- df[orient2_index, id_col1]
        guide_list[["orient2_id2"]] <- df[orient2_index, id_col2]
      }
    }
    
    # Handles exonic-exonic guides targeting two genes
    else {
      guide_list[["guide_type"]] <- "exonic_exonic"
      index1 <- all_index & df[[gene_col1]] == gene1 & df[[gene_col2]] == gene2
      index2 <- all_index & df[[gene_col1]] == gene2 & df[[gene_col2]] == gene1
      for (col in guide_cols) {
        guide_name1 <- paste0("orient1_", col)
        guide_name2 <- paste0("orient2_", col)
        guide_list[[guide_name1]] <- df[index1, col]
        guide_list[[guide_name2]] <- df[index2, col]
      }
      if (!is.null(id_col1)) {
        guide_list[["orient1_id1"]] <- df[index1, id_col1]
        guide_list[["orient1_id2"]] <- df[index1, id_col2]
      }
      if (!is.null(id_col2)) {
        guide_list[["orient2_id1"]] <- df[index2, id_col1]
        guide_list[["orient2_id2"]] <- df[index2, id_col2]
      }
    }
    
    # Appends guides to list and sets empty values to NA
    guide_list[which(guide_list == numeric())] <- NA
    guides[[i]] <- guide_list
  }
  return(guides)
}

#' Splits guides by type.
#' 
#' Splits guide returned from \code{retrieve_guides_by_gene} by targeting type, which is 
#' one of dual-gene exonic-exonic, single-gene exonic-exonic or exonic-intergenic.
#' 
#' @param guides Dataframe returned from \code{retrieve_guides_by_gene}.
#' @return A list of three separate guide lists. Dual-gene exonic-exonic guides are
#'   stored in the key "exonic_exonic", exonic-intergenic guides are stored in the 
#'   key "exonic-intergenic" and single-gene exonic-exonic guides are stored in the
#'   key "single_gene_dual_targeted". 
#' @export
split_guides_by_type <- function(guides) {
  dual_targeted <- guides[unlist(lapply(guides, function(x) x[["guide_type"]] == "single_gene_dual_targeted"))]
  single_targeted <- guides[unlist(lapply(guides, function(x) x[["guide_type"]] == "exonic_intergenic"))]
  pairwise_targeted <- guides[unlist(lapply(guides, function(x) x[["guide_type"]] == "exonic_exonic"))]
  return(list("single_gene_dual_targeted" = dual_targeted, 
              "exonic_intergenic" = single_targeted, 
              "exonic_exonic" = pairwise_targeted))
}

#' Log-normalizes reads.
#' 
#' Log2-normalizes reads with a given pseudocount and scaling factor, and also
#' depth-normalizes the data. 
#' 
#' @param data List of read counts
#' @param cf1 Scaling factor (default 1e6)
#' @param cf2 Pseudocount (default 1)
#' @return Log- and depth-normalized read counts
#' @export
normalize_reads <- function(data, cf1 = 1e6, cf2 = 1) {
  log2((data / sum(data, na.rm = TRUE)) * cf1 + cf2)
}

# Inner function to get sorted merge of two columns for a single row. 
# Maintains a special order for intergenic "---" values and "None"
sorted_merge <- function(row, col1, col2) {
  temp <- sort(c(row[[col1]], row[[col2]]))
  if (temp[1] == "None" | temp[1] == "---") {
    second_item <- temp[2]
    temp[2] <- temp[1]
    temp[1] <- second_item
  }
  temp <- paste0(temp[1], ";", temp[2])
  return(temp)
}