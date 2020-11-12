######
# PRE-PROCESSING CODE
######

#' Splits data into different types of guides
#' 
#' 
#' @param guides Dataframe returned from \code{retrieve_guides_by_gene}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param id_col1 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in the gene1 column (optional, default NULL).  
#' @param id_col2 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in the gene1 column (optional, default NULL).  
#' @return A list of three separate guide lists. Dual-gene exonic-exonic guides are
#'   stored in the key "exonic_exonic", exonic-intergenic guides are stored in the 
#'   key "exonic_intergenic" and single-gene exonic-exonic guides are stored in the
#'   key "single_gene_dual_targeted". 
split_guides <- function(guides, screens, id_col1 = NULL, id_col2 = NULL) {
  gene_pairs <- unique_gene_pairs(guides)
  guides <- retrieve_guides_by_gene(guides, gene_pairs, screens, id_col1, id_col2)
  return(split_guides_by_type(guides))
}

#' Get unique gene pairs
#' 
#' This function returns all unique gene pairs from two columns containing gene names in a dataframe.
#' This is necessary for scoring combinatorial data where each row corresponds to a single guide 
#' construct that targets two genes. NegControl regions must be listed as "NegControl" values, and guides
#' which only target one gene must have the non-targeted gene listed as "None". 
#' 
#' @param df A dataframe where each row corresponds to a single guide construct that targets 
#'   two genes, which contains the columns gene1 and gene2. 
#' @return A dataframe with two gene name columns where each row contains one unique gene pair. 
unique_gene_pairs <- function(df) {
  
  # Gets unique gene pairs
  index <- unlist(apply(df, 1, function(x) sorted_merge(x, "gene1", "gene2")))
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
#'   two genes, which contains the columns gene1, gene2 and optionally the columns specified by
#'   id_col1 and id_col2.
#' @param gene_pairs Dataframe returned by \code{unique_gene_pairs}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param id_col1 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in the column gene1 (optional, default NULL).  
#' @param id_col2 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in the column gene2 (optional, default NULL).  
#' @return A list indexed by unique gene pairs that contains the gene names in the keys 
#'   "gene1" and "gene2". Each column specified in guide_cols has its guide values 
#'   stored in two keys, one for each orientation. For instance, the column "DMSO_T15"
#'   will have its values stored in the keys "orient1_DMSO_T15" and "orient2_DMSO_T15".
retrieve_guides_by_gene <- function(df, gene_pairs, screens,
                                    id_col1 = NULL, id_col2 = NULL) {
  
  # Gets all columns
  guide_cols <- c()
  for (screen in screens) {
    guide_cols <- c(guide_cols, screen[["replicates"]]) 
  }
  
  # Constructs output list of guide-level values for the given columns
  guides <- list()
  for (i in 1:nrow(gene_pairs)) {
    
    # Gets guide values for both orientations
    gene1 <- gene_pairs$gene1[i]
    gene2 <- gene_pairs$gene2[i]
    all_index <- df$gene1 == gene1 & df$gene2 == gene2 | 
      df$gene1 == gene2 & df$gene2 == gene1
    guide_list <- list()
    guide_list[["gene1"]] <- gene1
    guide_list[["gene2"]] <- gene2
    
    # Handles guides targeting intergenic-intergenic regions and 
    # guides that target a single gene twice
    if ((gene1 == "NegControl" & gene2 == "NegControl") | gene2 == "None") {
      if (gene1 == "NegControl" & gene2 == "NegControl") {
        guide_list[["guide_type"]] <- "intergenic_intergenic"
      } else if (gene2 == "None") {
        guide_list[["guide_type"]] <- "single_gene_dual_targeted"
      }
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
    else if (gene1 == "NegControl" | gene2 == "NegControl") {
      guide_list[["guide_type"]] <- "exonic_intergenic"
      orient1_index <- all_index & df$gene2 == "NegControl"
      orient2_index <- all_index & df$gene1 == "NegControl"
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
    else if (gene1 != gene2) {
      guide_list[["guide_type"]] <- "exonic_exonic"
      index1 <- all_index & df$gene1 == gene1 & df$gene2 == gene2
      index2 <- all_index & df$gene1 == gene2 & df$gene2 == gene1
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
    
    # Throws a warning if no guides found for the current gene pair
    else {
      warning(paste("No unique guides found for", gene1, "and", gene2))
    }
    
    # Appends guides to list and sets empty values to NA
    guide_list[which(guide_list == numeric())] <- NA
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
#'   label_col1, label_col2 and replicate columns in screens.
#' @param gene_pairs Dataframe returned by \code{unique_gene_pairs}.
#' @param screens List of screens generated with \code{add_screens}.
#' @param label_col1 A column specify whether the guide in gene_col1 targeted an exonic or
#'   intergenic region, which must be labeled as "exonic" or "intergenic". 
#' @param label_col2 See above, but for gene_col2. 
#' @param id_col1 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in gene_col1 (optional, default NULL).  
#' @param id_col2 A column containing guide-specific indices, e.g. guide sequences or IDs, for
#'   guides in gene_col2 (optional, default NULL). 
#' @return A list indexed by unique gene pairs that contains the gene names in the keys 
#'   "gene1" and "gene2". Each column contained in each screen's replicates has its guide values 
#'   stored in two keys, one for each orientation. For instance, the column "DMSO_T15"
#'   will have its values stored in the keys "orient1_DMSO_T15" and "orient2_DMSO_T15".
#'   The guide type is additionally contained in the key "guide_type". 
retrieve_guides_by_label <- function(df, gene_pairs, screens,
                                     label_col1, label_col2,
                                     id_col1 = NULL, id_col2 = NULL) {
  
  # Gets all columns
  guide_cols <- c()
  for (screen in screens) {
   guide_cols <- c(guide_cols, screen[["replicates"]]) 
  }
  
  # Constructs output list of guide-level values for the given columns
  guides <- list()
  for (i in 1:nrow(gene_pairs)) {
    
    # Gets guide values for both orientations and all indices
    gene1 <- gene_pairs$gene1[i]
    gene2 <- gene_pairs$gene2[i]
    all_index <- df$gene1 == gene1 & df$gene2 == gene2 |
      df$gene1 == gene2 & df$gene2 == gene1
    
    # Creates empty guide list to fill depending on guide type
    guide_list <- list()
    guide_list[["gene1"]] <- gene1
    guide_list[["gene2"]] <- gene2
    
    # Handles guides targeting intergenic-intergenic regions
    if (gene1 == "NegControl" & gene2 == "NegControl") {
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
    else if (gene2 == "NegControl") {
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
      index1 <- all_index & df$gene1 == gene1 & df$gene2 == gene2
      index2 <- all_index & df$gene1 == gene2 & df$gene2 == gene1
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
#'   key "exonic_intergenic" and single-gene exonic-exonic guides are stored in the
#'   key "single_gene_dual_targeted". 
split_guides_by_type <- function(guides) {
  dual_targeted <- guides[unlist(lapply(guides, function(x) x[["guide_type"]] == "single_gene_dual_targeted"))]
  single_targeted <- guides[unlist(lapply(guides, function(x) x[["guide_type"]] == "exonic_intergenic"))]
  pairwise_targeted <- guides[unlist(lapply(guides, function(x) x[["guide_type"]] == "exonic_exonic"))]
  return(list("single_gene_dual_targeted" = dual_targeted, 
              "exonic_intergenic" = single_targeted, 
              "exonic_exonic" = pairwise_targeted))
}

#' Normalizes reads for given screens
#' 
#' Log2 and depth-normalizes reads between a given list of columns and a given
#' column of an earlier timepoint (e.g. T0) specified in the "normalize_name"
#' entry of each screen. Screens with NULL for their "normalize_name" entry
#' are log2 and depth-normalized, but not normalized to earlier timepoints.
#' If a screen to normalize against has multiple replicates, those replicates
#' are averaged before normalization. Multiple replicates for a screen being
#' normalized, however, are normalized separately against the provided early
#' timepoint.
#' 
#' @param df Reads dataframe.
#' @param screens List of screens generated with \code{add_screens}. 
#' @param filter_names List of screen names to filter based on read counts by. 
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @param min_reads Minimum number of reads to keep (default 30, anything
#'   below this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (default 10000, anything
#'   above this value will be filtered out).
#' @return Normalized dataframe.
#' @export 
normalize_screens <- function(df, screens, filter_names = NULL, cf1 = 1e6, cf2 = 1, 
                              min_reads = 30, max_reads = 10000) {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Flags guides with too few read counts
  all_names <- names(screens)
  for (name in filter_names) {
    if (!(name %in% all_names)) {
      cat(paste("WARNING: screen", name, "not found, data not filtered by this screen\n"))
    }
  }
  to_remove <- rep(FALSE, nrow(df))
  filter_cols <- sapply(screens[filter_names], "[[", "replicates")
  for (col in filter_cols) {
    to_remove[df[,col] < min_reads] <- TRUE
  }
  sum_low <- sum(to_remove)
  for (col in filter_cols) {
    to_remove[df[,col] > max_reads] <- TRUE
  }
  sum_high <- sum(to_remove) - sum_low
  removed_guides_ind <- which(to_remove)
  
  # Log2 and depth-normalizes every screen
  for (screen in screens) {
    for (col in screen[["replicates"]]) {
      df[,col] <- normalize_reads(df[,col], cf1, cf2)
    }
  }
  
  # Normalizes specified screens to earlier timepoints
  for (screen in screens) {
    normalize_name <- screen[["normalize_name"]]
    if (!is.null(normalize_name)) {
      if (normalize_name %in% all_names) {
        for (col in screen[["replicates"]]) {
          rep_cols <- screens[[normalize_name]][["replicates"]]
          rep_norm <- df[,rep_cols]
          if (length(rep_cols) > 1) {
            rep_norm <- rowMeans(rep_norm)
          }
          df[,col] <- df[,col] - rep_norm
        }
      } else {
        cat(paste("WARNING: screen", normalize_name, "not found.\n"))
      }
    }
  }
  
  # Removes flagged guides
  df <- df[!to_remove,]
  cat(paste("Excluded a total of", sum_low, "guides for low t0 representation\n"))
  cat(paste("Excluded a total of", sum_high, "guides for high t0 representation\n"))
  return(df)
}

#' Filters guides with too few read counts.
#' 
#' Filters guides out with too few read counts from a given reads dataframe
#' and a given set of columns (e.g. T0 columns).
#' 
#' @param df Reads dataframe.
#' @param cols Columns to filter by.
#' @param min_reads Minimum number of reads to keep (anything below
#'   this value will be filtered out).
#' @param max_reads Maximum number of reads to keep (anything above
#'   this value will be filtered out).
#' @return Filtered dataframe.
filter_reads <- function(df, cols, min_reads = 30, max_reads = 10000) {
  to_remove <- rep(FALSE, nrow(df))
  for (col in cols) {
    to_remove[df[,col] < min_reads] <- TRUE
  }
  sum_low <- sum(to_remove)
  for (col in cols) {
    to_remove[df[,col] > max_reads] <- TRUE
  }
  sum_high <- sum(to_remove) - sum_low
  removed_guides_ind <- which(to_remove)
  df <- df[!to_remove,]
  cat(paste("Excluded a total of", sum_low, "guides for low t0 representation\n"))
  cat(paste("Excluded a total of", sum_high, "guides for low t0 representation\n"))
  return(df)
}

#' Log-normalizes reads.
#' 
#' Log2-normalizes reads with a given pseudocount and scaling factor, and also
#' depth-normalizes the data. 
#' 
#' @param df List of read counts.
#' @param cf1 Scaling factor (default 1e6).
#' @param cf2 Pseudocount (default 1).
#' @return Log- and depth-normalized read counts.
#' @export
normalize_reads <- function(df, cf1 = 1e6, cf2 = 1) {
  log2((df / sum(df, na.rm = TRUE)) * cf1 + cf2)
}

# Inner function to get sorted merge of two columns for a single row. 
# Maintains a special order for intergenic "NegControl" values and "None"
sorted_merge <- function(row, col1, col2) {
  temp <- sort(c(row[[col1]], row[[col2]]))
  if (temp[1] == "None" | temp[1] == "NegControl") {
    second_item <- temp[2]
    temp[2] <- temp[1]
    temp[1] <- second_item
  }
  temp <- paste0(temp[1], ";", temp[2])
  return(temp)
}