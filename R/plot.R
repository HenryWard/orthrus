######
# PLOTTING CODE
######

#' Plot replicate comparisons.
#' 
#' Plots replicate comparisons for all replicates in a list of screens and outputs
#' plots to a given folder.
#' 
#' @param df Reads or lfc dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @param output_folder Folder to output plots to. 
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export 
plot_lfc_qc <- function(df, screens, output_folder, plot_type = "png") {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  # Builds dataframe of replicate PCCs
  pcc_df <- NULL
  
  # Compares replicates across all screens
  for (screen in screens) {
    rep_cols <- screen[["replicates"]]
    if (length(rep_cols) > 1) {
      pairs <- utils::combn(rep_cols, 2)
      for (i in 1:ncol(pairs)) {
        col1 <- pairs[1,i]
        col2 <- pairs[2,i]
        x_label <- paste0(col1, " log fold change")
        y_label <- paste0(col2, " log fold change")
        temp <- plot_samples(df, col1, col2, x_label, y_label, print_cor = TRUE)
        p <- temp[[1]]
        pcc <- temp[[2]]
        file_name <- paste0(col1, "_vs_", col2, "_replicate_comparison.", plot_type)
        ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)
        
        # Stores PCC in dataframe
        if (is.null(pcc_df)) {
          pcc_df <- data.frame(rep1 = col1, rep2 = col2, pcc = pcc,
                               stringsAsFactors = FALSE)
        } else {
          pcc_df <- rbind(pcc_df, c(col1, col2, pcc))
        }
      }   
    }
  }
  
  # Writes PCCs to file
  pcc_file <- file.path(output_folder, "replicate_pcc.tsv")
  utils::write.table(pcc_df, pcc_file, quote = FALSE, sep = "\t",
                     row.names = FALSE, col.names = TRUE)
}

#' Plot read counts for a screen.
#' 
#' Plots a histogram of read counts for each replicate of all screens. Also
#' plots total reads for all screens.
#' 
#' @param df Reads dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @param output_folder Folder to output plots to. 
#' @param log_scale If true, log-normalizes data.
#' @param pseudocount Pseudocounts to add to log-normalized data if specified (default 1).
#' @param display_numbers Whether or not to include PCC values in heatmap (default TRUE).
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export
plot_reads_qc <- function(df, screens, output_folder,
                          log_scale = TRUE, pseudocount = 1,
                          display_numbers = TRUE, 
                          plot_type = "png") {
  
  # Checks for input errors
  check_screen_params(df, screens)
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  # Plots read count histograms for all replicates of all screens and stores total reads 
  reads_df <- NULL
  all_cols <- c()
  col_groups <- c()
  all_coverage <- c()
  i <- 1
  for (screen_name in names(screens)) {
    screen <- screens[[screen_name]]
    for (col in screen[["replicates"]]) {
      total_reads <- sum(df[,col], na.rm = TRUE)
      if (is.null(reads_df)) {
        reads_df <- data.frame(rep = col, reads = total_reads,
                               stringsAsFactors = FALSE)
      } else {
        reads_df <- rbind(reads_df, c(col, total_reads))
      }
      all_coverage <- c(all_coverage, screen[["target_coverage"]])
      all_cols <- c(all_cols, col)
      p <- plot_reads(df, col, log_scale, pseudocount)
      file_name <- paste0(col, "_raw_reads_histogram.", plot_type)
      ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)  
      col_groups[i] <- screen_name
      i <- i + 1
    }
  }
  
  # Gets unique coverage breakpoints
  all_coverage <- unique(all_coverage)
  all_coverage <- all_coverage*nrow(df)
  
  # Plots total reads
  reads_df$reads <- as.numeric(reads_df$reads)
  p <- ggplot2::ggplot(reads_df, ggplot2::aes_string(x = "rep", y = "reads")) +
    ggplot2::geom_bar(stat = "identity", color = "Black", fill = "gray30")
  for (coverage in all_coverage) {
    p <- p + ggplot2::geom_hline(yintercept = coverage, linetype = 2, size = 1, alpha = 0.9, color = "Gray")
  }
  p <- p +
    ggplot2::xlab("Replicate") +
    ggplot2::ylab("Total reads") +
    ggplot2::scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = ggplot2:: element_text(angle = 45, hjust = 1))
  file_name <- paste0("total_reads.", plot_type)
  ggplot2::ggsave(file.path(output_folder, file_name), plot = p, width = 10, height = 7, dpi = 300)
  
  # Gets colors for different screens
  screen_colors <- NA
  n_colors <- length(unique(col_groups))
  if (n_colors < 10) {
    screen_colors <- list(group = RColorBrewer::brewer.pal(length(unique(col_groups)), "Set1"))
  } else {
    pal <- RColorBrewer::brewer.pal(9, "Set1")
    pal <- grDevices::colorRampPalette(pal)(n_colors)
    screen_colors <- list(group = pal)
  }
  names(screen_colors$group) <- unique(col_groups)
  
  # Gets PCCs for heatmap
  df <- df[,all_cols]
  if (log_scale) {
    for (col in colnames(df)) {
      df[,col] <- log2(df[,col] + 1)
    }
  }
  cor_mat <- data.matrix(stats::cor(df))
  
  # Gets annotation for heatmap
  col_groups <- data.frame("Screen" = col_groups)
  rownames(col_groups) <- colnames(df)
  colnames(cor_mat) <- colnames(df)
  
  # Gets color for heatmap values
  breaks <- seq(-1, 1, by = (1/150))
  pal <- grDevices::colorRampPalette(c("#7fbf7b", "#f7f7f7", "#af8dc3"))(n = length(breaks))
  
  # Plots heatmap of raw reads
  filename <- file.path(output_folder, paste0("reads_heatmap.", plot_type))
  pheatmap::pheatmap(cor_mat,
                     border_color = NA,
                     annotation_col = col_groups,
                     annotation_colors = screen_colors,
                     display_numbers = display_numbers,
                     color = pal, 
                     breaks = breaks,
                     filename = filename)
}

#' Plot read counts.
#'
#' Plots a histogram of read counts for a given column.
#'
#' @param df Reads dataframe.
#' @param col Name of column to plot.
#' @param log_scale If true, log-normalizes data.
#' @param pseudocount Pseudocounts to add to log-normalized data if specified (default 1).
#' @return A ggplot object.
plot_reads <- function(df, col, log_scale = TRUE, pseudocount = 1) {
  x_label <- paste(col, "log-normalized read counts")
  y_label <- "Number of read counts"
  if (log_scale) {
    df[,col] <- log2(df[,col] + 1)
    y_label <- "Number of log-normalized read counts"
  }
  p <- ggplot2::ggplot(df, ggplot2::aes_string(col)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::xlab(x_label) +
    ggplot2::ylab(y_label) +
    ggthemes::theme_tufte(base_size = 20)
  return(p)
}

#' Plots sample comparisons.
#'
#' Pretty-plots comparisons between two samples in a scatterplot.
#
#' @param df Reads or lfc dataframe.
#' @param xcol Name of column containing values to plot on the x-axis.
#' @param ycol Name of column containing values to plot on the y-axis.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param color_col Name of column to color points by (optional).
#' @param color_lab Name of color legend (optional, defaults to color_col).
#' @param print_cor If true, prints Pearson correlation between columns 
#'   (default FALSE).
#' @return A list of two elements. The first is a ggplot object and the
#'   second is the correlation between xcol and ycol, if print_cor = TRUE.
plot_samples <- function(df, xcol, ycol, xlab, ylab, 
                         color_col = NULL, color_lab = NULL,
                         print_cor = FALSE) {
  
  # Optionally prints Pearson correlation between given columns
  pcc <- NA
  if (print_cor) {
    pcc <- stats::cor(df[[xcol]], df[[ycol]])
    cat(paste("Pearson correlation between", xcol, "and", ycol, ":", pcc, "\n"))
  }
  
  # Makes plot
  p <- NULL
  if (is.null(color_col)) {
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = xcol, y = ycol)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2, size = 1) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggthemes::theme_tufte(base_size = 20) 
  } else {
    if (is.null(color_lab)) {
      color_lab <- color_col
    }
    p <- ggplot2::ggplot(df, ggplot2::aes_string(x = xcol, y = ycol, color = color_col)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2, size = 1) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7) +
      ggplot2::scale_color_gradientn(colors = c("blue", "gray"), name = color_lab) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggthemes::theme_tufte(base_size = 20)
  }
  return(list(p, pcc))
}

#' Plots drug response for scored data.
#' 
#' Pretty-plots response for data which does not use a derived null model (e.g. for directly
#' comparing drug response to DMSO response). Assumes that data was scored by 
#' \code{score_conditions_vs_control} and significant effects were called by 
#' \code{call_significant_response}.
#' 
#' @param scores Dataframe of scores returned from \code{call_significant_response}.
#' @param control_name Name of control passed to \code{call_significant_response}.
#' @param condition_name Name of condition passed to \code{call_significant_response}.
#' @param loess If true and data was loess-normalized, plots loess null model instead
#'   (default TRUE).
#' @param neg_type Label for significant effects with a negative differential effect
#'   passed to \code{call_significant_response} (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   passed to \code{call_significant_response} (default "Positive").
#' @return A ggplot object.
#' @export
plot_significant_response <- function(scores, control_name, condition_name,
                                      loess = TRUE, neg_type = "Negative", 
                                      pos_type = "Positive") {
  
  # Manually sets colors for plot
  scores <- scores[order(scores[[paste0("differential_", condition_name, "_vs_", control_name)]]),]
  colors <- c("Black", "Gray", "Black")
  fill <- c("Blue", "Gray", "Yellow")
  response_col <- paste0("effect_type_", condition_name)
  neg_ind <- scores[[paste0("differential_", condition_name, "_vs_", control_name)]] < 0 &
    scores[[response_col]] != "None"
  pos_ind <- scores[[paste0("differential_", condition_name, "_vs_", control_name)]] > 0 &
    scores[[response_col]] != "None"
  if (any(neg_ind) & !any(pos_ind)) {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c(neg_type, "None"))
    colors <- c("Blue", "Gray")
    fill <- c("Black", "Gray")
  } else if (!any(neg_ind) & any(pos_ind)) {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c("None", pos_type))
    colors <- c("Gray", "Black")
    fill <- c("Gray", "Yellow")
  } else if (!any(neg_ind) & !(any(pos_ind))) {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c("None"))
    colors <- c("Gray")
    fill <- c("Gray")
  } else {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c(neg_type, "None", pos_type))
  }

  # Builds basic plot
  p <- ggplot2::ggplot(scores, ggplot2::aes_string(x = paste0("mean_", control_name), 
                                          y = paste0("mean_", condition_name))) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray")
  
  # Appends choice of null model to plot
  if (loess) {
  } else {
    p <- p + 
      ggplot2::geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black")
  }
  
  # Finishes plot
  p <- p + 
    ggplot2::geom_point(ggplot2::aes_string(color = response_col, fill = response_col), shape = 21, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::xlab(paste0(control_name, " mean log FC")) +
    ggplot2::ylab(paste0(condition_name, " mean log FC")) +
    ggplot2::labs(fill = "Significant response") +
    ggplot2::guides(color = FALSE, size = FALSE) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = ggplot2:: element_text(color = "Black", size = 16),
                   axis.text.y = ggplot2:: element_text(color = "Black", size = 16),
                   legend.text = ggplot2:: element_text(size = 16))
  return(p)
}

#' Plots drug response for scored paired data.
#' 
#' Pretty-plots response for data which uses a derived null model (e.g. for comparing
#' dual-gene knockout effects to a multiplicative null model derived from single-gene
#' effects). Assumes that data was scored by \code{score_combn_vs_single} and 
#' significant effects were called by \code{call_significant_response_combn}.
#' 
#' @param scores Dataframe of scores returned from \code{call_significant_response_combn}.
#' @param condition_name Name of condition passed to \code{call_significant_response_combn}.
#' @param filter_names If a list of column names is given, calls points as non-significant 
#'   if they are significant in the provided columsn (e.g. to remove points significant 
#'   in control screens; default NULL).
#' @param loess If true and data was loess-normalized, plots loess null model instead
#'   (default TRUE).
#' @param neg_type Label for significant effects with a negative differential effect
#'   passed to \code{call_significant_response_combn} (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   passed to \code{call_significant_response_combn} (default "Positive").
#' @return A ggplot object.
#' @export
plot_significant_response_combn <- function(scores, condition_name, filter_names = NULL, 
                                            loess = TRUE, neg_type = "Negative", 
                                            pos_type = "Positive") {
  
  # If filter_names given, remove significant guides with the same effect as a control 
  # type of guides (e.g. DMSO) from plot
  response_col <- paste0("effect_type_", condition_name)
  if (!is.null(filter_names)) {
    for (name in filter_names) {
      control_response_col <- paste0("effect_type_", name)
      scores[[response_col]][
        scores[[response_col]] != "None" & 
          scores[[control_response_col]] == scores[[response_col]]] <-
        "None"
    }
  }
  
  # Manually sets colors for plot
  scores <- scores[order(scores[[paste0("differential_combn_vs_single_", condition_name)]]),]
  colors <- c("Black", "Gray", "Black")
  fill <- c("Blue", "Gray", "Yellow")
  neg_ind <- scores[[paste0("differential_combn_vs_single_", condition_name)]] < 0 &
    scores[[response_col]] != "None"
  pos_ind <- scores[[paste0("differential_combn_vs_single_", condition_name)]] > 0 &
    scores[[response_col]] != "None"
  if (any(neg_ind) & !any(pos_ind)) {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c(neg_type, "None"))
    colors <- c("Blue", "Gray")
    fill <- c("Black", "Gray")
  } else if (!any(neg_ind) & any(pos_ind)) {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c("None", pos_type))
    colors <- c("Gray", "Black")
    fill <- c("Gray", "Yellow")
  } else if (!any(neg_ind) & !(any(pos_ind))) {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c("None"))
    colors <- c("Gray")
    fill <- c("Gray")
  } else {
    scores[[response_col]] <- factor(scores[[response_col]], levels = c(neg_type, "None", pos_type))
  }
  
  # Plots data
  p <- ggplot2::ggplot(scores, ggplot2::aes_string(x = paste0("mean_single_", condition_name), 
                                          y = paste0("mean_combn_", condition_name))) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray")
  
  # Appends choice of null model to plot
  if (loess) {
  } else {
    p <- p + 
      ggplot2::geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black")
  }
  
  # Finishes plot
  p <- p +
    ggplot2::geom_point(ggplot2::aes_string(color = response_col, fill = response_col), shape = 21, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::xlab(paste0(condition_name, " mean expected single-targeted log FC")) +
    ggplot2::ylab(paste0(condition_name, " mean observed combinatorial-targeted log FC")) +
    ggplot2::labs(fill = "Significant response") +
    ggplot2::guides(color = FALSE, size = FALSE) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = ggplot2:: element_text(color = "Black", size = 16),
                   axis.text.y = ggplot2:: element_text(color = "Black", size = 16),
                   legend.text = ggplot2:: element_text(size = 16))
  return(p)
}

#' Plot LFCs for all gene pairs.
#' 
#' Plots replicate comparisons for all replicates in a list of screens and outputs
#' plots to a given folder. Works for data returned from \code{call_significant_response}.
#' 
#' @param scores Dataframe of scores returned from \code{call_significant_response}.
#' @param residuals Residuals returned with the return_residuals argument set to true
#'   from \code{call_significant_response}.
#' @param control_name Name of control passed to \code{call_significant_response}.
#' @param condition_name Name of condition passed to \code{call_significant_response}.
#' @param output_folder Folder to output plots to. 
#' @param neg_type Label for significant effects with a negative differential effect
#'   passed to \code{call_significant_response} (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   passed to \code{call_significant_response} (default "Positive").
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export 
plot_lfc <- function(scores, residuals, control_name, condition_name, output_folder,
                     neg_type = "Negative", pos_type = "Positive", plot_type = "png") {
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  # Makes output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  # Gets top hits
  response_col <- paste0("effect_type_", condition_name)
  control_col <- paste0("mean_", control_name)
  condition_col <- paste0("mean_", condition_name)
  diff_col <- paste0("differential_", condition_name, "_vs_", control_name)
  scores <- scores[scores[[response_col]] != "None",]
  residuals <- residuals[residuals$n %in% as.numeric(rownames(scores)),]
  residuals$lfc <- residuals[[condition_col]] - residuals[[control_col]]
  
  # Gets ranking of top hits
  neg_order <- order(scores[[diff_col]])
  scores$neg_rank <- NA
  scores$pos_rank <- NA
  scores$neg_rank[neg_order] <- 1:nrow(scores)
  scores$pos_rank[neg_order] <- nrow(scores):1
  
  # Makes LFC plots for all top hits
  for (i in unique(residuals$n)) {
    
    # Gets data and gene names
    df <- residuals[residuals$n == i,]
    ind <- which(as.numeric(rownames(scores)) == i)
    gene1 <- scores$gene1[ind]
    gene2 <- scores$gene2[ind]
    x_label <- paste0("Guides")
    y_label <- paste0("Average differential LFC across replicates")
    
    # Adds ID column for plotting
    df$ID <- paste("Guide", 1:nrow(df))
    
    # Plots data
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, size = 1, alpha = 0.75, color = "Yellow") +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 0.75, color = "Gray") +
      ggplot2::geom_hline(yintercept = -1, linetype = 2, size = 1, alpha = 0.75, color = "Blue") +
      ggplot2::xlab(x_label) +
      ggplot2::ylab(y_label) +
      ggplot2::geom_bar(ggplot2::aes_string(x = "ID", y = "lfc"), stat = "identity", color = "Black", 
                        fill = ggplot2::alpha(c("gray30"), 1)) +
      ggplot2::coord_flip() +
      ggthemes::theme_tufte(base_size = 20)
    
    # Gets type and rank of effect
    effect <- ""
    rank <- 0
    effect_type <- scores[[response_col]][ind]
    if (effect_type == neg_type) {
      effect <- "neg"
      rank <- scores$neg_rank[ind]
    } else {
      effect <- "pos"
      rank <- scores$pos_rank[ind]
    }
    
    # Saves to file
    file_name <- paste0(effect, "_", rank, "_", gene1, "_", gene2, ".", plot_type)
    ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)
  }
}

#' Plot LFCs for all gene pairs.
#' 
#' Plots replicate comparisons for all replicates in a list of screens and outputs
#' plots to a given folder. Works for data returned from \code{call_significant_response_combn}.
#' 
#' @param scores Dataframe of scores returned from \code{call_significant_response_combn}.
#' @param residuals Residuals returned with the return_residuals argument set to true
#'   from \code{call_significant_response_combn}.
#' @param condition_name Name of condition passed to \code{call_significant_response_combn}.
#' @param output_folder Folder to output plots to. 
#' @param neg_type Label for significant effects with a negative differential effect
#'   passed to \code{call_significant_response_combn} (default "Negative").
#' @param pos_type Label for significant effects with a positive differential effect
#'   passed to \code{call_significant_response_combn} (default "Positive").
#' @param plot_type Type of plot to output, one of "png" or "pdf" (default "png").
#' @export 
plot_lfc_combn <- function(scores, residuals, condition_name, output_folder,
                           neg_type = "Negative", pos_type = "Positive",
                           plot_type = "png") {
  
  # Checks plot type and converts to lowercase
  plot_type <- tolower(plot_type)
  if (plot_type != "png" & plot_type != "pdf") {
    stop("plot_type must be either png or pdf")
  }
  
  # Makes output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  # Gets top hits across both orientations
  response_col <- paste0("effect_type_", condition_name)
  single_col <- paste0("mean_single_", condition_name)
  combn_col <- paste0("mean_combn_", condition_name)
  diff_col <- paste0("differential_combn_vs_single_", condition_name)
  scores <- scores[scores[[response_col]] != "None",]
  residuals1 <- residuals[[1]]
  residuals2 <- residuals[[2]]
  residuals1 <- residuals1[residuals1$n %in% as.numeric(rownames(scores)),]
  residuals2 <- residuals2[residuals2$n %in% as.numeric(rownames(scores)),]
  residuals1$lfc <- residuals1[[combn_col]] - residuals1[[single_col]]
  residuals2$lfc <- residuals2[[combn_col]] - residuals2[[single_col]]
  residuals1$orientation <- "Orientation 1"
  residuals2$orientation <- "Orientation 2"
  
  # Gets ranking of top hits
  neg_order <- order(scores[[diff_col]])
  scores$neg_rank <- NA
  scores$pos_rank <- NA
  scores$neg_rank[neg_order] <- 1:nrow(scores)
  scores$pos_rank[neg_order] <- nrow(scores):1
  
  # Makes LFC plots for all top hits
  for (i in unique(residuals1$n)) {
    
    # Gets data
    df1 <- residuals1[residuals1$n == i,]
    df2 <- residuals2[residuals2$n == i,]
    df1$orientation <- paste0(df1$orientation[1], " (n = ", nrow(df1), ")")
    df2$orientation <- paste0(df2$orientation[1], " (n = ", nrow(df2), ")")
    df <- rbind(df1, df2)
    
    # Gets gene names and sets axis labels
    ind <- which(as.numeric(rownames(scores)) == i)
    gene1 <- scores$gene1[ind]
    gene2 <- scores$gene2[ind]
    x_label <- paste0("Guides")
    y_label <- paste0("Average differential LFC across replicates")
    
    # Adds ID column for plotting
    df$ID <- nrow(df):1
    
    # Plots data
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, size = 1, alpha = 0.75, color = "Yellow") +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 0.75, color = "Gray") +
      ggplot2::geom_hline(yintercept = -1, linetype = 2, size = 1, alpha = 0.75, color = "Blue") +
      ggplot2::xlab(x_label) +
      ggplot2::ylab(y_label) +
      ggplot2::geom_bar(ggplot2::aes_string(x = "ID", y = "lfc"), stat = "identity", color = "Black", 
                        fill = ggplot2::alpha(c("gray30"), 1)) +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(. ~ orientation) +
      ggthemes::theme_tufte(base_size = 20) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA, color = "black"))
    
    # Gets type and rank of effect
    effect <- ""
    rank <- 0
    effect_type <- scores[[response_col]][ind]
    if (effect_type == neg_type) {
      effect <- "neg"
      rank <- scores$neg_rank[ind]
    } else {
      effect <- "pos"
      rank <- scores$pos_rank[ind]
    }
      
    # Saves to file
    file_name <- paste0(effect, "_", rank, "_", gene1, "_", gene2, ".", plot_type)
    ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)
  }
}


