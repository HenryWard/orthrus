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
#' @export 
plot_rep_comparisons <- function(df, screens, output_folder) {
  for (screen in screens) {
    rep_cols <- screen[["replicates"]]
    if (length(rep_cols) > 1) {
      pairs <- combn(rep_cols, 2)
      for (i in 1:ncol(pairs)) {
        col1 <- pairs[1,i]
        col2 <- pairs[2,i]
        x_label <- paste0(col1, " log fold change")
        y_label <- paste0(col2, " log fold change")
        p <- plot_samples(df, col1, col2, x_label, y_label, print_cor = TRUE)
        file_name <- paste0(col1, "_vs_", col2, "_replicate_comparison.png")
        ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)
      }   
    }
  }
}

#' Plot read counts for a screen.
#' 
#' Plots a histogram of read counts for each replicate of all screens.
#' 
#' @param df Reads dataframes.
#' @param screens List of screens created with \code{add_screens}.
#' @param output_folder Folder to output plots to. 
#' @param log_scale If true, log-normalizes data.
#' @param pseudocount Pseudocounts to add to log-normalized data if specified (default 1).
#' @export
plot_screen_reads <- function(df, screens, output_folder, 
                              log_scale = TRUE, pseudocount = 1) {
  for (screen in screens) {
    for (col in screen[["replicates"]]) {
      p <- plot_reads(df, col, log_scale, pseudocount)
      file_name <- paste0(col, "_raw_reads_histogram.png")
      ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)  
    } 
  }
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
  p <- ggplot2::ggplot(df, aes_string(col)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::xlab(x_label) +
    ggplot2::ylab(y_label) +
    ggthemes::theme_tufte(base_size = 20)
  return(p)
}

#' Plots sample comparisons.
#'
#' Pretty-plots comparisons between two samples in a scatterplot.
#'
#' @param df Reads or lfc dataframe.
#' @param xcol Name of column containing values to plot on the x-axis.
#' @param ycol Name of column containing values to plot on the y-axis.
#' @param color_col Name of column to color points by (optional).
#' @param color_lab Name of color legend (optional, defaults to color_col).
#' @param print_cor If true, prints Pearson correlation between columns 
#'   (default FALSE).
#' @return A ggplot object.
#' @export
plot_samples <- function(df, xcol, ycol, xlab, ylab, 
                         color_col = NULL, color_lab = NULL,
                         print_cor = FALSE) {
  
  # Optionally prints Pearson correlation between given columns
  if (print_cor) {
    pcc <- cor(df[[xcol]], df[[ycol]])
    cat(paste("Pearson correlation between", xcol, "and", ycol, ":", pcc, "\n"))
  }
  
  # Makes plot
  p <- NULL
  if (is.null(color_col)) {
    p <- ggplot2::ggplot(df, aes_string(x = xcol, y = ycol)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2, size = 1) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggthemes::theme_tufte(base_size = 20) 
  } else {
    if (is.null(color_lab)) {
      color_lab <- color_col
    }
    p <- ggplot2::ggplot(df, aes_string(x = xcol, y = ycol, color = color_col)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black", linetype = 2, size = 1) +
      ggplot2::geom_point(size = 1.5, alpha = 0.7) +
      ggplot2::scale_color_gradientn(colors = c("blue", "gray"), name = color_lab) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggthemes::theme_tufte(base_size = 20)
  }
  return(p)
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
                                      loess = TRUE, neg_type = "Negative", pos_type = "Positive") {
  
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
  p <- ggplot2::ggplot(scores, aes_string(x = paste0("mean_", control_name), 
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
    ggplot2::geom_point(aes_string(color = response_col, fill = response_col), shape = 21, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::xlab(paste0(control_name, " mean log FC")) +
    ggplot2::ylab(paste0(condition_name, " mean log FC")) +
    ggplot2::labs(fill = "Significant response") +
    ggplot2::guides(color = FALSE, size = FALSE) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = element_text(color = "Black", size = 16),
                   axis.text.y = element_text(color = "Black", size = 16),
                   legend.text = element_text(size = 16))
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
                                            loess = TRUE, neg_type = "Negative", pos_type = "Positive") {
  
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
  p <- ggplot2::ggplot(scores, aes_string(x = paste0("mean_single_", condition_name), 
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
    ggplot2::geom_point(aes_string(color = response_col, fill = response_col), shape = 21, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::xlab(paste0(condition_name, " mean expected single-targeted log FC")) +
    ggplot2::ylab(paste0(condition_name, " mean observed combinatorial-targeted log FC")) +
    ggplot2::labs(fill = "Significant response") +
    ggplot2::guides(color = FALSE, size = FALSE) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = element_text(color = "Black", size = 16),
                   axis.text.y = element_text(color = "Black", size = 16),
                   legend.text = element_text(size = 16))
  return(p)
}

#' Plot LFCs for all gene pairs.
#' 
#' Plots replicate comparisons for all replicates in a list of screens and outputs
#' plots to a given folder.
#' 
#' @param scores Dataframe of scores returned from \code{call_significant_response}.
#' @param residuals Residuals returned with the return_residuals argument set to true
#'   from \code{call_significant_response}.
#' @param control_name Name of control passed to \code{call_significant_response}.
#' @param condition_name Name of condition passed to \code{call_significant_response}.
#' @param output_folder Folder to output plots to. 
#' @export 
plot_lfc <- function(scores, residuals, control_name, condition_name, output_folder) {
  
  # Gets top hits
  response_col <- paste0("effect_type_", condition_name)
  control_col <- paste0("mean_", control_name)
  condition_col <- paste0("mean_", condition_name)
  scores <- scores[scores[[response_col]] != "None",]
  residuals <- residuals[residuals$n %in% as.numeric(rownames(scores)),]
  residuals$lfc <- residuals[[condition_col]] - residuals[[control_col]]
  
  # Makes LFC plots for all top hits
  for (i in unique(residuals$n)) {
    
    # Gets data and gene names
    df <- residuals[residuals$n == i,]
    ind <- which(as.numeric(rownames(scores)) == i)
    gene1 <- scores$gene1[ind]
    gene2 <- scores$gene2[ind]
    x_label <- paste0("Guides")
    y_label <- paste0("Average LFC across replicates")
    
    # ADds ID column for plotting
    df$ID <- paste("Guide", 1:nrow(df))
    
    # Plots data
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_hline(yintercept = 1, linetype = 2, size = 1, alpha = 0.75, color = "Yellow") +
      ggplot2::geom_hline(yintercept = -1, linetype = 2, size = 1, alpha = 0.75, color = "Blue") +
      ggplot2::xlab(x_label) +
      ggplot2::ylab(y_label) +
      ggplot2::geom_bar(aes(x = ID, y = lfc), stat = "identity", color = "Black", fill = alpha(c("gray30"), .9)) +
      ggplot2::coord_flip() +
      ggthemes::theme_tufte(base_size = 20)
    
    # Saves to file
    file_name <- paste0(gene1, "_", gene2, "_lfc.png")
    ggplot2::ggsave(file.path(output_folder, file_name), width = 10, height = 7, dpi = 300)
  }
}


