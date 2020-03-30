######
# PLOTTING CODE
######

#' Plot read counts.
#'
#' Plots a histogram of read counts for a given column.
#'
#' @param df Reads dataframe.
#' @param col Name of column to plot.
#' @param log_scale If true, log-normalizes data.
#' @param pseudocount Pseudocounts to add to log-normalized data if specified (default 1).
#' @return A ggplot object.
#' @export
plot_reads <- function(df, col, log_scale = TRUE, pseudocount = 1) {
  y_label <- "Number of read counts"
  if (log_scale) {
    df[,col] <- log2(df[,col] + 1)
    y_label <- "Number of log-normalized read counts"
  }
  p <- ggplot2::ggplot(df, aes_string(col)) +
    ggplot2::geom_histogram(bins = 30) +
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
#' @return A ggplot object.
#' @export
plot_samples <- function(df, xcol, ycol, xlab, ylab, 
                         color_col = NULL, color_lab = NULL) {
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
#' @return A ggplot object.
#' @export
plot_significant_response <- function(scores, control_name, condition_name) {
  response_col <- paste0("effect_type_", condition_name)
  scores$response_factor <- ifelse(scores[[response_col]] != "None", 1, 0)
  scores <- scores[order(scores$response_factor, decreasing = FALSE),]
  p <- ggplot2::ggplot(scores, aes_string(x = paste0("mean_", control_name), 
                                          y = paste0("mean_", condition_name))) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black") +
    ggplot2::geom_point(aes_string(color = response_col, fill = response_col), shape = 21, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("Gray", "Black", "Black")) +
    ggplot2::scale_fill_manual(values = c("Gray", "Blue", "Yellow")) +
    ggplot2::xlab(paste(control_name, " mean log FC")) +
    ggplot2::ylab(paste(condition_name, " mean log FC")) +
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
#' effects). Assumes that data was scored by \code{score_dual_vs_single} and 
#' significant effects were called by \code{call_significant_response_dual}.
#' 
#' @param scores Dataframe of scores returned from \code{call_significant_response_dual}.
#' @param condition_name Name of condition passed to \code{call_significant_response_dual}.
#' @param filter_name If specified, calls points as non-significant if they are significant
#'   in this column (e.g. to remove points significant in DMSO screens). Default NULL.
#' @return A ggplot object.
#' @export
plot_significant_response_dual <- function(scores, condition_name, filter_name = NULL) {
  
  # If filter_name given, remove significant guides with the same effect as a control 
  # type of guides (e.g. DMSO) from plot
  condition_response_col <- paste0("effect_type_", condition_name)
  if (!is.null(filter_name)) {
    control_response_col <- paste0("effect_type_", filter_name)
    scores[[condition_response_col]][
      scores[[condition_response_col]] != "None" & 
        scores[[control_response_col]] == scores[[condition_response_col]]] <-
      "None"
  }
  
  # Plots data
  scores$response_factor <- ifelse(scores[[condition_response_col]] != "None", 1, 0)
  scores <- scores[order(scores$response_factor, decreasing = FALSE),]
  p <- ggplot2::ggplot(scores, aes_string(x = paste0("mean_single_", condition_name), 
                                          y = paste0("mean_dual_", condition_name))) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, size = 1, alpha = 1, color = "Gray") +
    ggplot2::geom_abline(slope = 1, intercept = 0, size = 1.5, alpha = 0.5, color = "Black") +
    ggplot2::geom_point(aes_string(color = condition_response_col, fill = condition_response_col), shape = 21, alpha = 0.7) +
    ggplot2::scale_color_manual(values = c("Gray", "Black", "Black")) +
    ggplot2::scale_fill_manual(values = c("Gray", "Blue", "Yellow")) +
    ggplot2::xlab(paste(condition_name, " mean expected single-targeted log FC")) +
    ggplot2::ylab(paste(condition_name, " mean dual-targeted log FC")) +
    ggplot2::labs(fill = "Significant response") +
    ggplot2::guides(color = FALSE, size = FALSE) +
    ggthemes::theme_tufte(base_size = 20) +
    ggplot2::theme(axis.text.x = element_text(color = "Black", size = 16),
                   axis.text.y = element_text(color = "Black", size = 16),
                   legend.text = element_text(size = 16))
  return(p)
}