######
# UTILITY FUNCTIONS
######

#' Adds a list of screens from a sample table file
#' 
#' For a sample table formatted as a .tsv file with columns for screen, replicates,
#' and the name of a screen to compute log fold-changes against, adds each screen
#' to a list and returns that list.
#' 
#' @param table_file A tab-separated sample table with three columns, described above: 
#'   Screen, Replicates, and NormalizeTo. Screen is the name of the screen, replicates 
#'   are each technical replicate separated by semicolons, and NormalizeTo is either 
#'   the screen to normalize against or NA if unnecessary (e.g. for T0 screens).
#' @return A named list corresponding to provided screen names, where each sub-list 
#'   contains a list of the replicate columns (in "replicates") and the screen to 
#'   normalize against (in "normalize_name").
#' @export
add_screens_from_table <- function(table_file) {
  table <- utils::read.csv(table_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Checks columns of table
  cols <- c("Screen", "Replicates", "NormalizeTo")
  for (col in cols) {
    if (!(col %in% colnames(table))) {
      stop(paste("column", col, "not in sample table"))
    } else {
      if (col != "NormalizeTo") {
        if (any(is.na(table[[col]])))
          stop(paste("column", col, "contains NA values"))
      }
    }
  }
  
  # Adds all screens in table to list and returns
  screens <- NULL
  for (i in 1:nrow(table)) {
    name <- table[["Screen"]][i]
    replicates <- unlist(strsplit(table[["Replicates"]][i], ";"))
    normalize_name <- table[["NormalizeTo"]][i]
    if (is.na(normalize_name)) {
      normalize_name <- NULL
    }
    if (i == 1) {
      screens <- add_screen(name = name, replicates = replicates, normalize_name = normalize_name)
    } else {
      screens <- add_screen(screens, name, replicates, normalize_name)
    }
  }
  return(screens)
}

#' Adds a new screen to a list of screens
#' 
#' Makes a list containing info for a given screen and optionally appends it to a given list. 
#' 
#' @param screen_list An existing list of screens (optional, default NULL).
#' @param name A name for the screen (e.g. "RPE1_T18").
#' @param replicates A list of columns containing replicate data for the given screen.
#' @param normalize_name A name of another screen to normalize data in this screen to
#'   (e.g. "RPE1_T0", optional, default NULL). 
#' @return A named list corresponding to provided screen names, where each sub-list 
#'   contains a list of the replicate columns (in "replicates") and the screen to 
#'   normalize against (in "normalize_name").
#' @export
add_screen <- function(screen_list = NULL, name, replicates, 
                       normalize_name = NULL) {
  
  # Checks arguments
  if (is.na(name) | !is.character(name)) {
    stop("name must be a string")
  } else if (any(is.na(replicates)) | length(replicates) < 1 | !is.character(replicates)) {
    stop("replicates must be a list of one or more column names containing technical replicate readcounts")
  }
  
  # Adds arguments to list
  screen <- list()
  screen[["replicates"]] <- replicates
  screen[["normalize_name"]] <- normalize_name
  if (!is.null(screen_list)) {
    screen_list[[name]] <- screen
  } else {
    screen_list <- list()
    screen_list[[name]] <- screen
  }
  return(screen_list)
}

#' Removes a screen from a list of screens
#' 
#' Removes a screen with a given name from a list of screens. If multiple 
#' screens match the provided name, they are all removed.
#' 
#' @param screen_list An existing list of screens.
#' @param name A name for the screen to remove (e.g. "RPE1_T18").
#' @return \code{screen_list} with the named screen removed.
#' @export
remove_screen <- function(screen_list, name) {
  
  # Checks input
  if (length(name) > 1) {
    stop("name must be a single screen name contained in screen_list")
  } else if (is.name(name) | is.null(name)) {
    stop("name must be a single screen name contained in screen_list")
  }
  
  # Checks that screen is in screen_list
  ind <- names(screen_list) == name
  if (sum(ind)) {
    screen_list <- screen_list[!ind]
  } else {
    warning(paste("screen", name, "not found in screen_list"))
  }
  
  # Re-formats screen_list as empty list instead of empty named list
  if (length(screen_list) == 0) {
    screen_list <- list()
    message("No screens remaining in list. Returning empty list")
  }
  return(screen_list)
}

#' Checks input parameters
#' 
#' Checks to make sure that the given screens match the given dataframe.
#' 
#' @param df Reads or LFC dataframe.
#' @param screens List of screens created with \code{add_screens}.
#' @return TRUE.
#' @keywords internal
check_screen_params <- function(df, screens) {
  
  # Checks that the dataframe contains two columns named "Gene1" and "Gene2"
  if (!("gene1" %in% colnames(df))) {
    stop(paste("gene ID column gene1 not in dataframe column names"))
  } 
  if (!("gene2" %in% colnames(df))) {
    stop(paste("gene ID column gene2 not in dataframe column names"))
  }
  
  # Gets all replicates
  reps <- c()
  for (screen in screens) {
    reps <- c(reps, screen[["replicates"]])
  }
  
  # Checks input
  if (!all(reps %in% colnames(df))) {
    ind <- which(!(reps %in% colnames(df)))
    rep_name <- reps[ind[1]]
    screen_name <- ""
    for (name in names(screens)) {
      if (rep_name %in% screens[[name]][["replicates"]]) {
        screen_name <- name
      }
    }
    stop(paste("replicate", rep_name, "not in df, remove screen", screen_name, "with remove_screens"))
  }
  return(TRUE)
}


