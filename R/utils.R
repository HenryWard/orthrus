######
# UTILITY FUNCTIONS
######

#' Adds a new screen to a list of screens
#' 
#' Makes a list containing info for a given screen and optionally appends it to a given list. 
#' 
#' @param screen_list An existing list of screens (optional, default NULL).
#' @param name A name for the screen (e.g. "RPE1_T18").
#' @param replicates A list of columns containing replicate data for the given screen.
#' @param normalize_name A name of another screen to normalize data in this screen to
#'   (e.g. "RPE1_T0", optional, default NULL). 
#' @param target_coverage Target coverage of the given screen (optional, default 100).
#' @return A named list corresponding to provided screen names, where each sub-list 
#'   contains a list of the replicate columns (in "replicates") and the screen to 
#'   normalize against (in "normalize_name").
#' @export
add_screen <- function(screen_list = NULL, name, replicates, 
                       normalize_name = NULL, target_coverage = 200) {
  
  # Checks arguments
  if (is.null(target_coverage) | is.na(target_coverage) | is.character(target_coverage)) {
    stop("ERROR: target_coverage must be a numeric value (e.g. 100 for 100x coverage)")
  }
  
  # Adds arguments to list
  screen <- list()
  screen[["replicates"]] <- replicates
  screen[["normalize_name"]] <- normalize_name
  screen[["target_coverage"]] <- target_coverage
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
  ind <- names(screen_list) == name
  if (sum(ind)) {
    screen_list <- screen_list[!ind]
  }
  return(screen_list)
}
