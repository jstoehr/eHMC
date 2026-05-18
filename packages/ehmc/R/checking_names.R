checking_names <- function(names_input, list_name, default = NULL, required = NULL) {
  # Check for additional elements not in the default list
  unknown <- names_input[!names_input %in% c(default, required)]
  if (length(unknown)) {
    warning("ignored elements in ", list_name, ": ", paste(unknown, collapse = ", "))
  }

  # Check for missing required elements
  missing <- required[!required %in% names_input]
  if (length(missing)) {
    stop("missing elements in ", list_name, " (required): ", paste(missing, collapse = ", "))
  }
  return(names_input[names_input %in% default])
}
