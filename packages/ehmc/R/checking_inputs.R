checking_inputs <- function(names_input, input, default) {
  # --- Check for list inputs (general types) compared to default values
  for (name in names_input) {
    input_value <- input[[name]]
    default_value <- default[[name]]
    if (is.numeric(default_value)) {
      # --- Check for numeric (double) inputs
      if (!is.numeric(input_value) || length(input_value) != 1L || input_value < 0) {
        stop(name, "' must be a positive numeric value")
      }
    } else if (is.logical(default_value)) {
      # --- Check for boolean inputs
      if (!is.logical(input_value) || length(input_value) != 1L) {
        stop(name, "' must be a boolean (TRUE or FALSE)")
      }
    } else if (is.matrix(default_value)) {
      # --- Existing matrix and vector checks (same as before)
      if (is.null(dim(input_value))) {
        if (length(default_value) == 1L && length(input_value) == 1L) {
          input[[name]] <- matrix(input[[name]], nrow = 1)
          input_value <- input[[name]]
        } else {
          stop(name, "' must be a matrix")
        }
      }
      if (nrow(input_value) != nrow(default_value) || ncol(input_value) != ncol(default_value)) {
        stop(name, "' dimensions are incompatible. Expected ",
             nrow(default_value), "x", ncol(default_value),
             " but got ", nrow(input_value), "x", ncol(input_value))
      }
    } else if (length(default_value) > 1) {
      # --- Check for vector inputs (non-matrix)
      if (length(input_value) != length(default_value)) {
        stop(name, "' length is incompatible. Expected length",
             length(default_value), " but got length ", length(input_value))
      }
    }
  }
}
