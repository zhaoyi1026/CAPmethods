# =============================================================================
# Method registry
# -----------------------------------------------------------------------------
# The website is a "plugin" host. Every CAP-family method registers itself by
# calling register_method() with a single list describing:
#   - identity / display metadata + the manuscript PDF
#   - the data files the user uploads (data_inputs)
#   - the tunable parameters (params)        -> drives the input form
#   - example()   : returns a built-in demo dataset
#   - run()       : runs the method, returns a standardized result object
#   - summarize() : returns named list of data.frames (tables)
#   - plots()     : returns named list of {title, plot} (plotly/ggplot/base)
#   - explain     : path to a markdown writeup
#
# Adding a new method == dropping one file in methods/<id>/ that calls
# register_method(). No edits to app.R or the generic UI are required.
# =============================================================================

# Container, populated at load time by each method file.
.METHOD_REGISTRY <- new.env(parent = emptyenv())
.METHOD_REGISTRY$methods <- list()
.METHOD_REGISTRY$order   <- character(0)

#' Register a CAP-family method with the website.
#' @param spec a list, see field documentation above.
register_method <- function(spec) {
  required <- c("id", "name", "full_name", "run")
  missing <- setdiff(required, names(spec))
  if (length(missing) > 0) {
    stop(sprintf("Method '%s' is missing required field(s): %s",
                 spec$id %||% "<unknown>", paste(missing, collapse = ", ")))
  }
  # sensible defaults for optional fields
  spec$data_inputs <- spec$data_inputs %||% list()
  spec$params      <- spec$params %||% list()
  spec$tags        <- spec$tags %||% character(0)
  spec$status      <- spec$status %||% "ready"   # "ready" | "beta" | "planned"

  .METHOD_REGISTRY$methods[[spec$id]] <- spec
  if (!spec$id %in% .METHOD_REGISTRY$order) {
    .METHOD_REGISTRY$order <- c(.METHOD_REGISTRY$order, spec$id)
  }
  invisible(spec)
}

#' All registered methods, in registration order.
list_methods <- function() {
  .METHOD_REGISTRY$methods[.METHOD_REGISTRY$order]
}

#' Fetch one method spec by id.
get_method <- function(id) {
  .METHOD_REGISTRY$methods[[id]]
}

# Small null-coalescing helper used throughout the app.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
