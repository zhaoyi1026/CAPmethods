# =============================================================================
# UI helpers
# -----------------------------------------------------------------------------
# Turn a method's declarative `params` / `data_inputs` specs into Shiny inputs,
# and collect their values back out. Keeps the generic method module agnostic
# to any particular method's parameters.
# =============================================================================

#' Build a single Shiny control from a parameter spec entry.
#' Supported types: numeric, integer, select, checkbox, slider, text.
make_param_control <- function(ns, p) {
  inputId <- ns(paste0("param_", p$id))
  lab <- p$label %||% p$id
  help <- p$help
  ctrl <- switch(
    p$type %||% "numeric",
    numeric = numericInput(inputId, lab, value = p$default,
                           min = p$min, max = p$max, step = p$step %||% NA),
    integer = numericInput(inputId, lab, value = p$default,
                           min = p$min, max = p$max, step = 1),
    select  = selectInput(inputId, lab, choices = p$choices,
                          selected = p$default),
    checkbox = checkboxInput(inputId, lab, value = isTRUE(p$default)),
    slider  = sliderInput(inputId, lab, min = p$min, max = p$max,
                          value = p$default, step = p$step %||% 1),
    text    = textInput(inputId, lab, value = p$default %||% ""),
    numericInput(inputId, lab, value = p$default)
  )
  if (!is.null(help)) {
    ctrl <- tagList(ctrl, tags$small(class = "form-text text-muted", help))
  }
  ctrl
}

#' Build the full parameter panel for a method.
make_param_panel <- function(ns, params) {
  if (length(params) == 0) {
    return(tags$p(class = "text-muted", "This method has no tunable parameters."))
  }
  lapply(params, function(p) make_param_control(ns, p))
}

#' Read the current values of all params back into a named list, coercing types.
collect_params <- function(input, params) {
  vals <- list()
  for (p in params) {
    v <- input[[paste0("param_", p$id)]]
    if (is.null(v)) v <- p$default
    if ((p$type %||% "numeric") %in% c("integer")) v <- as.integer(v)
    vals[[p$id]] <- v
  }
  vals
}

#' Turn a display label into a safe file-name stem (handles greek letters).
safe_filename <- function(x) {
  x <- gsub("β", "beta", x); x <- gsub("γ", "gamma", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (nchar(x) == 0) "data" else x
}

#' A small status badge for a method (ready / beta / planned).
status_badge <- function(status) {
  cls <- switch(status,
                ready = "text-bg-success",
                beta = "text-bg-warning",
                planned = "text-bg-secondary",
                "text-bg-secondary")
  lab <- switch(status,
                ready = "Ready", beta = "Beta", planned = "Planned", status)
  span(class = paste("badge", cls), lab)
}

#' Render a data.frame as a DT datatable with sensible defaults.
nice_table <- function(df, digits = 4, ...) {
  num <- vapply(df, is.numeric, logical(1))
  dt <- DT::datatable(
    df, rownames = !is.null(rownames(df)) && !identical(rownames(df), as.character(seq_len(nrow(df)))),
    options = list(dom = "tip", pageLength = 15, scrollX = TRUE),
    class = "compact stripe hover", ...
  )
  if (any(num)) dt <- DT::formatRound(dt, columns = which(num), digits = digits)
  dt
}
