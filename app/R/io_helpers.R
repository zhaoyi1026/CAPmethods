# =============================================================================
# Data I/O helpers
# -----------------------------------------------------------------------------
# Reusable parsers that translate the files a user uploads (a native R list of
# subject matrices, or tidy CSVs) into the internal structures the original CAP
# R code expects (a list of n matrices, each T_i x p, plus a covariate matrix).
# =============================================================================

#' Read any uploaded file into an R object.
#'
#' - .rds                  -> the stored object (e.g. a list of matrices)
#' - .RData / .rda         -> the single stored object, or a named list of all
#'                            objects if the file holds several
#' - .csv / .tsv / .txt    -> a data.frame
read_upload <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("rdata", "rda")) {
    e <- new.env()
    load(path, envir = e)
    objs <- ls(e)
    if (length(objs) == 1) return(get(objs[1], envir = e))
    return(mget(objs, envir = e))
  }
  read_table_any(path, name)
}

#' Coerce an uploaded response object into a list of subject matrices.
#'
#' Accepts: a long-format data.frame (id + p response columns), a list of
#' Ti x p matrices, or a named container (.RData) holding such a list under
#' one of Y / Y_list / Ylist / response.
#' @return list(Y, ids, Tvec, p, var_names) -- same shape as parse_Y_long().
coerce_Y_input <- function(obj) {
  if (is.data.frame(obj)) return(parse_Y_long(obj))

  # unwrap a named container that holds the list under a known key
  is_matrix_list <- is.list(obj) &&
    all(vapply(obj, function(m) is.matrix(m) || is.data.frame(m), logical(1)))
  if (is.list(obj) && !is.null(names(obj)) && !is_matrix_list) {
    key <- intersect(c("Y", "Y_list", "Ylist", "response", "Ymat"), names(obj))
    if (length(key) == 0)
      stop("Uploaded .RData has multiple objects but none named Y / Y_list.")
    obj <- obj[[key[1]]]
  }

  if (!is.list(obj) ||
      !all(vapply(obj, function(m) is.matrix(m) || is.data.frame(m), logical(1)))) {
    stop("Response must be a list of n matrices (each T_i x p), or a long CSV.")
  }
  Y <- lapply(obj, function(m) {
    m <- as.matrix(m); storage.mode(m) <- "double"; m
  })
  if (length(unique(vapply(Y, ncol, integer(1)))) != 1L)
    stop("All subject matrices must have the same number of columns (p).")
  Tvec <- vapply(Y, nrow, integer(1))
  if (any(Tvec < 2))
    stop("Each subject needs at least 2 samples (matrix rows).")

  ids <- names(Y)
  if (is.null(ids) || any(ids == "")) ids <- paste0("S", seq_along(Y))
  names(Y) <- ids
  p <- ncol(Y[[1]])
  var_names <- colnames(Y[[1]])
  if (is.null(var_names)) var_names <- paste0("V", seq_len(p))
  list(Y = Y, ids = ids, Tvec = Tvec, p = p, var_names = var_names)
}

#' Coerce an uploaded covariate object into a matrix aligned to subject order.
#'
#' Accepts: a data.frame with a leading id column (aligned by id), a plain
#' matrix/data.frame with one row per subject in list order, or a named .RData
#' container holding it under X / covariates / Xcov.
coerce_X_input <- function(obj, ids, add_intercept = TRUE) {
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj)) {
    key <- intersect(c("X", "covariates", "Xcov", "Z"), names(obj))
    if (length(key) == 0)
      stop("Could not find a covariate object (X) in the uploaded .RData.")
    obj <- obj[[key[1]]]
  }
  # data.frame whose first column lists the subject ids -> align by id
  if (is.data.frame(obj) && all(as.character(ids) %in% as.character(obj[[1]]))) {
    return(parse_X(obj, ids = ids, add_intercept = add_intercept))
  }
  if (is.data.frame(obj)) {
    nm <- names(obj)
    obj <- as.matrix(sapply(obj, function(x) suppressWarnings(as.numeric(x))))
    colnames(obj) <- nm
  }
  X <- as.matrix(obj); storage.mode(X) <- "double"
  if (anyNA(X)) stop("Covariates contain non-numeric or missing values.")
  if (nrow(X) != length(ids))
    stop(sprintf("Covariates have %d rows but there are %d subjects (rows must be in list order).",
                 nrow(X), length(ids)))
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))
  if (add_intercept) X <- cbind(Intercept = 1, X)
  X
}

#' Read an uploaded file (csv/tsv) into a data.frame.
read_table_any <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  df <- utils::read.csv(path, sep = sep, header = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < 2) {
    # maybe wrong separator guess; retry with the other one
    alt <- if (sep == ",") "\t" else ","
    df2 <- utils::read.csv(path, sep = alt, header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df2) > ncol(df)) df <- df2
  }
  df
}

#' Parse a "long format" response file into a list of subject matrices.
#'
#' Expected layout: first column = subject id, remaining columns = the p
#' response variables. Rows are observations stacked across subjects.
#'
#' @return list(Y = list of Ti x p matrices, ids = subject ids in order,
#'              Tvec = obs per subject, p = #responses, var_names = colnames)
parse_Y_long <- function(df, id_col = NULL) {
  if (is.null(id_col)) id_col <- names(df)[1]
  if (!id_col %in% names(df)) {
    stop(sprintf("Subject id column '%s' not found in response file.", id_col))
  }
  ids_raw <- df[[id_col]]
  ymat <- df[, setdiff(names(df), id_col), drop = FALSE]

  # coerce responses to numeric
  ymat[] <- lapply(ymat, function(x) suppressWarnings(as.numeric(x)))
  if (any(vapply(ymat, function(x) all(is.na(x)), logical(1)))) {
    stop("One or more response columns are non-numeric.")
  }

  ids <- unique(ids_raw)
  Y <- lapply(ids, function(i) as.matrix(ymat[ids_raw == i, , drop = FALSE]))
  names(Y) <- as.character(ids)
  Tvec <- vapply(Y, nrow, integer(1))

  if (any(Tvec < 2)) {
    stop("Each subject needs at least 2 observations to estimate a covariance.")
  }
  list(Y = Y, ids = ids, Tvec = Tvec, p = ncol(ymat), var_names = names(ymat))
}

#' Parse a subject-level covariate file into a matrix aligned to Y subject order.
#'
#' Expected layout: first column = subject id, remaining columns = covariates,
#' exactly one row per subject. If `add_intercept` is TRUE a leading column of
#' 1s named "Intercept" is prepended.
#'
#' @param ids subject id order to align to (from parse_Y_long)
parse_X <- function(df, ids, id_col = NULL, add_intercept = TRUE) {
  if (is.null(id_col)) id_col <- names(df)[1]
  if (!id_col %in% names(df)) {
    stop(sprintf("Subject id column '%s' not found in covariate file.", id_col))
  }
  rownames(df) <- as.character(df[[id_col]])
  xvars <- setdiff(names(df), id_col)
  if (length(xvars) == 0) stop("Covariate file has no covariate columns.")

  missing_ids <- setdiff(as.character(ids), rownames(df))
  if (length(missing_ids) > 0) {
    stop(sprintf("Covariate file is missing %d subject(s): %s",
                 length(missing_ids),
                 paste(utils::head(missing_ids, 5), collapse = ", ")))
  }

  Xdf <- df[as.character(ids), xvars, drop = FALSE]
  X <- as.matrix(sapply(Xdf, function(x) suppressWarnings(as.numeric(x))))
  if (any(is.na(X))) stop("Covariate file contains non-numeric or missing values.")
  X <- matrix(X, nrow = length(ids))
  colnames(X) <- xvars

  if (add_intercept) {
    X <- cbind(Intercept = 1, X)
  }
  X
}

#' Turn a list of subject matrices back into a long data.frame (for previews
#' and for serving example data as downloadable CSV).
Y_list_to_long <- function(Y, id_name = "id", var_names = NULL) {
  if (is.null(var_names)) var_names <- paste0("V", seq_len(ncol(Y[[1]])))
  out <- do.call(rbind, lapply(seq_along(Y), function(i) {
    m <- as.data.frame(Y[[i]])
    names(m) <- var_names
    cbind(stats::setNames(data.frame(names(Y)[i] %||% i), id_name), m)
  }))
  rownames(out) <- NULL
  out
}
