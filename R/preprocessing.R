# ==============================================================================
# 00_preprocessing.R
# ------------------------------------------------------------------------------
# Argument-validation routines for the EXPORTED functions of PEAXAI.
#
# Scope of this file: these functions check ONLY that the arguments supplied by
# the user are well specified (presence, type, length, admissible values and
# cross-argument consistency). They do NOT inspect the contents of the data
# matrix; that is the responsibility of preprocessing(), in R/preprocessing.R.
#
# Admissible values follow Sections 2.2.1-2.2.3 of the SoftwareX manuscript.
#
# Contents:
#   1. peaxai_require_arg()                       presence of a mandatory argument
#   2. validate_parametes_PEAXAI_fitting()        PEAXAI_fitting()
#   3. shared helpers                             data / index / model / model variables /
#                                                 thresholds / targets / importance_method / seed
#   4. validate_parametes_PEAXAI_global_importance()   PEAXAI_global_importance()
#   5. validate_parametes_PEAXAI_local_importance()    PEAXAI_local_importance()
#   6. validate_parametes_PEAXAI_counterfactuals()     PEAXAI_counterfactuals()
#   7. validate_parametes_PEAXAI_ranking()        PEAXAI_ranking()
#   8. validate_parametes_PEAXAI_peer()           PEAXAI_peer()
#   9. validate_parametes_PEAXAI_predict()        PEAXAI_predict()
#  10. validate_parametes_reffcy()                reffcy()
#
# R/preprocessing.R keeps only validate_parametes_label_efficiency() (internal,
# not exported) and preprocessing() itself.
#
# Naming/collation note: this file is prefixed with "00_" so that it is the first
# one collated. Do not re-create any of the validators above in preprocessing.R:
# R collates alphabetically, so that file is sourced AFTER this one and its
# definition would silently overwrite the one here.
# ==============================================================================


# ==============================================================================
# Internal helper: presence of a mandatory argument
# ------------------------------------------------------------------------------
# PEAXAI_fitting() forwards every argument by name, e.g.
#   validate_parametes_PEAXAI_fitting(trControl = trControl, ...)
# so from inside the validator missing(trControl) is always FALSE: an argument
# *was* supplied, it just happens to be a promise pointing at a missing argument
# of the caller. Forcing that promise raises R's default
#   'argument "trControl" is missing, with no default'
# which points at the validator instead of at the user's call.
#
# This helper forces the promise inside tryCatch() so that we can emit our own
# message. It is only used for the arguments of PEAXAI_fitting() that have no
# default value: data, x, y, methods and trControl.
# ==============================================================================

peaxai_require_arg <- function(value, arg_name, hint = NULL) {

  err <- NULL

  supplied <- tryCatch(
    {
      force(value)
      TRUE
    },
    error = function(e) {
      err <<- conditionMessage(e)
      FALSE
    }
  )

  if (!supplied) {

    if (grepl("is missing, with no default", err, fixed = TRUE)) {
      msg <- sprintf("'%s' is required and has no default value.", arg_name)
      if (!is.null(hint)) {
        msg <- paste0(msg, " ", hint)
      }
      stop(msg, call. = FALSE)
    }

    stop(sprintf("'%s' could not be evaluated: %s", arg_name, err), call. = FALSE)
  }

  invisible(TRUE)
}


#' @title Validate PEAXAI_fitting input arguments
#'
#' @description
#' Checks that every argument supplied to \code{PEAXAI_fitting()} is well
#' specified, argument by argument, in the order in which they appear in the
#' function signature. Only argument specification is verified here: the
#' contents of \code{data} (missing values, admissible ranges) are checked by
#' \code{preprocessing()}.
#'
#' @param data A \code{data.frame} or \code{matrix}. Rows are DMUs.
#' @param x Integer vector with the column indices of the input variables.
#' @param y Integer vector with the column indices of the output variables.
#' @param RTS Returns-to-scale assumption: \code{"vrs"} / \code{1} or
#'   \code{"crs"} / \code{3}.
#' @param imbalance_rate \code{NULL} (no oversampling) or a numeric vector of
#'   target proportions of DMUs labelled as \code{'efficient'}, strictly in (0,1).
#' @param methods Named list of classification algorithms and their tuning grids.
#' @param trControl List defining the stratified cross-validation scheme:
#'   \code{list(method = "cv", number = <folds>)}.
#' @param metric_priority Character vector with the ordered performance criteria
#'   used to select the best model.
#' @param hold_out \code{NULL} or a single proportion in (0,1) reserved for the
#'   out-of-sample assessment.
#' @param seed Single integer-valued number used for reproducibility.
#' @param verbose Single logical controlling progress messages.
#'
#' @details
#' The admissible values follow Section 2.2.1 of the SoftwareX manuscript:
#' \itemize{
#'   \item \code{RTS}: "The available returns-to-scale assumptions are constant
#'         returns to scale ('crs') and variable returns to scale ('vrs')."
#'   \item \code{trControl}: "specifies the stratified cross-validation scheme
#'         and the number of folds".
#'   \item \code{hold_out}: "defines the proportion of data reserved through a
#'         stratified split for an additional out-of-sample assessment".
#'   \item \code{imbalance_rate}: "a numeric vector that allows users to evaluate
#'         multiple target proportions (...). When imbalance_rate = NULL, no
#'         oversampling is performed."
#' }
#'
#' Arguments belonging to unreleased versions of the package (\code{z_numeric},
#' \code{z_factor}, \code{B}, \code{alpha}, \code{m}) are deliberately absent:
#' they are forced to \code{NULL} in \code{PEAXAI_fitting()} and are not part of
#' the published interface.
#'
#' @return Invisibly \code{NULL}. Any violation raises \code{stop()} with a
#'   message identifying the offending argument.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_fitting <- function(
    data,
    x,
    y,
    RTS,
    imbalance_rate,
    methods,
    trControl,
    metric_priority,
    hold_out,
    seed,
    verbose
) {

  # ============================================================================
  # 1. data
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_fitting().
  # Expected: a non-empty data.frame or matrix. Rows are DMUs, columns are
  # variables. Manuscript, l. 113: "DMUs are provided as a data.frame or matrix".
  # ============================================================================

  peaxai_require_arg(
    data, "data",
    "It must be a data.frame or a matrix whose rows are DMUs."
  )

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop(
      sprintf("'data' must be a data.frame or a matrix; got an object of class '%s'.",
              paste(class(data), collapse = "/")),
      call. = FALSE
    )
  }

  if (nrow(data) == 0L) {
    stop("'data' has zero rows: at least one DMU is required.", call. = FALSE)
  }

  if (ncol(data) == 0L) {
    stop(
      "'data' must contain columns representing the input and output variables; an object with zero columns was provided.",
      call. = FALSE
    )
  }

  # Number of available columns, reused by the checks on 'x' and 'y'.
  n_col <- ncol(data)
  n_dmu <- nrow(data)


  # ============================================================================
  # 2. x  (input variable indices)
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_fitting().
  # Expected: a non-empty numeric vector of whole numbers, without NA, without
  # repetitions, and within 1:ncol(data).
  # ============================================================================

  peaxai_require_arg(
    x, "x",
    "It must be a vector with the column indices of the input variables, e.g. x = 1:3."
  )

  if (!is.numeric(x) || !is.null(dim(x))) {
    stop(
      sprintf("'x' must be a numeric vector of column indices corresponding to the input variables; got an object of class '%s'.",
              paste(class(x), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(x) == 0L) {
    stop("'x' is empty: at least one input variable is required.", call. = FALSE)
  }

  if (anyNA(x)) {
    stop("'x' contains NA values.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(
      "'x' contains non-finite values.",
      call. = FALSE
    )
  }

  # Whole numbers only: x = 1.5 would be silently truncated by data[, x].
  non_integer <- x != floor(x)

  if (any(non_integer)) {
    stop(
      sprintf(
        "'x' must contain integer-valued column indices corresponding to the input variables; got: %s.",
        paste(x[non_integer], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (any(duplicated(x))) {
    stop(
      sprintf("'x' contains repeated indices: %s.",
              paste(unique(x[duplicated(x)]), collapse = ", ")),
      call. = FALSE
    )
  }

  out_of_range <- x < 1L | x > n_col

  if (any(out_of_range)) {
    stop(
      sprintf(
        "'x' contains indices outside the valid range 1:%d: %s.",
        n_col,
        paste(x[out_of_range], collapse = ", ")
      ),
      call. = FALSE
    )
  }


  # ============================================================================
  # 3. y  (output variable indices)
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_fitting().
  # Same requirements as 'x', plus: it must not overlap with 'x'.
  # ============================================================================

  peaxai_require_arg(
    y, "y",
    "It must be a vector with the column indices of the output variables, e.g. y = 4."
  )

  if (!is.numeric(y) || !is.null(dim(y))) {
    stop(
      sprintf("'y' must be a numeric vector of column indices corresponding to the output variables; got an object of class '%s'.",
              paste(class(y), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(y) == 0L) {
    stop("'y' is empty: at least one output variable is required.", call. = FALSE)
  }

  if (anyNA(y)) {
    stop("'y' contains NA values.", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop(
      "'y' contains non-finite values.",
      call. = FALSE
    )
  }

  # Whole numbers only: y = 1.5 would be silently truncated by data[, y].
  non_integer <- y != floor(y)

  if (any(non_integer)) {
    stop(
      sprintf(
        "'y' must contain integer-valued column indices corresponding to the output variables; got: %s.",
        paste(y[non_integer], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (any(duplicated(y))) {
    stop(
      sprintf("'y' contains repeated indices: %s.",
              paste(unique(y[duplicated(y)]), collapse = ", ")),
      call. = FALSE
    )
  }

  out_of_range <- y < 1L | y > n_col

  if (any(out_of_range)) {
    stop(
      sprintf(
        "'y' contains indices outside the valid range 1:%d: %s.",
        n_col,
        paste(y[out_of_range], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # A variable cannot be an input and an output at the same time.
  if (any(x %in% y)) {
    stop(
      sprintf("'x' and 'y' must not overlap; shared indices: %s.",
              paste(intersect(x, y), collapse = ", ")),
      call. = FALSE
    )
  }


  # ============================================================================
  # 4. RTS
  # ----------------------------------------------------------------------------
  # Default in PEAXAI_fitting(): "vrs".
  # Manuscript, l. 135-136: "The available returns-to-scale assumptions are
  # constant returns to scale ('crs') and variable returns to scale ('vrs')."
  # The numeric codes are those of Benchmarking::dea.add(): vrs = 1, crs = 3.
  # ============================================================================

  RTS_available_char <- c("vrs", "crs")
  RTS_available_num  <- c(1, 3)          # 1 = vrs, 3 = crs

  if (is.null(RTS)) {
    stop("'RTS' cannot be NULL. Use \"vrs\" or \"crs\".", call. = FALSE)
  }

  if (length(RTS) != 1L) {
    stop(
      sprintf("'RTS' must be a single value; got a vector of length %d.", length(RTS)),
      call. = FALSE
    )
  }

  if (is.na(RTS)) {
    stop("'RTS' is NA. Use \"vrs\" or \"crs\".", call. = FALSE)
  }

  # Explicit branching on the two accepted types, with a final 'else' so that
  # logicals, factors and lists cannot slip through unchecked.
  if (is.character(RTS)) {

    if (!(RTS %in% RTS_available_char)) {
      stop(
        sprintf("'RTS' = \"%s\" is not recognized. Available: %s.",
                RTS, paste(sprintf("\"%s\"", RTS_available_char), collapse = ", ")),
        call. = FALSE
      )
    }

  } else if (is.numeric(RTS)) {

    if (!(RTS %in% RTS_available_num)) {
      stop(
        sprintf("'RTS' = %s is not recognized. Available codes: 1 (\"vrs\") and 3 (\"crs\").",
                format(RTS)),
        call. = FALSE
      )
    }

  } else {
    stop(
      sprintf("'RTS' must be a character string or a numeric code; got an object of class '%s'.",
              paste(class(RTS), collapse = "/")),
      call. = FALSE
    )
  }


  # ============================================================================
  # 5. imbalance_rate
  # ----------------------------------------------------------------------------
  # Default in PEAXAI_fitting(): NULL.
  # Manuscript, l. 147-158: "a numeric vector that allows users to evaluate
  # multiple target proportions of DMUs labeled as 'efficient' (...). When
  # imbalance_rate = NULL, no oversampling is performed."
  # Expected: NULL, or a numeric vector strictly inside (0,1).
  # ============================================================================

  if (!is.null(imbalance_rate)) {

    if (!is.numeric(imbalance_rate)) {
      stop(
        sprintf("'imbalance_rate' must be NULL or a numeric vector; got an object of class '%s'.",
                paste(class(imbalance_rate), collapse = "/")),
        call. = FALSE
      )
    }

    if (length(imbalance_rate) == 0L) {
      stop("'imbalance_rate' is empty. Use NULL to disable oversampling.", call. = FALSE)
    }

    if (anyNA(imbalance_rate)) {
      stop("'imbalance_rate' contains NA values.", call. = FALSE)
    }

    if (any(!is.finite(imbalance_rate))) {
      stop("'imbalance_rate' contains non-finite values.", call. = FALSE)
    }

    if (any(imbalance_rate <= 0 | imbalance_rate >= 1)) {
      stop(
        sprintf("'imbalance_rate' values must lie strictly in (0,1); offending values: %s.",
                paste(imbalance_rate[imbalance_rate <= 0 | imbalance_rate >= 1], collapse = ", ")),
        call. = FALSE
      )
    }

    # Repeated targets would fit the very same configuration twice.
    if (any(duplicated(imbalance_rate))) {
      stop(
        sprintf("'imbalance_rate' contains repeated values: %s.",
                paste(unique(imbalance_rate[duplicated(imbalance_rate)]), collapse = ", ")),
        call. = FALSE
      )
    }
  }


  # ============================================================================
  # 6. methods
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_fitting().
  # Manuscript, l. 161-163: "The argument methods specifies the classification
  # algorithms to be evaluated, together with the hyperparameter grids used for
  # model tuning. These algorithms are implemented via caret".
  # Expected: a named list; EVERY name must be a supported algorithm.
  # ============================================================================

  methods_available <- c("nnet", "svmPoly", "svmRadial", "rf", "glm")

  peaxai_require_arg(
    methods, "methods",
    "It must be a named list, e.g. methods = list(\"nnet\" = list(...))."
  )

  if (!is.list(methods)) {
    stop(
      sprintf("'methods' must be a list; got an object of class '%s'.",
              paste(class(methods), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(methods) == 0L) {
    stop("'methods' is empty: at least one classification algorithm is required.", call. = FALSE)
  }

  methods_names <- names(methods)

  if (is.null(methods_names) || any(!nzchar(methods_names))) {
    stop(
      sprintf("'methods' must be a NAMED list. Supported names: %s.",
              paste(methods_available, collapse = ", ")),
      call. = FALSE
    )
  }

  # all(), not any(): a single valid name must not validate the whole list.
  if (!all(methods_names %in% methods_available)) {
    stop(
      sprintf("'methods' contains unsupported algorithms: %s. Supported: %s.",
              paste(setdiff(methods_names, methods_available), collapse = ", "),
              paste(methods_available, collapse = ", ")),
      call. = FALSE
    )
  }

  if (any(duplicated(methods_names))) {
    stop(
      sprintf("'methods' contains repeated algorithms: %s.",
              paste(unique(methods_names[duplicated(methods_names)]), collapse = ", ")),
      call. = FALSE
    )
  }


  # ============================================================================
  # 7. trControl
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_fitting().
  # Manuscript, l. 163-165: "The training procedure is controlled through
  # trControl, which specifies the stratified cross-validation scheme and the
  # number of folds".
  # Expected: list(method = "cv", number = <integer >= 2>).
  # ============================================================================

  trControl_methods_available <- "cv"

  peaxai_require_arg(
    trControl, "trControl",
    "It must be a list, e.g. trControl = list(method = \"cv\", number = 5)."
  )

  if (!is.list(trControl)) {
    stop(
      sprintf("'trControl' must be a list; got an object of class '%s'. Example: list(method = \"cv\", number = 5).",
              paste(class(trControl), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(trControl) == 0L) {
    stop("'trControl' is empty. Example: list(method = \"cv\", number = 5).", call. = FALSE)
  }

  # --- trControl$method -------------------------------------------------------
  tr_method <- trControl[["method"]]

  if (is.null(tr_method)) {
    stop("'trControl$method' is missing. It must be \"cv\".", call. = FALSE)
  }

  if (!is.character(tr_method) || length(tr_method) != 1L || is.na(tr_method)) {
    stop("'trControl$method' must be a single character string.", call. = FALSE)
  }

  if (!(tr_method %in% trControl_methods_available)) {
    stop(
      sprintf("'trControl$method' = \"%s\" is not supported. PEAXAI uses stratified cross-validation: \"cv\".",
              tr_method),
      call. = FALSE
    )
  }

  # --- trControl$number (number of folds) -------------------------------------
  tr_number <- trControl[["number"]]

  if (is.null(tr_number)) {
    stop("'trControl$number' is missing. It must be the number of cross-validation folds, e.g. 5.",
         call. = FALSE)
  }

  if (!is.numeric(tr_number) || length(tr_number) != 1L || is.na(tr_number)) {
    stop("'trControl$number' must be a single numeric value.", call. = FALSE)
  }

  if (tr_number != as.integer(tr_number)) {
    stop(
      sprintf("'trControl$number' must be a whole number; got %s.", format(tr_number)),
      call. = FALSE
    )
  }

  if (tr_number < 2L) {
    stop(
      sprintf("'trControl$number' must be at least 2; got %d.", as.integer(tr_number)),
      call. = FALSE
    )
  }

  # Cross-argument consistency: more folds than DMUs is not feasible.
  if (tr_number > n_dmu) {
    stop(
      sprintf("'trControl$number' (%d folds) cannot exceed the number of DMUs in 'data' (%d).",
              as.integer(tr_number), n_dmu),
      call. = FALSE
    )
  }


  # ============================================================================
  # 8. metric_priority
  # ----------------------------------------------------------------------------
  # Default in PEAXAI_fitting(): "Balanced_Accuracy".
  # Manuscript, l. 174-175: "Candidate models are compared and selected according
  # to a performance criterion introduced in metric_priority."
  # It is an ORDERED vector: the first metric decides, the rest break ties
  # (manuscript example: c("Balanced_Accuracy", "F1", "ROC_AUC")).
  # ============================================================================

  metric_priority_available <- c(
    "Accuracy", "Kappa",
    "Recall", "Specificity", "Precision",
    "F1", "Balanced_Accuracy", "G_mean",
    "ROC_AUC", "PR_AUC"
  )

  if (is.null(metric_priority)) {
    stop(
      sprintf("'metric_priority' cannot be NULL. Available metrics: %s.",
              paste(metric_priority_available, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.character(metric_priority)) {
    stop(
      sprintf("'metric_priority' must be a character vector; got an object of class '%s'.",
              paste(class(metric_priority), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(metric_priority) == 0L) {
    stop("'metric_priority' is empty: at least one metric is required.", call. = FALSE)
  }

  if (anyNA(metric_priority)) {
    stop("'metric_priority' contains NA values.", call. = FALSE)
  }

  # all(), not any(): a single valid metric must not validate a misspelled one.
  if (!all(metric_priority %in% metric_priority_available)) {
    stop(
      sprintf("'metric_priority' contains unsupported metrics: %s. Available: %s.",
              paste(setdiff(metric_priority, metric_priority_available), collapse = ", "),
              paste(metric_priority_available, collapse = ", ")),
      call. = FALSE
    )
  }

  # Repeated criteria are meaningless as tie-breakers.
  if (any(duplicated(metric_priority))) {
    stop(
      sprintf("'metric_priority' contains repeated metrics: %s.",
              paste(unique(metric_priority[duplicated(metric_priority)]), collapse = ", ")),
      call. = FALSE
    )
  }


  # ============================================================================
  # 9. hold_out
  # ----------------------------------------------------------------------------
  # Default in PEAXAI_fitting(): NULL.
  # Manuscript, l. 165-166: "hold_out defines the proportion of data reserved
  # through a stratified split for an additional out-of-sample assessment."
  # Expected: NULL, or a single proportion strictly inside (0,1).
  # ============================================================================

  if (!is.null(hold_out)) {

    if (!is.numeric(hold_out)) {
      stop(
        sprintf("'hold_out' must be NULL or a single numeric proportion; got an object of class '%s'.",
                paste(class(hold_out), collapse = "/")),
        call. = FALSE
      )
    }

    if (length(hold_out) != 1L) {
      stop(
        sprintf("'hold_out' must be a single value; got a vector of length %d.", length(hold_out)),
        call. = FALSE
      )
    }

    if (is.na(hold_out) || !is.finite(hold_out)) {
      stop("'hold_out' must be a finite numeric value.", call. = FALSE)
    }

    if (hold_out <= 0 || hold_out >= 1) {
      stop(
        sprintf("'hold_out' must lie strictly in (0,1); got %s.", format(hold_out)),
        call. = FALSE
      )
    }
  }


  # ============================================================================
  # 10. seed
  # ----------------------------------------------------------------------------
  # Default in PEAXAI_fitting(): 314.
  # Manuscript, l. 175: "reproducibility is ensured through seed".
  # Expected: a single whole number, as required by set.seed().
  # ============================================================================

  if (is.null(seed)) {
    stop("'seed' cannot be NULL. Provide a single integer, e.g. seed = 314.", call. = FALSE)
  }

  if (!is.numeric(seed)) {
    stop(
      sprintf("'seed' must be a single numeric value; got an object of class '%s'.",
              paste(class(seed), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(seed) != 1L) {
    stop(
      sprintf("'seed' must be a single value; got a vector of length %d.", length(seed)),
      call. = FALSE
    )
  }

  if (is.na(seed) || !is.finite(seed)) {
    stop("'seed' must be a finite numeric value.", call. = FALSE)
  }

  if (seed != as.integer(seed)) {
    stop(
      sprintf("'seed' must be a whole number; got %s.", format(seed)),
      call. = FALSE
    )
  }


  # ============================================================================
  # 11. verbose
  # ----------------------------------------------------------------------------
  # Default in PEAXAI_fitting(): TRUE.
  # Manuscript, l. 176: "the verbose argument controls whether the training
  # progress is displayed."
  # Expected: a single non-NA logical.
  # ============================================================================

  if (is.null(verbose)) {
    stop("'verbose' cannot be NULL. Use TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("'verbose' must be a single logical value: TRUE or FALSE.", call. = FALSE)
  }


  # ============================================================================
  # 12. Cross-argument advisory (warning, not an error)
  # ----------------------------------------------------------------------------
  # Supplementary Material, Section 4.2, citing Cooper, Seiford & Tone (2007):
  # to limit the curse of dimensionality DEA requires
  #   n >= max{ m * s , 3 * (m + s) }
  # Falling below it does not prevent PEAXAI from running, but DEA tends to
  # label nearly every DMU as efficient, which is the failure mode documented in
  # Supplementary Section 4.1. Remove this block if you prefer no advisories.
  # ============================================================================

  n_inputs  <- length(x)
  n_outputs <- length(y)
  n_min     <- max(n_inputs * n_outputs, 3 * (n_inputs + n_outputs))

  if (n_dmu < n_min) {
    warning(
      sprintf(
        paste0("'data' has %d DMUs for %d inputs and %d outputs. The DEA literature ",
               "recommends at least %d DMUs (Cooper, Seiford & Tone, 2007). With fewer, ",
               "the DEA labelling stage may assign almost every DMU to the 'efficient' class."),
        n_dmu, n_inputs, n_outputs, n_min
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}


# ==============================================================================
# SHARED HELPERS
# ------------------------------------------------------------------------------
# Every exported function receives the same handful of argument shapes: a data
# container, index vectors, a fitted caret model, probability thresholds and a
# targets object. They are factored out here so that a fix lands in one place
# instead of in seven validators.
#
# NOTE: blocks 1-3 of validate_parametes_PEAXAI_fitting() in 00_preprocessing.R
# still carry their own inline copy of the data/index logic. Once you are happy
# with the helpers below, replace those blocks with calls to
# peaxai_check_data() / peaxai_check_index() so there is a single implementation.
# ==============================================================================


#' @title Validate a data container argument
#'
#' @description
#' Checks that a data argument is a non-empty \code{data.frame} or \code{matrix}.
#'
#' @param data Object to validate.
#' @param arg_name Name of the argument, used in the error messages.
#'
#' @return Invisibly a list with \code{n_row} and \code{n_col}.
#'
#' @keywords internal
#' @noRd

peaxai_check_data <- function(data, arg_name = "data") {

  if (!is.data.frame(data) && !is.matrix(data)) {
    stop(
      sprintf("'%s' must be a data.frame or a matrix whose rows are DMUs; got an object of class '%s'.",
              arg_name, paste(class(data), collapse = "/")),
      call. = FALSE
    )
  }

  if (nrow(data) == 0L) {
    stop(
      sprintf("'%s' has zero rows: at least one DMU is required.", arg_name),
      call. = FALSE
    )
  }

  if (ncol(data) == 0L) {
    stop(
      sprintf("'%s' must contain columns representing the input and output variables; an object with zero columns was provided.",
              arg_name),
      call. = FALSE
    )
  }

  invisible(list(n_row = nrow(data), n_col = ncol(data)))
}


#' @title Validate a vector of column indices
#'
#' @description
#' Checks that an index vector is a non-empty numeric vector of whole numbers,
#' without NA, without non-finite values, without repetitions, and within
#' \code{1:n_col}. Same logic used for \code{x} and \code{y} everywhere.
#'
#' @param index Object to validate.
#' @param arg_name Name of the argument, used in the error messages.
#' @param role Human-readable role, e.g. "input" or "output".
#' @param n_col Number of columns of the associated data container.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_index <- function(index, arg_name, role, n_col) {

  if (!is.numeric(index) || !is.null(dim(index))) {
    stop(
      sprintf("'%s' must be a numeric vector of column indices corresponding to the %s variables; got an object of class '%s'.",
              arg_name, role, paste(class(index), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(index) == 0L) {
    stop(
      sprintf("'%s' is empty: at least one %s variable is required.", arg_name, role),
      call. = FALSE
    )
  }

  if (anyNA(index)) {
    stop(sprintf("'%s' contains NA values.", arg_name), call. = FALSE)
  }

  if (any(!is.finite(index))) {
    stop(sprintf("'%s' contains non-finite values.", arg_name), call. = FALSE)
  }

  # Whole numbers only: 1.5 would be silently truncated by data[, index].
  non_integer <- index != floor(index)

  if (any(non_integer)) {
    stop(
      sprintf("'%s' must contain integer-valued column indices corresponding to the %s variables; got: %s.",
              arg_name, role, paste(index[non_integer], collapse = ", ")),
      call. = FALSE
    )
  }

  if (any(duplicated(index))) {
    stop(
      sprintf("'%s' contains repeated indices: %s.",
              arg_name, paste(unique(index[duplicated(index)]), collapse = ", ")),
      call. = FALSE
    )
  }

  out_of_range <- index < 1L | index > n_col

  if (any(out_of_range)) {
    stop(
      sprintf("'%s' contains indices outside the valid range 1:%d: %s.",
              arg_name, n_col, paste(index[out_of_range], collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' @title Validate that x and y do not overlap
#'
#' @description
#' A variable cannot be an input and an output at the same time.
#'
#' @param x Input indices, already validated.
#' @param y Output indices, already validated.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_xy_overlap <- function(x, y) {

  if (any(x %in% y)) {
    stop(
      sprintf("'x' and 'y' must not overlap; shared indices: %s.",
              paste(intersect(x, y), collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' @title Validate a fitted classifier argument
#'
#' @description
#' Checks that the model handed over to the efficiency-analysis functions is the
#' object returned by \code{PEAXAI_fitting()}, i.e. a \code{caret::train} fit.
#'
#' @param final_model Object to validate.
#' @param arg_name Name of the argument, used in the error messages.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_final_model <- function(final_model, arg_name = "final_model") {

  peaxai_require_arg(
    final_model, arg_name,
    "It is the fitted classifier returned by PEAXAI_fitting(), e.g. models[[\"best_model_fit\"]][[\"nnet\"]]."
  )

  if (is.null(final_model)) {
    stop(
      sprintf("'%s' cannot be NULL. Provide the fitted classifier returned by PEAXAI_fitting().", arg_name),
      call. = FALSE
    )
  }

  if (!inherits(final_model, "train")) {
    stop(
      sprintf("'%s' must be a caret::train() object (class 'train'); got an object of class '%s'.",
              arg_name, paste(class(final_model), collapse = "/")),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' @title Validate the variables of the data against those of the fitted model
#'
#' @description
#' Checks that the DMUs handed over to an efficiency-analysis or explanation
#' function carry every variable the classifier was trained on.
#'
#' @details
#' Variables are matched by their \strong{exact} name. \code{caret} sanitises
#' column names with \code{make.names()} when it builds its training data, so a
#' variable the user calls \code{"Total assets"} is stored inside
#' \code{final_model} as \code{"Total.assets"}. PEAXAI does \emph{not} try to
#' undo that transformation: \code{make.names()} is not injective (\code{"a b"}
#' and \code{"a.b"} both become \code{"a.b"}), so matching on sanitised names can
#' pair a column with the wrong variable and score the DMUs against the wrong
#' predictor without anyone noticing. Instead the variable is reported as
#' missing, under the name the model actually carries, and the user renames the
#' column.
#'
#' \code{data} may carry any number of additional columns: only the variables
#' the model was trained on are kept, in the order the model expects. The
#' positions returned are the ones the caller must use for that selection.
#'
#' If the fit does not carry its training columns the check is skipped, so it
#' can never turn a working call into a failing one.
#'
#' @param final_model A fitted \code{caret::train} object, already validated.
#' @param data Data container already validated.
#' @param x,y Index vectors already validated. When both are \code{NULL} every
#'   column of \code{data} is compared (used by \code{PEAXAI_predict()}).
#' @param arg_name Name of the data argument, used in the error messages.
#'
#' @return Invisibly an integer vector with the column positions of \code{data}
#'   that correspond to the model variables, in the order expected by
#'   \code{final_model} and named after them; \code{NULL} when the check is
#'   skipped.
#'
#' @keywords internal
#' @noRd

peaxai_check_model_variables <- function(final_model, data, x = NULL, y = NULL,
                                         arg_name = "data") {

  # ----------------------------------------------------------------------------
  # Variables the model was trained on -----------------------------------------
  # ----------------------------------------------------------------------------
  # trainingData first: it holds the actual training columns. xNames and
  # coefnames are fallbacks and may carry dummy-coded names instead.
  model_variables <- NULL

  training_data <- final_model[["trainingData"]]

  if (!is.null(training_data)) {
    model_variables <- setdiff(colnames(training_data), ".outcome")
  }

  if (length(model_variables) == 0L) {
    model_variables <- final_model[["finalModel"]][["xNames"]]
  }

  if (length(model_variables) == 0L) {
    model_variables <- final_model[["coefnames"]]
  }

  # The fit does not carry its training columns: nothing to compare against.
  if (length(model_variables) == 0L) {
    return(invisible(NULL))
  }

  # ----------------------------------------------------------------------------
  # Variables supplied by the user ---------------------------------------------
  # ----------------------------------------------------------------------------
  data_names     <- colnames(as.data.frame(data))
  data_positions <- seq_along(data_names)

  if (is.null(data_names)) {
    return(invisible(NULL))
  }

  if (!is.null(x) || !is.null(y)) {
    selected       <- c(x, y)
    data_names     <- data_names[selected]
    data_positions <- data_positions[selected]
  }

  # ----------------------------------------------------------------------------
  # Comparison on exact names --------------------------------------------------
  # ----------------------------------------------------------------------------
  # No make.names() on either side: see the @details of this helper. A variable
  # whose name does not match character by character is a missing variable.

  # Ambiguous names must not be matched silently.
  if (anyDuplicated(model_variables)) {
    stop(
      sprintf(
        "'final_model' was trained on repeated predictor names: %s.",
        paste(
          sprintf("'%s'", unique(model_variables[duplicated(model_variables)])),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }

  duplicated_data <- unique(data_names[duplicated(data_names)])
  ambiguous_names <- intersect(duplicated_data, model_variables)

  if (length(ambiguous_names) > 0L) {
    stop(
      sprintf(
        "'%s' contains more than one column named %s, so the variables of 'final_model' cannot be identified unambiguously.",
        arg_name,
        paste(sprintf("'%s'", ambiguous_names), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  # Positions of the data columns, in the order expected by final_model.
  position <- match(model_variables, data_names)

  missing_variables <- model_variables[is.na(position)]

  if (length(missing_variables) > 0L) {

    # A column that would match only after make.names() is the commonest cause:
    # point at it, but do not accept it. Renaming is the user's decision.
    near_matches <- data_names[
      make.names(data_names) %in% missing_variables &
        !(data_names %in% model_variables)
    ]

    hint <- if (length(near_matches) > 0L) {
      sprintf(
        paste0(
          " '%s' does contain %s, whose sanitised name would match. caret applies ",
          "make.names() to the column names when it trains the model, so rename the ",
          "column(s) to the exact name reported above."
        ),
        arg_name,
        paste(sprintf("'%s'", near_matches), collapse = ", ")
      )
    } else {
      ""
    }

    stop(
      sprintf(
        paste0(
          "'%s' does not contain every variable 'final_model' was trained on. Missing: %s. ",
          "Variables are matched by their exact name.%s"
        ),
        arg_name,
        paste(sprintf("'%s'", missing_variables), collapse = ", "),
        hint
      ),
      call. = FALSE
    )
  }

  # Original data positions, ordered and named as expected by final_model.
  # Any other column of 'data' is simply not selected.
  matched_positions <- data_positions[position]
  names(matched_positions) <- model_variables

  invisible(matched_positions)
}


#' @title Validate the efficiency-probability thresholds
#'
#' @description
#' Manuscript, Section 2.2.3: \code{efficiency_thresholds} "defines the
#' probability levels used in threshold-based efficiency analysis". They are
#' probabilities, hence strictly inside (0,1).
#'
#' @param efficiency_thresholds Object to validate.
#' @param arg_name Name of the argument, used in the error messages.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_thresholds <- function(efficiency_thresholds,
                                    arg_name = "efficiency_thresholds") {

  peaxai_require_arg(
    efficiency_thresholds, arg_name,
    "It must be a numeric vector of probabilities, e.g. seq(0.75, 0.95, 0.1)."
  )

  if (is.null(efficiency_thresholds)) {
    stop(
      sprintf("'%s' cannot be NULL. Provide a numeric vector of probabilities, e.g. seq(0.75, 0.95, 0.1).",
              arg_name),
      call. = FALSE
    )
  }

  if (!is.numeric(efficiency_thresholds) || !is.null(dim(efficiency_thresholds))) {
    stop(
      sprintf("'%s' must be a numeric vector of probabilities; got an object of class '%s'.",
              arg_name, paste(class(efficiency_thresholds), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(efficiency_thresholds) == 0L) {
    stop(
      sprintf("'%s' is empty: at least one probability threshold is required.", arg_name),
      call. = FALSE
    )
  }

  if (anyNA(efficiency_thresholds)) {
    stop(sprintf("'%s' contains NA values.", arg_name), call. = FALSE)
  }

  if (any(!is.finite(efficiency_thresholds))) {
    stop(sprintf("'%s' contains non-finite values.", arg_name), call. = FALSE)
  }

  outside <- efficiency_thresholds <= 0 | efficiency_thresholds >= 1

  if (any(outside)) {
    stop(
      sprintf("'%s' values must lie strictly in (0,1); offending values: %s.",
              arg_name, paste(efficiency_thresholds[outside], collapse = ", ")),
      call. = FALSE
    )
  }

  if (any(duplicated(efficiency_thresholds))) {
    stop(
      sprintf("'%s' contains repeated values: %s.",
              arg_name,
              paste(unique(efficiency_thresholds[duplicated(efficiency_thresholds)]), collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' @title Validate a targets object produced by PEAXAI_counterfactuals()
#'
#' @description
#' Manuscript, Section 2.2.3: \code{targets} "refers to the counterfactual
#' projections returned by \code{PEAXAI_counterfactuals()}". Each element must
#' carry a \code{counterfactual_dataset} whose variables match
#' \code{colnames(data)[c(x, y)]}.
#'
#' @param targets Object to validate.
#' @param data Data container already validated.
#' @param x,y Index vectors already validated.
#' @param arg_name Name of the argument, used in the error messages.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_targets <- function(targets, data, x, y, arg_name = "targets") {

  if (!is.list(targets) || is.data.frame(targets)) {
    stop(
      sprintf("'%s' must be the list returned by PEAXAI_counterfactuals(); got an object of class '%s'.",
              arg_name, paste(class(targets), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(targets) == 0L) {
    stop(sprintf("'%s' is an empty list.", arg_name), call. = FALSE)
  }

  expected_names <- colnames(as.data.frame(data))[c(x, y)]

  for (i in seq_along(targets)) {

    element_label <- if (!is.null(names(targets)) && nzchar(names(targets)[i])) {
      names(targets)[i]
    } else {
      paste0("[[", i, "]]")
    }

    counterfactual_dataset <- targets[[i]][["counterfactual_dataset"]]

    if (is.null(counterfactual_dataset)) {
      stop(
        sprintf("Element '%s' of '%s' has no 'counterfactual_dataset'. Pass the object returned by PEAXAI_counterfactuals().",
                element_label, arg_name),
        call. = FALSE
      )
    }

    # identical(), not ==: with different lengths '==' recycles and warns
    # instead of failing cleanly.
    if (!identical(colnames(as.data.frame(counterfactual_dataset)), expected_names)) {
      stop(
        sprintf("The variables of element '%s' of '%s' do not match colnames(data)[c(x, y)]: expected %s.",
                element_label, arg_name, paste(expected_names, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  invisible(NULL)
}


#' @title Validate the importance_method list
#'
#' @description
#' Shared by \code{PEAXAI_global_importance()} and \code{PEAXAI_local_importance()}.
#' Manuscript, Section 2.2.2, l. 223-245: the list carries the name of the
#' technique and its hyperparameters, and "if hyperparameters are omitted,
#' PEAXAI applies standard default values". Every hyperparameter is therefore
#' optional; only \code{name} is mandatory.
#'
#' @param importance_method Object to validate.
#' @param names_available Character vector with the techniques admitted by the
#'   calling function: \code{c("SA","SHAP","PI")} for global explanations and
#'   \code{c("SA","SHAP","LIME")} for local ones.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_importance_method <- function(importance_method, names_available) {

  peaxai_require_arg(
    importance_method, "importance_method",
    "It must be a named list, e.g. list(name = \"SHAP\", bg_n = 200)."
  )

  if (!is.list(importance_method)) {
    stop(
      sprintf("'importance_method' must be a list; got an object of class '%s'.",
              paste(class(importance_method), collapse = "/")),
      call. = FALSE
    )
  }

  # --- name -------------------------------------------------------------------
  name <- importance_method[["name"]]

  if (is.null(name)) {
    stop(
      sprintf("'importance_method$name' is missing. Available techniques: %s.",
              paste(names_available, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("'importance_method$name' must be a single character string.", call. = FALSE)
  }

  if (!(name %in% names_available)) {
    stop(
      sprintf("'importance_method$name' = \"%s\" is not available here. Available techniques: %s.",
              name, paste(names_available, collapse = ", ")),
      call. = FALSE
    )
  }

  # --- SA, via rminer ---------------------------------------------------------
  if (identical(name, "SA")) {

    sa_methods_available  <- c("1D-SA", "sens", "DSA", "MSA", "CSA", "GSA")
    sa_measures_available <- c("AAD", "gradient", "variance", "range")

    if (!is.null(importance_method[["method"]])) {
      sa_method <- importance_method[["method"]]
      if (!is.character(sa_method) || length(sa_method) != 1L ||
          !(sa_method %in% sa_methods_available)) {
        stop(
          sprintf("'importance_method$method' (SA/rminer) must be one of: %s.",
                  paste(sa_methods_available, collapse = ", ")),
          call. = FALSE
        )
      }
    }

    if (!is.null(importance_method[["measures"]])) {
      sa_measures <- importance_method[["measures"]]
      if (!is.character(sa_measures) ||
          !all(sa_measures %in% sa_measures_available)) {
        stop(
          sprintf("'importance_method$measures' (SA/rminer) must be in: %s.",
                  paste(sa_measures_available, collapse = ", ")),
          call. = FALSE
        )
      }
    }

    if (!is.null(importance_method[["levels"]])) {
      sa_levels <- importance_method[["levels"]]
      if (!is.numeric(sa_levels) || length(sa_levels) != 1L || is.na(sa_levels) ||
          sa_levels < 2 || sa_levels != floor(sa_levels)) {
        stop("'importance_method$levels' (SA/rminer) must be a single whole number >= 2.",
             call. = FALSE)
      }
    }

    if (!is.null(importance_method[["baseline"]])) {
      sa_baseline <- importance_method[["baseline"]]
      baseline_ok <- (is.character(sa_baseline) && length(sa_baseline) == 1L &&
                        sa_baseline %in% c("mean", "median")) ||
        is.data.frame(sa_baseline)
      if (!baseline_ok) {
        stop("'importance_method$baseline' (SA/rminer) must be \"mean\", \"median\", or a one-row data.frame with a specific DMU.",
             call. = FALSE)
      }
    }
  }

  # --- SHAP, via kernelshap ---------------------------------------------------
  if (identical(name, "SHAP")) {

    bg_n <- importance_method[["bg_n"]]

    if (!is.null(bg_n)) {
      if (!is.numeric(bg_n) || length(bg_n) != 1L || is.na(bg_n) ||
          bg_n < 1 || bg_n != floor(bg_n)) {
        stop("'importance_method$bg_n' (SHAP/kernelshap) must be a single whole number >= 1. Default: 200.",
             call. = FALSE)
      }
    }
  }

  # --- PI, via iml (global only) ----------------------------------------------
  if (identical(name, "PI")) {

    n_repetitions <- importance_method[["n.repetitions"]]

    if (!is.null(n_repetitions)) {
      if (!is.numeric(n_repetitions) || length(n_repetitions) != 1L || is.na(n_repetitions) ||
          n_repetitions < 1 || n_repetitions != floor(n_repetitions)) {
        stop("'importance_method$n.repetitions' (PI/iml) must be a single whole number >= 1. Default: 5.",
             call. = FALSE)
      }
    }
  }

  # --- LIME, via lime (local only) --------------------------------------------
  if (identical(name, "LIME")) {

    n_permutations <- importance_method[["n_permutations"]]

    if (!is.null(n_permutations)) {
      if (!is.numeric(n_permutations) || length(n_permutations) != 1L || is.na(n_permutations) ||
          n_permutations < 1 || n_permutations != floor(n_permutations)) {
        stop("'importance_method$n_permutations' (LIME/lime) must be a single whole number >= 1. Default: 5000.",
             call. = FALSE)
      }
    }

    if (!is.null(importance_method[["feature_select"]])) {
      feature_select <- importance_method[["feature_select"]]
      feature_select_available <- c("auto", "none", "forward_selection",
                                    "highest_weights", "lasso_path", "tree")
      if (!is.character(feature_select) || length(feature_select) != 1L ||
          !(feature_select %in% feature_select_available)) {
        stop(
          sprintf("'importance_method$feature_select' (LIME/lime) must be one of: %s. Default: \"auto\".",
                  paste(feature_select_available, collapse = ", ")),
          call. = FALSE
        )
      }
    }

    if (!is.null(importance_method[["bin_continuous"]])) {
      bin_continuous <- importance_method[["bin_continuous"]]
      if (!is.logical(bin_continuous) || length(bin_continuous) != 1L || is.na(bin_continuous)) {
        stop("'importance_method$bin_continuous' (LIME/lime) must be a single logical value. Default: TRUE.",
             call. = FALSE)
      }
    }

    if (!is.null(importance_method[["n_bins"]])) {
      n_bins <- importance_method[["n_bins"]]
      if (!is.numeric(n_bins) || length(n_bins) != 1L || is.na(n_bins) ||
          n_bins < 1 || n_bins != floor(n_bins)) {
        stop("'importance_method$n_bins' (LIME/lime) must be a single whole number >= 1. Default: 4.",
             call. = FALSE)
      }
    }
  }

  invisible(NULL)
}


#' @title Validate a seed argument
#'
#' @description
#' A single whole number, as required by \code{set.seed()}.
#'
#' @param seed Object to validate.
#' @param arg_name Name of the argument, used in the error messages.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

peaxai_check_seed <- function(seed, arg_name = "seed") {

  if (is.null(seed)) {
    stop(
      sprintf("'%s' cannot be NULL. Provide a single integer, e.g. %s = 314.", arg_name, arg_name),
      call. = FALSE
    )
  }

  if (!is.numeric(seed)) {
    stop(
      sprintf("'%s' must be a single numeric value; got an object of class '%s'.",
              arg_name, paste(class(seed), collapse = "/")),
      call. = FALSE
    )
  }

  if (length(seed) != 1L) {
    stop(
      sprintf("'%s' must be a single value; got a vector of length %d.", arg_name, length(seed)),
      call. = FALSE
    )
  }

  if (is.na(seed) || !is.finite(seed)) {
    stop(sprintf("'%s' must be a finite numeric value.", arg_name), call. = FALSE)
  }

  if (seed != floor(seed)) {
    stop(
      sprintf("'%s' must be a whole number; got %s.", arg_name, format(seed)),
      call. = FALSE
    )
  }

  invisible(NULL)
}


#' @title Validate PEAXAI_global_importance input arguments
#'
#' @description
#' Checks the arguments of \code{PEAXAI_global_importance()}, argument by
#' argument, in the order of the signature published in Section 2.2.2 of the
#' SoftwareX manuscript (l. 201-202).
#'
#' @param final_model The fitted classifier to be explained, of class
#'   \code{"train"}.
#' @param x Integer vector with the column indices of the input variables in
#'   \code{explain_data}.
#' @param y Integer vector with the column indices of the output variables in
#'   \code{explain_data}.
#' @param explain_data Target observations for which the explanations are
#'   computed, usually the observed, non-synthetic DMUs.
#' @param reference_data Background training dataset used as reference, which
#'   may include synthetic DMUs.
#' @param importance_method Named list selecting the explanation technique.
#'   Globally: \code{"SHAP"} (kernelshap), \code{"PI"} (iml) or \code{"SA"}
#'   (rminer).
#' @param seed Single integer-valued number used for reproducibility.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_global_importance <- function(
    final_model,
    x,
    y,
    explain_data,
    reference_data,
    importance_method,
    seed = 314
) {

  # ============================================================================
  # 1. final_model
  # ============================================================================

  peaxai_check_final_model(final_model, "final_model")

  # ============================================================================
  # 2. explain_data
  # ----------------------------------------------------------------------------
  # Validated before 'x' and 'y' because it defines the column range they must
  # fall into. Manuscript, l. 209-211: "contains the target observations
  # (usually the observed, non-synthetic DMUs) for which the explanations are
  # computed".
  # ============================================================================

  peaxai_require_arg(
    explain_data, "explain_data",
    "It must be a data.frame with the DMUs to be explained, e.g. explain_data = data."
  )

  explain_dim <- peaxai_check_data(explain_data, "explain_data")
  n_col       <- explain_dim[["n_col"]]

  # ============================================================================
  # 3. x  (input variable indices)
  # ============================================================================

  peaxai_require_arg(
    x, "x",
    "It must be a vector with the column indices of the input variables, e.g. x = 1:3."
  )

  peaxai_check_index(x, "x", "input", n_col)

  # ============================================================================
  # 4. y  (output variable indices)
  # ============================================================================

  peaxai_require_arg(
    y, "y",
    "It must be a vector with the column indices of the output variables, e.g. y = 4."
  )

  peaxai_check_index(y, "y", "output", n_col)
  peaxai_check_xy_overlap(x, y)

  # ============================================================================
  # 5. reference_data
  # ----------------------------------------------------------------------------
  # Manuscript, l. 211-213: "defines the background training dataset used as
  # reference for the explanation, which may include the augmented training set
  # with synthetic DMUs". It is usually final_model$trainingData, which carries
  # the outcome column as well, so the variables of 'explain_data' must be a
  # SUBSET of its columns, not identical to them.
  # ============================================================================

  peaxai_require_arg(
    reference_data, "reference_data",
    "It must be the training dataset of the fitted model, e.g. reference_data = final_model$trainingData."
  )

  peaxai_check_data(reference_data, "reference_data")

  explain_names   <- colnames(as.data.frame(explain_data))[c(x, y)]
  reference_names <- colnames(as.data.frame(reference_data))

  if (!is.null(explain_names) && !is.null(reference_names)) {

    # Exact names on both sides, as everywhere else in PEAXAI. caret sanitises
    # the column names when it builds trainingData ("Total assets" is stored as
    # "Total.assets"), so 'explain_data' must already use the sanitised name.
    missing_names <- explain_names[!(explain_names %in% reference_names)]

    if (length(missing_names) > 0L) {
      stop(
        sprintf("'reference_data' does not contain the variables of 'explain_data'[c(x, y)]; missing: %s.",
                paste(missing_names, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  # Every variable the model was trained on must be present in the DMUs to be
  # explained, under the exact name the model carries.
  peaxai_check_model_variables(final_model, explain_data, x, y, "explain_data")

  # ============================================================================
  # 6. importance_method
  # ----------------------------------------------------------------------------
  # Manuscript, l. 216-219: globally PEAXAI uses aggregated SHAP values,
  # Permutation Feature Importance (PI) and Sensitivity Analysis (SA).
  # ============================================================================

  peaxai_check_importance_method(
    importance_method,
    names_available = c("SHAP", "PI", "SA")
  )

  # ============================================================================
  # 7. seed
  # ============================================================================

  peaxai_check_seed(seed, "seed")

  invisible(NULL)
}


#' @title Validate PEAXAI_local_importance input arguments
#'
#' @description
#' Checks the arguments of \code{PEAXAI_local_importance()}. Identical to the
#' global validator except for the set of admissible techniques: locally PEAXAI
#' offers SHAP, LIME and SA, and not PI (manuscript, l. 219-222).
#'
#' @param final_model The fitted classifier to be explained, of class
#'   \code{"train"}.
#' @param x Integer vector with the column indices of the input variables in
#'   \code{explain_data}.
#' @param y Integer vector with the column indices of the output variables in
#'   \code{explain_data}.
#' @param explain_data DMUs for which the local explanations are computed.
#' @param reference_data Background training dataset used as reference.
#' @param importance_method Named list selecting the explanation technique.
#'   Locally: \code{"SHAP"} (kernelshap), \code{"LIME"} (lime) or \code{"SA"}
#'   (rminer).
#' @param seed Single integer-valued number used for reproducibility.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_local_importance <- function(
    final_model,
    x,
    y,
    explain_data,
    reference_data,
    importance_method,
    seed = 314
) {

  # ============================================================================
  # 1. final_model
  # ============================================================================

  peaxai_check_final_model(final_model, "final_model")

  # ============================================================================
  # 2. explain_data
  # ============================================================================

  peaxai_require_arg(
    explain_data, "explain_data",
    "It must be a data.frame with the DMUs to be explained, e.g. explain_data = data."
  )

  explain_dim <- peaxai_check_data(explain_data, "explain_data")
  n_col       <- explain_dim[["n_col"]]

  # ============================================================================
  # 3. x  (input variable indices)
  # ============================================================================

  peaxai_require_arg(
    x, "x",
    "It must be a vector with the column indices of the input variables, e.g. x = 1:3."
  )

  peaxai_check_index(x, "x", "input", n_col)

  # ============================================================================
  # 4. y  (output variable indices)
  # ============================================================================

  peaxai_require_arg(
    y, "y",
    "It must be a vector with the column indices of the output variables, e.g. y = 4."
  )

  peaxai_check_index(y, "y", "output", n_col)
  peaxai_check_xy_overlap(x, y)

  # ============================================================================
  # 5. reference_data
  # ============================================================================

  peaxai_require_arg(
    reference_data, "reference_data",
    "It must be the training dataset of the fitted model, e.g. reference_data = final_model$trainingData."
  )

  peaxai_check_data(reference_data, "reference_data")

  explain_names   <- colnames(as.data.frame(explain_data))[c(x, y)]
  reference_names <- colnames(as.data.frame(reference_data))

  if (!is.null(explain_names) && !is.null(reference_names)) {

    # Exact names on both sides, as everywhere else in PEAXAI. caret sanitises
    # the column names when it builds trainingData ("Total assets" is stored as
    # "Total.assets"), so 'explain_data' must already use the sanitised name.
    missing_names <- explain_names[!(explain_names %in% reference_names)]

    if (length(missing_names) > 0L) {
      stop(
        sprintf("'reference_data' does not contain the variables of 'explain_data'[c(x, y)]; missing: %s.",
                paste(missing_names, collapse = ", ")),
        call. = FALSE
      )
    }
  }

  # Every variable the model was trained on must be present in the DMUs to be
  # explained, under the exact name the model carries.
  peaxai_check_model_variables(final_model, explain_data, x, y, "explain_data")

  # ============================================================================
  # 6. importance_method
  # ----------------------------------------------------------------------------
  # Manuscript, l. 219-222: locally, explanations rely on SHAP via kernelshap,
  # LIME via lime, and SA via rminer. PI is a global-only technique.
  # ============================================================================

  peaxai_check_importance_method(
    importance_method,
    names_available = c("SHAP", "LIME", "SA")
  )

  # ============================================================================
  # 7. seed
  # ============================================================================

  peaxai_check_seed(seed, "seed")

  invisible(NULL)
}


#' @title Validate PEAXAI_counterfactuals input arguments
#'
#' @description
#' Checks the arguments of \code{PEAXAI_counterfactuals()}, in the order of the
#' signature published in Section 2.2.3 of the SoftwareX manuscript (l. 259-261).
#'
#' @param data A \code{data.frame} or \code{matrix} with the DMUs under analysis.
#' @param x Integer vector with the column indices of the input variables.
#' @param y Integer vector with the column indices of the output variables.
#' @param final_model The fitted classifier used to estimate efficiency
#'   probabilities, of class \code{"train"}.
#' @param efficiency_thresholds Numeric vector of probability levels, strictly
#'   in (0,1).
#' @param directional_vector List with \code{relative_importance} and
#'   \code{baseline}, defining the direction used in the counterfactual search.
#' @param n_expand Single non-negative number controlling how the initial search
#'   interval is expanded.
#' @param n_grid Single whole number >= 2 with the candidate counterfactuals
#'   evaluated at each iteration.
#' @param max_y Single non-negative number bounding the output expansion.
#' @param min_x Single non-negative number bounding the input contraction.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_counterfactuals <- function(
    data,
    x,
    y,
    final_model,
    efficiency_thresholds,
    directional_vector,
    n_expand,
    n_grid,
    max_y = 1,
    min_x = 1
) {

  # ============================================================================
  # 1. data
  # ============================================================================

  peaxai_require_arg(
    data, "data",
    "It must be a data.frame or a matrix whose rows are DMUs."
  )

  data_dim <- peaxai_check_data(data, "data")
  n_col    <- data_dim[["n_col"]]

  # ============================================================================
  # 2. x  (input variable indices)
  # ============================================================================

  peaxai_require_arg(
    x, "x",
    "It must be a vector with the column indices of the input variables, e.g. x = 1:3."
  )

  peaxai_check_index(x, "x", "input", n_col)

  # ============================================================================
  # 3. y  (output variable indices)
  # ============================================================================

  peaxai_require_arg(
    y, "y",
    "It must be a vector with the column indices of the output variables, e.g. y = 4."
  )

  peaxai_check_index(y, "y", "output", n_col)
  peaxai_check_xy_overlap(x, y)

  # ============================================================================
  # 4. final_model
  # ============================================================================

  peaxai_check_final_model(final_model, "final_model")

  # Every variable the model was trained on must be present in 'data', matched
  # by its exact name: a column caret renamed at training time must already
  # carry the sanitised name here.
  peaxai_check_model_variables(final_model, data, x, y, "data")

  # ============================================================================
  # 5. efficiency_thresholds
  # ============================================================================

  peaxai_check_thresholds(efficiency_thresholds, "efficiency_thresholds")

  # ============================================================================
  # 6. directional_vector
  # ----------------------------------------------------------------------------
  # Manuscript, l. 271-283: a list defining the direction used in the
  # counterfactual search, with two fields:
  #   - relative_importance: the relative importance of the variables, obtained
  #     from PEAXAI_global_importance() / PEAXAI_local_importance() or supplied
  #     by the user;
  #   - baseline: the unit of the direction vector in the DDF.
  # ============================================================================

  peaxai_require_arg(
    directional_vector, "directional_vector",
    "It must be a list, e.g. list(relative_importance = ri, baseline = \"mean\")."
  )

  if (!is.list(directional_vector) || is.data.frame(directional_vector)) {
    stop(
      sprintf("'directional_vector' must be a list; got an object of class '%s'.",
              paste(class(directional_vector), collapse = "/")),
      call. = FALSE
    )
  }

  # --- directional_vector$relative_importance ---------------------------------
  # Two admissible shapes, with DIFFERENT rules, because the two importance
  # functions return different things:
  #
  #   * GLOBAL (one row, or a plain numeric vector). PEAXAI_global_importance()
  #     returns unsigned magnitudes, so PEAXAI_counterfactuals() is the one that
  #     imposes the economic sign (inputs down, outputs up). The weights must
  #     therefore be non-negative, and they add up to 1 over the whole vector.
  #
  #   * LOCAL (one row per DMU). PEAXAI_local_importance() re-imposes the sign of
  #     the underlying SHAP/LIME value, so the weights are SIGNED and already
  #     carry the direction in which each variable must move for that DMU. What
  #     is normalised to 1 is the sum of the ABSOLUTE values of each row, not the
  #     plain sum, which may be negative or zero for a perfectly valid matrix.
  #
  # The scope is decided exactly as PEAXAI_counterfactuals() decides it: by the
  # number of rows.
  relative_importance <- directional_vector[["relative_importance"]]

  if (is.null(relative_importance)) {
    stop("'directional_vector$relative_importance' is missing. It is the object returned by PEAXAI_global_importance() or PEAXAI_local_importance().",
         call. = FALSE)
  }

  # Normalise the container to a numeric matrix. A plain vector is read as a
  # single (global) row.
  if (is.data.frame(relative_importance)) {
    relative_importance <- as.matrix(relative_importance)
  } else if (is.numeric(relative_importance) && is.null(dim(relative_importance))) {
    relative_importance <- matrix(relative_importance, nrow = 1L)
  }

  if (!is.matrix(relative_importance) || !is.numeric(relative_importance)) {
    stop(
      sprintf("'directional_vector$relative_importance' must be a numeric vector, matrix or data.frame; got an object of class '%s'.",
              paste(class(directional_vector[["relative_importance"]]), collapse = "/")),
      call. = FALSE
    )
  }

  n_variables <- length(c(x, y))

  if (ncol(relative_importance) != n_variables) {
    stop(
      sprintf("'directional_vector$relative_importance' must have one weight per variable: expected %d columns (length(c(x, y))), got %d.",
              n_variables, ncol(relative_importance)),
      call. = FALSE
    )
  }

  if (anyNA(relative_importance) || any(!is.finite(relative_importance))) {
    stop("'directional_vector$relative_importance' contains NA or non-finite values.", call. = FALSE)
  }

  if (nrow(relative_importance) == 1L) {

    # --- global importance: unsigned magnitudes ------------------------------
    if (any(relative_importance < 0)) {
      stop(
        sprintf("'directional_vector$relative_importance' must be non-negative when a single (global) row is supplied; offending values: %s.",
                paste(relative_importance[relative_importance < 0], collapse = ", ")),
        call. = FALSE
      )
    }

    if (all(relative_importance == 0)) {
      stop("'directional_vector$relative_importance' cannot be the all-zero vector.", call. = FALSE)
    }

    if (abs(sum(relative_importance) - 1) > 1e-8) {
      warning(
        sprintf("'directional_vector$relative_importance' adds up to %s instead of 1; proceeding without normalization.",
                format(sum(relative_importance))),
        call. = FALSE
      )
    }

  } else {

    # --- local importance: signed weights, one row per DMU -------------------
    if (nrow(relative_importance) != nrow(data)) {
      stop(
        sprintf("'directional_vector$relative_importance' must have either 1 row (global importance) or one row per DMU (local importance); 'data' has %d DMUs and %d rows were supplied.",
                nrow(data), nrow(relative_importance)),
        call. = FALSE
      )
    }

    # Signs are admissible here: they encode the direction of movement. What
    # cannot happen is a DMU with no weight at all, because then no direction
    # can be built for it.
    row_totals <- rowSums(abs(relative_importance))
    empty_rows <- which(row_totals == 0)

    if (length(empty_rows) > 0L) {
      stop(
        sprintf("'directional_vector$relative_importance' has %d row(s) whose weights are all zero, so no direction can be built for those DMUs; first ones: %s.",
                length(empty_rows), paste(head(empty_rows, 5L), collapse = ", ")),
        call. = FALSE
      )
    }

    # PEAXAI_local_importance() normalises |weights| to 1 within each row.
    unnormalised <- abs(row_totals - 1) > 1e-8

    if (any(unnormalised)) {
      warning(
        sprintf("the absolute weights of %d row(s) of 'directional_vector$relative_importance' do not add up to 1; proceeding without normalization; first ones: %s.",
                sum(unnormalised), paste(head(which(unnormalised), 5L), collapse = ", ")),
        call. = FALSE
      )
    }
  }

  # --- directional_vector$baseline --------------------------------------------
  baseline_available <- c("mean", "median", "self", "ones")
  baseline           <- directional_vector[["baseline"]]

  if (is.null(baseline)) {
    stop(
      sprintf("'directional_vector$baseline' is missing. Available: %s.",
              paste(sprintf("\"%s\"", baseline_available), collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.character(baseline) || length(baseline) != 1L || is.na(baseline) ||
      !(baseline %in% baseline_available)) {
    stop(
      sprintf("'directional_vector$baseline' must be one of: %s.",
              paste(sprintf("\"%s\"", baseline_available), collapse = ", ")),
      call. = FALSE
    )
  }

  # ============================================================================
  # 7. n_expand
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_counterfactuals().
  # ============================================================================

  peaxai_require_arg(
    n_expand, "n_expand",
    "It controls how the initial search interval is expanded, e.g. n_expand = 0.5."
  )

  if (!is.numeric(n_expand) || length(n_expand) != 1L || is.na(n_expand) ||
      !is.finite(n_expand) || n_expand < 0) {
    stop("'n_expand' must be a single finite non-negative numeric value.", call. = FALSE)
  }

  # ============================================================================
  # 8. n_grid
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_counterfactuals().
  # ============================================================================

  peaxai_require_arg(
    n_grid, "n_grid",
    "It is the number of candidate counterfactuals evaluated at each iteration, e.g. n_grid = 50."
  )

  if (!is.numeric(n_grid) || length(n_grid) != 1L || is.na(n_grid) ||
      !is.finite(n_grid) || n_grid != floor(n_grid) || n_grid < 2) {
    stop("'n_grid' must be a single whole number >= 2.", call. = FALSE)
  }

  # ============================================================================
  # 9. max_y
  # ============================================================================

  if (!is.numeric(max_y) || length(max_y) != 1L || is.na(max_y) ||
      !is.finite(max_y) || max_y < 0) {
    stop("'max_y' must be a single finite non-negative numeric value.", call. = FALSE)
  }

  # ============================================================================
  # 10. min_x
  # ============================================================================

  if (!is.numeric(min_x) || length(min_x) != 1L || is.na(min_x) ||
      !is.finite(min_x) || min_x < 0) {
    stop("'min_x' must be a single finite non-negative numeric value.", call. = FALSE)
  }

  invisible(NULL)
}


#' @title Validate PEAXAI_ranking input arguments
#'
#' @description
#' Checks the arguments of \code{PEAXAI_ranking()}, in the order of the
#' signature published in Section 2.2.3 of the SoftwareX manuscript (l. 291-292).
#'
#' @param data A \code{data.frame} or \code{matrix} with the DMUs under analysis.
#' @param x Integer vector with the column indices of the input variables.
#' @param y Integer vector with the column indices of the output variables.
#' @param final_model The fitted classifier, of class \code{"train"}.
#' @param efficiency_thresholds Numeric vector of probability levels, strictly
#'   in (0,1). Only required when \code{rank_basis = "attainable"}.
#' @param targets Counterfactual projections returned by
#'   \code{PEAXAI_counterfactuals()}. Required when
#'   \code{rank_basis = "attainable"}.
#' @param rank_basis Either \code{"predicted"} or \code{"attainable"}.
#'
#' @details
#' Manuscript, l. 293-298: with \code{rank_basis = "predicted"} the function
#' ranks units solely by their predicted probability, so neither
#' \code{efficiency_thresholds} nor \code{targets} is needed. With
#' \code{rank_basis = "attainable"} both become mandatory. \code{rank_basis} is
#' therefore validated FIRST, because it decides what the rest must look like.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_ranking <- function(
    data,
    x,
    y,
    final_model,
    efficiency_thresholds,
    targets = NULL,
    rank_basis
) {

  # ============================================================================
  # 1. data
  # ============================================================================

  peaxai_require_arg(
    data, "data",
    "It must be a data.frame or a matrix whose rows are DMUs."
  )

  data_dim <- peaxai_check_data(data, "data")
  n_col    <- data_dim[["n_col"]]

  # ============================================================================
  # 2. x  (input variable indices)
  # ============================================================================

  peaxai_require_arg(
    x, "x",
    "It must be a vector with the column indices of the input variables, e.g. x = 1:3."
  )

  peaxai_check_index(x, "x", "input", n_col)

  # ============================================================================
  # 3. y  (output variable indices)
  # ============================================================================

  peaxai_require_arg(
    y, "y",
    "It must be a vector with the column indices of the output variables, e.g. y = 4."
  )

  peaxai_check_index(y, "y", "output", n_col)
  peaxai_check_xy_overlap(x, y)

  # ============================================================================
  # 4. final_model
  # ============================================================================

  peaxai_check_final_model(final_model, "final_model")

  # Every variable the model was trained on must be present in 'data', matched
  # by its exact name: a column caret renamed at training time must already
  # carry the sanitised name here.
  peaxai_check_model_variables(final_model, data, x, y, "data")

  # ============================================================================
  # 5. rank_basis
  # ----------------------------------------------------------------------------
  # Validated BEFORE 'efficiency_thresholds' and 'targets' because it determines
  # whether those two are required. No default in PEAXAI_ranking().
  # ============================================================================

  rank_basis_available <- c("predicted", "attainable")

  peaxai_require_arg(
    rank_basis, "rank_basis",
    "It must be \"predicted\" or \"attainable\"."
  )

  if (is.null(rank_basis)) {
    stop(
      sprintf("'rank_basis' cannot be NULL. Available: %s.",
              paste(sprintf("\"%s\"", rank_basis_available), collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.character(rank_basis) || length(rank_basis) != 1L || is.na(rank_basis) ||
      !(rank_basis %in% rank_basis_available)) {
    stop(
      sprintf("'rank_basis' must be one of: %s.",
              paste(sprintf("\"%s\"", rank_basis_available), collapse = ", ")),
      call. = FALSE
    )
  }

  # ============================================================================
  # 6. efficiency_thresholds
  # ----------------------------------------------------------------------------
  # Only used when the ranking is threshold-specific.
  # ============================================================================

  if (identical(rank_basis, "attainable")) {
    peaxai_check_thresholds(efficiency_thresholds, "efficiency_thresholds")
  }

  # ============================================================================
  # 7. targets
  # ----------------------------------------------------------------------------
  # Mandatory with rank_basis = "attainable", optional otherwise.
  # ============================================================================

  if (identical(rank_basis, "attainable")) {

    if (is.null(targets)) {
      stop("'targets' is required when 'rank_basis' is \"attainable\". Pass the object returned by PEAXAI_counterfactuals().",
           call. = FALSE)
    }

    peaxai_check_targets(targets, data, x, y, "targets")

  } else if (!is.null(targets)) {

    peaxai_check_targets(targets, data, x, y, "targets")
  }

  invisible(NULL)
}


#' @title Validate PEAXAI_peer input arguments
#'
#' @description
#' Checks the arguments of \code{PEAXAI_peer()}, in the order of the signature
#' published in Section 2.2.3 of the SoftwareX manuscript (l. 284-285).
#'
#' @param data A \code{data.frame} or \code{matrix} with the DMUs under analysis.
#' @param x Integer vector with the column indices of the input variables.
#' @param y Integer vector with the column indices of the output variables.
#' @param final_model The fitted classifier, of class \code{"train"}.
#' @param efficiency_thresholds Numeric vector of probability levels, strictly
#'   in (0,1).
#' @param targets Counterfactual projections returned by
#'   \code{PEAXAI_counterfactuals()}.
#' @param weighted Single logical. If \code{TRUE}, peer distances use the
#'   feature-importance-weighted Euclidean distance.
#' @param relative_importance One-row \code{data.frame} with the weights used
#'   when \code{weighted = TRUE}; ignored otherwise.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_peer <- function(
    data,
    x,
    y,
    final_model,
    efficiency_thresholds,
    targets,
    weighted = FALSE,
    relative_importance = NULL
) {

  # ============================================================================
  # 1. data
  # ============================================================================

  peaxai_require_arg(
    data, "data",
    "It must be a data.frame or a matrix whose rows are DMUs."
  )

  data_dim <- peaxai_check_data(data, "data")
  n_col    <- data_dim[["n_col"]]

  # ============================================================================
  # 2. x  (input variable indices)
  # ============================================================================

  peaxai_require_arg(
    x, "x",
    "It must be a vector with the column indices of the input variables, e.g. x = 1:3."
  )

  peaxai_check_index(x, "x", "input", n_col)

  # ============================================================================
  # 3. y  (output variable indices)
  # ============================================================================

  peaxai_require_arg(
    y, "y",
    "It must be a vector with the column indices of the output variables, e.g. y = 4."
  )

  peaxai_check_index(y, "y", "output", n_col)
  peaxai_check_xy_overlap(x, y)

  # ============================================================================
  # 4. final_model
  # ============================================================================

  peaxai_check_final_model(final_model, "final_model")

  # Every variable the model was trained on must be present in 'data', matched
  # by its exact name: a column caret renamed at training time must already
  # carry the sanitised name here.
  peaxai_check_model_variables(final_model, data, x, y, "data")

  # ============================================================================
  # 5. efficiency_thresholds
  # ============================================================================

  peaxai_check_thresholds(efficiency_thresholds, "efficiency_thresholds")

  # ============================================================================
  # 6. targets
  # ----------------------------------------------------------------------------
  # No default in PEAXAI_peer(). Manuscript, l. 285-289: peers are identified
  # for each counterfactual target, so the targets object is required.
  # ============================================================================

  peaxai_require_arg(
    targets, "targets",
    "It must be the object returned by PEAXAI_counterfactuals()."
  )

  if (is.null(targets)) {
    stop("'targets' cannot be NULL. Pass the object returned by PEAXAI_counterfactuals().",
         call. = FALSE)
  }

  peaxai_check_targets(targets, data, x, y, "targets")

  # ============================================================================
  # 7. weighted
  # ============================================================================

  if (is.null(weighted)) {
    stop("'weighted' cannot be NULL. Use TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(weighted) || length(weighted) != 1L || is.na(weighted)) {
    stop("'weighted' must be a single logical value: TRUE or FALSE.", call. = FALSE)
  }

  # ============================================================================
  # 8. relative_importance
  # ----------------------------------------------------------------------------
  # Only meaningful when weighted = TRUE; ignored otherwise. One weight per
  # variable in c(x, y).
  # ============================================================================

  if (isTRUE(weighted)) {

    if (is.null(relative_importance)) {
      stop("'relative_importance' is required when 'weighted = TRUE'. It is the object returned by PEAXAI_global_importance().",
           call. = FALSE)
    }

    if (!is.data.frame(relative_importance) && !is.matrix(relative_importance) &&
        !is.numeric(relative_importance)) {
      stop(
        sprintf("'relative_importance' must be a one-row data.frame, a matrix or a numeric vector; got an object of class '%s'.",
                paste(class(relative_importance), collapse = "/")),
        call. = FALSE
      )
    }

    if ((is.data.frame(relative_importance) || is.matrix(relative_importance)) &&
        nrow(relative_importance) != 1L) {
      stop(
        sprintf("'relative_importance' must have exactly one row (global weights only); got %d.",
                nrow(relative_importance)),
        call. = FALSE
      )
    }

    variable_names <- colnames(as.data.frame(data))[c(x, y)]
    weight_names   <- colnames(relative_importance)

    if (!is.null(weight_names) && any(nzchar(weight_names))) {

      missing_names <- setdiff(variable_names, weight_names)
      unknown_names <- setdiff(weight_names, variable_names)

      if (length(unknown_names) > 0L) {
        stop(
          sprintf("'relative_importance' has columns that are not inputs or outputs in 'data': %s.",
                  paste(unknown_names, collapse = ", ")),
          call. = FALSE
        )
      }

      if (length(missing_names) > 0L) {
        stop(
          sprintf("'relative_importance' must provide a weight for every input and output in colnames(data)[c(x, y)]; missing: %s.",
                  paste(missing_names, collapse = ", ")),
          call. = FALSE
        )
      }

    } else {

      # Unnamed: matched by position.
      n_weights <- if (is.null(dim(relative_importance))) {
        length(relative_importance)
      } else {
        ncol(relative_importance)
      }

      if (n_weights != length(c(x, y))) {
        stop(
          sprintf("Unnamed 'relative_importance' must have one weight per variable: expected %d (length(c(x, y))), got %d.",
                  length(c(x, y)), n_weights),
          call. = FALSE
        )
      }
    }

    weights <- as.numeric(as.matrix(relative_importance))

    if (anyNA(weights) || any(!is.finite(weights))) {
      stop("'relative_importance' contains NA or non-finite weights.", call. = FALSE)
    }

    if (any(weights < 0)) {
      stop(
        sprintf("'relative_importance' weights must be non-negative; offending values: %s.",
                paste(weights[weights < 0], collapse = ", ")),
        call. = FALSE
      )
    }

    if (all(weights == 0)) {
      stop("'relative_importance' cannot be the all-zero vector.", call. = FALSE)
    }

    if (abs(sum(weights) - 1) > 1e-8) {
      warning("'relative_importance' does not sum to 1; proceeding without normalization.",
              call. = FALSE)
    }
  }

  invisible(NULL)
}


#' @title Validate PEAXAI_predict input arguments
#'
#' @description
#' Checks the arguments of \code{PEAXAI_predict()}, which returns the predicted
#' probability of belonging to the \code{'efficient'} class for each DMU.
#'
#' @param data A \code{data.frame} or \code{matrix} with the DMUs to be scored,
#'   already restricted to the variables the model was trained on.
#' @param final_model The fitted classifier, of class \code{"train"}.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_PEAXAI_predict <- function(
    data,
    final_model
) {

  # ============================================================================
  # 1. data
  # ============================================================================

  peaxai_require_arg(
    data, "data",
    "It must be a data.frame or a matrix whose rows are the DMUs to be scored."
  )

  peaxai_check_data(data, "data")

  # ============================================================================
  # 2. final_model
  # ============================================================================

  peaxai_check_final_model(final_model, "final_model")

  # ============================================================================
  # 3. Cross-argument consistency
  # ----------------------------------------------------------------------------
  # Every variable the model was trained on must be present in 'data', matched
  # by its exact name. Here there are no 'x' and 'y', so every column of 'data'
  # is a candidate; the ones the model does not use are simply not selected.
  # The returned positions are what PEAXAI_predict() uses to subset and reorder.
  # ============================================================================

  model_positions <- peaxai_check_model_variables(
    final_model,
    data,
    arg_name = "data"
  )

  invisible(model_positions)
}


#' @title Validate reffcy input arguments
#'
#' @description
#' Checks the arguments of \code{reffcy()}, the exported generator of simulated
#' samples used in the Supplementary Material of the SoftwareX article.
#'
#' @param DGP Character string naming the data-generating process. Currently
#'   only \code{"cobb_douglas_XnY1"} is implemented.
#' @param parms Named list with the parameters of the data-generating process:
#'   \code{N}, \code{nX} and, optionally, \code{RTS}, \code{min_x}, \code{max_x},
#'   \code{sd_u} and \code{sd_v}.
#'
#' @return Invisibly \code{NULL}.
#'
#' @keywords internal
#' @noRd

validate_parametes_reffcy <- function(
    DGP,
    parms
) {

  # ============================================================================
  # 1. DGP
  # ============================================================================

  DGP_available <- "cobb_douglas_XnY1"

  peaxai_require_arg(
    DGP, "DGP",
    "It must name the data-generating process, e.g. DGP = \"cobb_douglas_XnY1\"."
  )

  if (is.null(DGP)) {
    stop(
      sprintf("'DGP' cannot be NULL. Available: %s.",
              paste(sprintf("\"%s\"", DGP_available), collapse = ", ")),
      call. = FALSE
    )
  }

  if (!is.character(DGP) || length(DGP) != 1L || is.na(DGP)) {
    stop("'DGP' must be a single character string.", call. = FALSE)
  }

  if (!(DGP %in% DGP_available)) {
    stop(
      sprintf("'DGP' = \"%s\" is not available. Available: %s.",
              DGP, paste(sprintf("\"%s\"", DGP_available), collapse = ", ")),
      call. = FALSE
    )
  }

  # ============================================================================
  # 2. parms
  # ============================================================================

  peaxai_require_arg(
    parms, "parms",
    "It must be a named list, e.g. list(N = 100, nX = 3, RTS = \"vrs\")."
  )

  if (!is.list(parms) || is.data.frame(parms)) {
    stop(
      sprintf("'parms' must be a named list; got an object of class '%s'.",
              paste(class(parms), collapse = "/")),
      call. = FALSE
    )
  }

  parms_available <- c("N", "nX", "RTS", "min_x", "max_x", "sd_u", "sd_v")
  parms_names     <- names(parms)

  if (is.null(parms_names) || any(!nzchar(parms_names))) {
    stop(
      sprintf("'parms' must be a NAMED list. Accepted names: %s.",
              paste(parms_available, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!all(parms_names %in% parms_available)) {
    stop(
      sprintf("'parms' contains unknown parameters: %s. Accepted: %s.",
              paste(setdiff(parms_names, parms_available), collapse = ", "),
              paste(parms_available, collapse = ", ")),
      call. = FALSE
    )
  }

  # --- parms$N (sample size) --------------------------------------------------
  if (is.null(parms[["N"]])) {
    stop("'parms$N' is missing. It is the sample size, e.g. N = 100.", call. = FALSE)
  }

  N <- parms[["N"]]

  if (!is.numeric(N) || length(N) != 1L || is.na(N) || !is.finite(N) ||
      N != floor(N) || N < 1) {
    stop("'parms$N' must be a single whole number >= 1.", call. = FALSE)
  }

  # --- parms$nX (number of inputs) --------------------------------------------
  if (is.null(parms[["nX"]])) {
    stop("'parms$nX' is missing. It is the number of inputs, e.g. nX = 3.", call. = FALSE)
  }

  nX           <- parms[["nX"]]
  nX_available <- c(1, 3, 6, 9, 12, 15)

  if (!is.numeric(nX) || length(nX) != 1L || is.na(nX) || !(nX %in% nX_available)) {
    stop(
      sprintf("'parms$nX' must be one of: %s.", paste(nX_available, collapse = ", ")),
      call. = FALSE
    )
  }

  # --- parms$RTS --------------------------------------------------------------
  if (!is.null(parms[["RTS"]])) {

    RTS <- parms[["RTS"]]

    if (!is.character(RTS) || length(RTS) != 1L || is.na(RTS) ||
        !(RTS %in% c("vrs", "crs"))) {
      stop("'parms$RTS' must be \"vrs\" or \"crs\".", call. = FALSE)
    }
  }

  # --- parms$min_x and parms$max_x (input bounds) -----------------------------
  for (bound_name in c("min_x", "max_x")) {

    if (!is.null(parms[[bound_name]])) {

      bound <- parms[[bound_name]]

      if (!is.numeric(bound) || length(bound) != 1L || is.na(bound) || !is.finite(bound)) {
        stop(
          sprintf("'parms$%s' must be a single finite numeric value.", bound_name),
          call. = FALSE
        )
      }
    }
  }

  if (!is.null(parms[["min_x"]]) && !is.null(parms[["max_x"]])) {
    if (parms[["min_x"]] >= parms[["max_x"]]) {
      stop("'parms$min_x' must be strictly smaller than 'parms$max_x'.", call. = FALSE)
    }
  }

  # --- parms$sd_u and parms$sd_v (standard deviations) ------------------------
  for (sd_name in c("sd_u", "sd_v")) {

    if (!is.null(parms[[sd_name]])) {

      sd_value <- parms[[sd_name]]

      if (!is.numeric(sd_value) || length(sd_value) != 1L || is.na(sd_value) ||
          !is.finite(sd_value) || sd_value < 0) {
        stop(
          sprintf("'parms$%s' must be a single finite non-negative numeric value.", sd_name),
          call. = FALSE
        )
      }
    }
  }

  invisible(NULL)
}


#' @title Prepare Data and Handle Errors
#'
#' @description This function arranges the data in the required format and displays some error messages.
#'
#' @param data A \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Column indexes of input variables in \code{data}.
#' @param y Column indexes of output variables in \code{data}.
#'
#' @importFrom caret trainControl
#'
#' @return It returns a \code{matrix} in the required format and displays some error messages.

preprocessing <- function (
    data, x, y
) {

  # x and y well / bad introduced

  cols <- 1:length(data)
  if (!(all(x %in% cols))) {
    stop("x index(es) are not in data.")

    if (!(all(y %in% cols))) {
      stop("y index(es) are not in data.")
    }
  }

  # (i) data.frame, (ii) list with variables, (iii) matrix

  # data.frame format to deal with classes
  if (is.list(data) && !is.data.frame(data)) {

    # data names?
    ifelse(is.null(names(data)),
           var_names <- 1:length(data), # if not 1:x
           var_names <- names(data)
    )

    data <- matrix(unlist(data), ncol = length(var_names), byrow = F)
    colnames(data) <- var_names

  } else if (is.matrix(data) || is.data.frame(data)) {
    data <- data.frame(data)
  }

  # Classes
  varClass <- unlist(sapply(data, class))

  # Output classes
  outClass <- varClass[y] %in% c("numeric", "double", "integer")

  # Error
  if (!all(outClass)){
    stop(paste(names(data)[y][!outClass][1], "is not a numeric or integer vector"))
  }

  # Input classes
  # Ordered --> numeric
  for (i in x){
    if (is.ordered(data[, i])) {
      data[, i] <- as.numeric(data[, i])
    }
  }

  # Define classes again
  varClass <- unlist(sapply(data, class))

  inpClass <- varClass[x] %in% c("numeric", "double", "integer")

  # Error
  if (!all(inpClass)){
    stop(paste(names(data)[x][!inpClass][1], "is not a numeric, integer or ordered vector"))
  }

  data <- data[, c(x, y)]

  return(as.matrix(data))
}
