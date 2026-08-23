#' @title Create New SMOTE Units to Balance Data combinations of m + s
#'
#' @description This function creates new DMUs to address data imbalances. If the majority class is efficient, it generates new inefficient DMUs by worsering the observed units. Conversely, if the majority class is inefficient, it projects inefficient DMUs to the frontier. Finally, a random selection if performed to keep a proportion of 0.65 for the majority class and 0.35 for the minority class.
#'
#' @param data A \code{data.frame} containing the variables used in the model.
#' @param x Column indexes of the input variables in the \code{data}.
#' @param y Column indexes of the output variables in the \code{data}.
#' @param RTS Text string or number defining the underlying DEA technology /
#'   returns-to-scale assumption (default: \code{"vrs"}). Accepted values:
#'   \describe{
#'     \item{\code{1} / \code{"vrs"}}{Variable returns to scale, convexity and free disposability.}
#'     \item{\code{3} / \code{"crs"}}{Constant returns to scale, convexity and free disposability.}
#'   }
#' @param balance_data Indicate level of efficient units to achive and the number of efficient and not efficient units.
#' @param seed  Integer. Seed for reproducibility.
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default \code{FALSE}).
#'
#' @return It returns a \code{data.frame} with the newly created set of DMUs incorporated.
#'
#' @keywords internal

SMOTE_data <- function (
    data, x, y, RTS = "vrs", balance_data, seed, verbose = TRUE
) {

  # first, determine the efficient facets
  facets <- convex_facets(
    data = data,
    x = x,
    y = y,
    RTS = RTS,
    verbose = verbose
  )

  # second, populate the efficient facets
  balance_datasets <- get_SMOTE_DMUs(
    data = data,
    facets = facets,
    x = x,
    y = y,
    balance_data = balance_data,
    seed = seed,
  )

  return(balance_datasets)
}
