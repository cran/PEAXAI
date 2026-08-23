#' @title Random Sample for Efficiency Analysis
#'
#' @description This function is used to simulate a \code{data.frame} with some given Data Generation Process.
#'
#' @param DGP Data Generation Process:
#' \itemize{
#'    \item{\code{"cobb_douglas_XnY1"}}: Cobb-Douglas Data Generation Process for a XnY1 scenario. Check \code{help("cobb_douglas_XnY1")}.
#'  }
#' @param parms \code{list} with the parameters of the simulation functions.
#' \itemize{
#'    \item{\code{"N"}}: Sample size.
#'    \item{\code{"nX"}}: Number of inputs. Available for the \code{cobb_douglas_XnY1(N, nX, RTS)} function.
#'    \item{\code{"RTS"}}: Returns to scale, \code{"vrs"} (default) or \code{"crs"}. Available for the \code{cobb_douglas_XnY1(N, nX, RTS)} function.
#'    \item{\code{"min_x"}}: Lower bound of the (uniformly distributed) inputs (default: \code{1}).
#'    \item{\code{"max_x"}}: Upper bound of the (uniformly distributed) inputs (default: \code{10}).
#'    \item{\code{"sd_u"}}: Standard deviation of the inefficiency term (default: \code{0.4}).
#'    \item{\code{"sd_v"}}: Standard deviation of the noise term (default: \code{0}, i.e. a deterministic frontier). Set to a value greater than \code{0} to obtain a stochastic frontier.
#'  }
#'
#' @details
#' Please refer to the help manuals of the mentioned functions for usage examples and more detailed information on the parameters.
#'
#'
#' @return \code{data.frame} simulated with the selected Data Generation Process.
#'
#' @examples
#' set.seed(1)
#'
#' # 100 DMUs with 3 inputs and 1 output under variable returns to scale
#' simulated_data <- reffcy(
#'   DGP   = "cobb_douglas_XnY1",
#'   parms = list(
#'     N   = 100,
#'     nX  = 3,
#'     RTS = "vrs"
#'   )
#' )
#'
#' # x1 ... x3 are the inputs, 'y' the observed output and 'yD' the
#' # theoretical frontier
#' head(simulated_data)
#'
#' # A stochastic frontier is obtained by setting sd_v > 0
#' noisy_data <- reffcy(
#'   DGP   = "cobb_douglas_XnY1",
#'   parms = list(N = 100, nX = 1, RTS = "crs", sd_u = 0.4, sd_v = 0.1)
#' )
#'
#' str(noisy_data)
#'
#' @export
reffcy <- function (
    DGP, parms
    ) {

  # check if parameters are well introduced
  validate_parametes_reffcy(
    DGP   = DGP,
    parms = parms
  )

  data <- do.call(cobb_douglas_XnY1, parms)

  return(data)

}

#' @title 1 output ~ nX inputs Cobb-Douglas Data Generation Process
#'
#' @description This function is used to simulate a 1 output ~ nX inputs \code{data.frame} with a Cobb-Douglas Data Generation Process.
#'
#' @param N Sample size.
#' @param nX Number of inputs. Possible values: \code{1}, \code{3}, \code{6}, \code{9}, \code{12} and \code{15}.
#' @param RTS Returns to scale assumed for the theoretical frontier (default: \code{"vrs"}). Accepted values:
#'   \describe{
#'     \item{\code{"vrs"}}{Variable returns to scale. The output elasticities sum to \code{0.5}.}
#'     \item{\code{"crs"}}{Constant returns to scale. The output elasticities are doubled with respect to \code{"vrs"} so they sum to \code{1}.}
#'   }
#' @param min_x,max_x Lower and upper bound of the uniformly distributed inputs (default: \code{1} and \code{10}).
#' @param sd_u Standard deviation of the (half-normal) inefficiency term (default: \code{0.4}).
#' @param sd_v Standard deviation of the (normal) noise term (default: \code{0}). With the default, the frontier
#'   is deterministic (\code{y = yD * exp(-u)}); setting \code{sd_v > 0} yields a stochastic frontier
#'   (\code{y = yD * exp(v - u)}).
#'
#' @importFrom stats runif rnorm
#'
#' @return \code{data.frame} with the simulated data: nX inputs, 1 output (y) and the theoretical frontier (yD).
#'
#' @keywords internal

cobb_douglas_XnY1 <- function (
    N, nX, RTS = "vrs",
    min_x = 1, max_x = 10,
    sd_u = 0.4, sd_v = 0
    ) {

  if (!(nX %in% c(1, 3, 6, 9, 12, 15))) {
    stop(paste(nX, "is not allowed"))
  }

  if (!(RTS %in% c("vrs", "crs"))) {
    stop(paste(RTS, "is not allowed"))
  }

  # Output elasticities under VRS (they sum to 0.5 in every scenario)
  elasticities <- list(
    "1"  = c(0.5),
    "3"  = c(0.05, 0.15, 0.3),
    "6"  = c(0.05, 0.001, 0.004, 0.045, 0.1, 0.3),
    "9"  = c(0.005, 0.001, 0.004, 0.005, 0.001, 0.004, 0.08, 0.1, 0.3),
    "12" = c(0.005, 0.001, 0.004, 0.005, 0.001, 0.004, 0.08, 0.05, 0.05, 0.075, 0.025, 0.2),
    "15" = c(0.005, 0.001, 0.004, 0.005, 0.001, 0.004, 0.08, 0.05, 0.05, 0.05, 0.025, 0.025, 0.025, 0.025, 0.15)
  )[[as.character(nX)]]

  # Under CRS the elasticities are doubled so they sum to 1
  if (RTS == "crs") {
    elasticities <- elasticities * 2
  }

  colnames <- c(paste("x", 1:nX, sep = ""), "y")

  data <- as.data.frame (
    matrix (
      ncol = length(colnames),
      nrow = N,
      dimnames = list(NULL, colnames)
      )
    )

  # Input generation
  for (x in 1:nX){
    data[, x] <- runif(n = N, min = min_x, max = max_x)
  }

  # Theoretical frontier: yD = 3 * prod(x_i ^ elasticity_i)
  log_inputs <- log(as.matrix(data[, 1:nX, drop = FALSE]))
  data[, "yD"] <- 3 * exp(as.vector(log_inputs %*% elasticities))

  # Inefficiency generation
  u <- abs(rnorm(n = N, mean = 0, sd = sd_u))

  # Noise generation (sd_v = 0 by default, i.e. no noise / deterministic frontier)
  v <- rnorm(n = N, mean = 0, sd = sd_v)

  # Output generation
  data[, "y"] <- data[, "yD"] * exp(v - u)

  return(data)
}
