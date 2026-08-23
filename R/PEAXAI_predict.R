#' @title Predict Probability of Efficiency Using a Fitted Model
#'
#' @description
#' Predicts probabilities for new decision-making units (DMUs) using a fitted
#' \pkg{caret} classification model. If \code{calibration_model} is provided,
#' the raw classifier probabilities are post-processed to obtain calibrated
#' probability estimates (e.g., Platt scaling via logistic regression or
#' isotonic calibration via a monotone mapping).
#'
#' @param data A \code{data.frame} or \code{matrix} with the DMUs to be scored. It may
#'   carry any number of additional columns: only the variables \code{final_model} was
#'   trained on are kept, in the order the model expects. Those variables are matched by
#'   their \strong{exact} name, so a variable \pkg{caret} renamed when the model was
#'   trained (\code{"Total assets"} becomes \code{"Total.assets"}) must already carry the
#'   sanitised name in \code{data}; otherwise it is reported as missing.
#' @param final_model A fitted \pkg{caret} model provided by PEAXAI_fitting().
#'
#' @return A numeric vector of predicted probabilities (raw or calibrated), one per
#'   row of \code{data}.
#'
#' @examples
#' \donttest{
#'   data("firms", package = "PEAXAI")
#'
#'   # Firms of the Valencian Community: four inputs and one output
#'   data <- subset(
#'     firms,
#'     autonomous_community == "Comunidad Valenciana"
#'   )
#'
#'   x <- 1:4
#'   y <- 5
#'
#'   models <- PEAXAI_fitting(
#'     data            = data,
#'     x               = x,
#'     y               = y,
#'     RTS             = "vrs",
#'     imbalance_rate  = NULL,
#'     methods         = list("glm" = list(weights = "dinamic")),
#'     trControl       = list(method = "cv", number = 3),
#'     metric_priority = c("Balanced_Accuracy", "ROC_AUC"),
#'     seed            = 1,
#'     verbose         = FALSE
#'   )
#'
#'   final_model <- models[["best_model_fit"]][["glm"]]
#'
#'   # 'data' must contain every variable the model was trained on
#'   probabilities <- PEAXAI_predict(
#'     data        = data[, c(x, y)],
#'     final_model = final_model
#'   )
#'
#'   head(probabilities)
#' }
#'
#' @export

PEAXAI_predict <- function (
    data, final_model
) {

  # check if parameters are well introduced
  model_positions <- validate_parametes_PEAXAI_predict(
    data        = data,
    final_model = final_model
  )

  data <- as.data.frame(data)

  # Keep only the variables final_model was trained on, in the order it expects.
  # Any other column of 'data' is dropped here. The names are already exact:
  # peaxai_check_model_variables() matched them character by character and
  # stopped otherwise, so nothing has to be renamed at this point.
  if (!is.null(model_positions)) {

    data <- data[, unname(model_positions), drop = FALSE]

    # Post-condition. It can only fail if model_positions is built wrongly, so
    # it is an assertion on PEAXAI itself, not on what the user supplied.
    if (!identical(names(data), names(model_positions))) {
      stop(
        sprintf(
          paste0(
            "Internal error: the columns of 'data' could not be aligned with the ",
            "variables 'final_model' was trained on. After the selection 'data' has %s ",
            "but %s was expected."
          ),
          paste(sprintf("'%s'", names(data)), collapse = ", "),
          paste(sprintf("'%s'", names(model_positions)), collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  name_model <- final_model$method

  calibration_model <- NULL
  if(is.null(calibration_model)) {

    # --------------------------------------------------------------------------
    # No calibration -----------------------------------------------------------
    # --------------------------------------------------------------------------
    predictions <- predict(final_model, newdata = data, type = "prob")[,1]

  } else {

    # --------------------------------------------------------------------------
    # Calibration --------------------------------------------------------------
    # --------------------------------------------------------------------------

    # Original probabilities
    s_new <- predict(final_model,
                     newdata = data,
                     type = "prob")[,1]

    if (calibration_model$method == "glm.fit") {

      # platt
      predictions <- predict(calibration_model,
                             newdata = data.frame(s = s_new),
                             type = "response")
    } else {

      # isotonic
      predictions <- calibration_model[[name_model]](s_new)

    }

  }

  return(predictions)
}
