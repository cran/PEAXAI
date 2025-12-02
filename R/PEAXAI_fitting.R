#' @title Training Classification Models to Estimate Efficiency
#'
#' @description
#' Trains one or multiple classification algorithms to identify Pareto-efficient
#' decision-making units (DMUs). It jointly searches model hyperparameters and the
#' class-balancing level (synthetic samples via SMOTE) using k-fold cross-
#' validation or a train/validation/test split, selecting the configuration that
#' maximizes the specified metric(s). Returns, for each technique, the best fitted
#' model together with training summaries, performance metrics, and the selected
#' balancing level.
#'
#' @param data A \code{data.frame} or \code{matrix} containing the variables in the model.
#' @param x Integer vector with column indices of input variables in \code{data}.
#' @param y Integer vector with column indices of output variables in \code{data}.
#' @param RTS Text string or number defining the underlying DEA technology /
#'   returns-to-scale assumption (default: \code{"vrs"}). Accepted values:
#'   \describe{
#'     \item{\code{0} / \code{"fdh"}}{Free disposability hull, no convexity assumption.}
#'     \item{\code{1} / \code{"vrs"}}{Variable returns to scale, convexity and free disposability.}
#'     \item{\code{2} / \code{"drs"}}{Decreasing returns to scale, convexity, down-scaling and free disposability.}
#'     \item{\code{3} / \code{"crs"}}{Constant returns to scale, convexity and free disposability.}
#'     \item{\code{4} / \code{"irs"}}{Increasing returns to scale (up-scaling, not down-scaling), convexity and free disposability.}
#'     \item{\code{5} / \code{"add"}}{Additivity (scaling up and down, but only with integers), and free disposability.}
#'   }
#' @param trControl A \code{caret::trainControl}-like list that specifies the resampling
#'   strategy; recognized values for \code{$method} include \code{"cv"}, \code{"test_set"},
#'   and \code{"none"}. See \pkg{caret} documentation.
#' @param methods A \code{list} of selected machine learning models and their hyperparameters.
#' @param metric_priority A \code{string} specifying the summary metric for classification to select the optimal model. Default includes \code{"Balanced_Accuracy"} due to (normally) unbalanced data.
#' @param hold_out Numeric proportion in (0,1) for validation split (default \code{NULL}).
#'   If \code{NULL}, training and validation use the same data.
#' @param imbalance_rate Optional target(s) for class balance via SMOTE. If \code{NULL},
#'   no synthetic balancing is performed.
#' @param seed  Integer. Seed for reproducibility.
#' @param verbose Logical; if \code{TRUE}, prints progress messages (default \code{FALSE}).
#'
#' @importFrom caret train createDataPartition confusionMatrix createFolds
#' @importFrom dplyr %>% arrange across filter select summarise group_by first semi_join all_of
#' @importFrom pROC roc
#' @importFrom PRROC pr.curve
#' @importFrom stats binomial family median na.omit runif setNames
#' @importFrom utils combn tail
#'
#' @return A \code{"PEAXAI"} (list) with the best technique, best fitted models and their performance and the results by fold.
#'
#' @examples
#' \donttest{
#'   data("firms", package = "PEAXAI")
#'
#'   data <- subset(
#'     firms,
#'     autonomous_community == "Comunidad Valenciana"
#'   )
#'
#'   trControl <- list(
#'     method = "cv",
#'     number = 3
#'   )
#'
#'   # glm method
#'   methods <- list(
#'     "glm" = list(
#'         weights = "dinamic"
#'      )
#'   )
#'
#'   models <- PEAXAI_fitting(
#'     data = data,
#'     x = c(1:4),
#'     y = 5,
#'     RTS = "vrs",
#'     imbalance_rate = NULL,
#'     methods = methods,
#'     trControl = trControl,
#'     metric_priority = c("Balanced_Accuracy", "ROC_AUC"),
#'     seed = 1,
#'     verbose = FALSE
#'   )
#' }
#'
#' @export

PEAXAI_fitting <- function (
    data, x, y, RTS = "vrs", imbalance_rate = NULL,
    trControl, methods, metric_priority = "Balanced_Accuracy",
    hold_out = NULL, seed = NULL, verbose = TRUE
    ) {

  # ----------------------------------------------------------------------------
  # pre-processing -------------------------------------------------------------
  # ----------------------------------------------------------------------------
  # check if parameters are well introduced
  validate_parametes_PEAXAI_fitting(
    data = data,
    x = x,
    y = y,
    RTS = RTS,
    trControl = trControl,
    methods = methods,
    metric_priority = metric_priority,
    imbalance_rate = imbalance_rate,
    hold_out = hold_out,
    verbose = verbose
  )

  if (is.null(seed)) {
    seed <- 314
  }

  # reorder index 'x' and 'y' in data
  data <- data[, c(x,y)]
  x <- 1:(ncol(data) - length(y))
  y <- (length(x) + 1):ncol(data)

  # Check data
  data <- preprocessing(
    data = data,
    x = x,
    y = y
  )

  # save a copy before add class_efficiency
  copy_data_no_label <- data

  # ----------------------------------------------------------------------------
  # Step 1: Data labeling ------------------------------------------------------
  # ----------------------------------------------------------------------------
  data <- label_efficiency(
    data = data,
    x = x,
    y = y,
    RTS = RTS
  )

  if (isTRUE(verbose)) {
    info_imb_rate <- round(prop.table(table(data$class_efficiency))[1], 4) * 100
    message(
      paste0("The dataset has an imbalance rate of ", info_imb_rate, "%.")
    )
  }

  # ----------------------------------------------------------------------------
  # Step 2: Create validation set ----------------------------------------------
  # ----------------------------------------------------------------------------
  # create datasets: (train + test) and (validation)
  if (is.null(hold_out)) {

    # all data is use to train and for validation
    valid_data <- data
    train_data <- data
    all_data <- data

  } else {

    # reproducibility.
    set.seed(seed)

    # validation set
    valid_index <- createDataPartition(
     data$class_efficiency,
      p = hold_out,
      list = FALSE)

    # split data
    valid_data <- data[valid_index, c(x,y)]
    train_data <- data[-valid_index, c(x,y)]

    # label DMUs by technology observed in train
    valid_data <- label_efficiency(
      data = valid_data,
      REF = train_data,
      x = x,
      y = y,
      RTS = RTS
    )

    train_data <- label_efficiency(
      data = train_data,
      REF = train_data,
      x = x,
      y = y,
      RTS = RTS
    )

    all_data <- data

    save_performance_validation <- vector("list", length = length(methods))
    names(save_performance_validation) <- names(methods)

  }

  # ----------------------------------------------------------------------------
  # Step 3: ML model training --------------------------------------------------
  # ----------------------------------------------------------------------------
  # the original train dataset
  real_balance <- prop.table(table(train_data$class_efficiency))[["efficient"]]
  datasets_to_train <- list(train_data)
  names(datasets_to_train)[1] <- paste0(as.character(round(real_balance, 4)),"*")
  real_balance_stable <- paste0(as.character(round(real_balance, 4)),"*")

  # ----------------------------------------------------------------------------
  # Step 3.1: Addressing imbalance rate (if there is convexity) ----------------
  # ----------------------------------------------------------------------------
  RTS_available <- c("1", "vrs", "0", "crs")

  if (!is.null(imbalance_rate)) {

    if (as.character(RTS) %in% RTS_available) {

      # determine complete facets
      train_data_SMOTE <- SMOTE_data(
        data = train_data,
        x = x,
        y = y,
        RTS = RTS,
        balance_data = imbalance_rate,
        seed = seed
      )

      datasets_to_train <- append(train_data_SMOTE, datasets_to_train, after = 0)

    } else {
      message(
        "Without the convexity assumption, SMOTE units cannot be created on the efficient frontier."
      )
    }

  }

  # ----------------------------------------------------------------------------
  # Step 3.2: Performance with different hyperparameters -----------------------
  # ----------------------------------------------------------------------------
  # save model
  best_models <- vector("list", length(methods))
  names(best_models) <- names(methods)

  # save performance
  best_performance  <- vector("list", length(methods))
  names(best_performance) <- names(methods)

  # more information about fitting (grid)
  performance_train_all_by_dataset <- vector("list", length(datasets_to_train))
  names(performance_train_all_by_dataset) <- names(datasets_to_train)

  # only the best by imbalance and ML method
  best_performance_train_all_by_dataset <- vector("list", length(datasets_to_train))
  names(best_performance_train_all_by_dataset) <- names(datasets_to_train)

  # metrics evaluation
  metric <- metric_priority[1]

  # save results each dataset
  for (dataset in names(datasets_to_train)) {
    if (isTRUE(verbose)) {
      message(paste0("Case: ", dataset, " imbalance rate"))
    }

    # create k-folds. Same folds for every ML method
    folds <- createFolds(
      datasets_to_train[[dataset]]$class_efficiency,
      k = trControl[["number"]]
    )

    # method
    performance_train_all_by_method <- vector("list", length(methods))
    names(performance_train_all_by_method) <- names(methods)

    best_performance_train_all_by_method <- vector("list", length(methods))
    names(best_performance_train_all_by_method) <- names(methods)

    for (method_i in names(methods)) {
      if (isTRUE(verbose)) {
        message(paste0("Training ", method_i, " method."))
      }

      # train model
      model_fit <- train_PEAXAI(
        data = datasets_to_train[[dataset]],
        method = method_i,
        parameters = methods[[method_i]],
        trControl = trControl,
        metric_priority = metric_priority,
        seed = seed
      )

      results_fit <- model_fit[["results"]]

      cols_out <- names(methods[[method_i]]$tuneGrid)

      results_fit <- results_fit %>%
        select(-all_of(cols_out))

      results_fit <- round(results_fit, 2)

      results_fit <- cbind(model_fit[["results"]][cols_out], results_fit)

      # # nuevo !!!!
      # if (method_i != "glm") {
      #
      #   if(nrow(results_fit > 1) && any(results_fit$cross_entropy < 0.01)) {
      #     results_fit <- results_fit %>%
      #       filter(cross_entropy > 0)
      #   }
      #
      # }

      # # first, eliminate Recall == 0
      # if(nrow(results_fit) > 1 & any(results_fit$Recall > 0.1)) {
      # results_fit <- results_fit %>%
      #   filter(cross_entropy  > 0.1)
      # }

      # e <- 0.01
      # # # second, eliminate Specificity == e
      # # if(nrow(results_fit) > 1 & any(results_fit$Specificity < e)) {
      # # results_fit <- results_fit %>%
      # #   filter(Specificity < e)
      # # }
      #
      # maxSpec <- max(results_fit$Specificity)
      #
      # results_fit <- results_fit %>%
      #   filter(Specificity >= maxSpec - e)

      # save results cv
      performance_train_all_by_method[[method_i]] <- results_fit

# calibrating end

      # results_fit <- results_fit %>%
      #   arrange(cross_entropy, desc(Recall), desc(`ROC-AUC`), cross_entropySD)

#       if (method_i == "glm") {
#
#         best <- model_fit[["results"]]
#
#         if(trControl[["method"]] == "cv") {
#
#           pred_best <- model_fit$pred
#
#           # only observed DMUs; no SMOTE
#           pred_best <- pred_best %>%
#             filter(rowIndex %in% 1:nrow(copy_data_no_label))
#
#         } else if (trControl[["method"]] == "repeatedcv") {
#
#           pred_best <- model_fit$pred %>%
#             group_by(rowIndex) %>%
#             summarise(
#               efficient = mean(efficient),
#               obs  = first(obs),
#               .groups = "drop"
#             ) %>%
#             arrange(rowIndex)
#
#           # only observed DMUs; no SMOTE
#           pred_best <- pred_best[1:nrow(copy_data_no_label),]
#
#         }
#
#       } else {
#
#         if(trControl[["method"]] == "cv") {
#
#           best <- results_fit[1, names(methods[[method_i]]$tuneGrid)]
#
#           best
#           pred_all <- model_fit$pred
#
#           idx_best <- apply(pred_all[, names(best), drop = FALSE], 1, function(z){
#             all(z == best)
#           })
#
#           pred_best <- pred_all[idx_best,]
#
#           # only observed DMUs; no SMOTE
#           pred_best <- pred_best %>%
#             filter(rowIndex %in% 1:nrow(copy_data_no_label))
#
#         } else if (trControl[["method"]] == "repeatedcv") {
#
#           best <- results_fit[1, names(methods[[method_i]]$tuneGrid)]
#           best
#
#           pred_best <- model_fit$pred %>%
#             semi_join(best, by = names(best)) %>%
#             group_by(rowIndex) %>%
#             summarise(
#               efficient = mean(efficient),
#               obs = first(obs),
#               .groups = "drop"
#             ) %>%
#             arrange(rowIndex)
#
#           # only observed DMUs; no SMOTE
#           pred_best <- pred_best[1:nrow(copy_data_no_label),]
#
#         }
#
#       }
#
#       calibrating_datasets_by_method[[method_i]] <- pred_best
# # calibrating

      # order to find the best hyperparameter combination
      idx_order <- do.call(
        order,
        c(
          as.list(model_fit[["results"]][, metric_priority, drop = FALSE]),
          list(decreasing = TRUE, na.last = TRUE)
        )
      )

      results_order <- model_fit[["results"]][idx_order, ]

      best_performance_train_all_by_method[[method_i]] <- results_fit[1,]

    }

    performance_train_all_by_dataset[[dataset]] <- performance_train_all_by_method
    best_performance_train_all_by_dataset[[dataset]] <- best_performance_train_all_by_method

  } # end addressing imbalance rate

  # ----------------------------------------------------------------------------
  # Step 3.3: Train the best models and their performance ----------------------
  # ----------------------------------------------------------------------------
  # get the best in each imbalance rate
  get_model <- function(lst, model = names(methods)) {
    model <- match.arg(model)
    out <- lapply(names(lst), function(k) {
      x <- lst[[k]][[model]]
      if (is.null(x)) return(NULL)
      cbind(model = model, imbalance = k, x, row.names = NULL)
    })
    do.call(rbind, out)
  }

  # train the best models without cross-validation
  save_performance_train <- vector("list", length = length(methods))
  names(save_performance_train) <- names(methods)

  save_best_model_fit <- vector("list", length = length(methods))
  names(save_best_model_fit) <- names(methods)

  for (method_i in names(methods)) {

    # first, select the best imbalance rate
    # dataframe with results of CV of method_i
    best_results_train <- get_model(best_performance_train_all_by_dataset, method_i)
    best_results_train <- best_results_train[,2:ncol(best_results_train)]

    best_results_train <- as.data.frame(best_results_train)

    best_results_train$imbalance <- as.numeric(gsub("\\*", "", best_results_train$imbalance))
    best_results_train$diff_imbalance <- best_results_train$imbalance - round(real_balance, 4)

    # how similar is the dataset to original dataset
    best_results_train$diff_imbalance <- 1 - best_results_train$diff_imbalance

    # add * in the orginal dataset name
    best_results_train$imbalance[best_results_train$diff_imbalance == 1] <- paste0(best_results_train$imbalance[best_results_train$diff_imbalance == 1], "*")

    # add difference imbalance
    metric_priority <- c(metric_priority, "diff_imbalance")

    # arrange(cross_entropy, desc(Recall), desc(`ROC-AUC`), cross_entropySD)
    # # filter rows if has 0 en Recall
    # best_results_train <- best_results_train %>%
    #   filter(Recall > 0)
    #
    # best_results_train <- best_results_train %>%
    #   filter(Specificity < 1)

    best_results_train <- best_results_train[do.call(order, c(best_results_train[metric_priority], list(decreasing = TRUE))), ]

    best_balance_i <- best_results_train[1, "imbalance"]

    # copy of methods
    best_methods <- methods

    # only the best hyperparameters
    names_tuneGrid <- names(best_methods[[method_i]][["tuneGrid"]])

    if (method_i != "glm") {
      best_methods[[method_i]][["tuneGrid"]] <- best_results_train[1, names_tuneGrid]
    }

    metric_priority <- metric_priority[-length(metric_priority)]

    # train the best model
    # NO VALIDATION SET
    # All data used, performance is obteined of cv results
    model_fit <- train_PEAXAI(
      data = datasets_to_train[[as.character(best_balance_i)]],
      method = method_i,
      parameters = best_methods[[method_i]],
      trControl = trainControl(method = "none", classProbs = TRUE),
      metric_priority = metric_priority,
      seed = seed
    )

    performance <- best_results_train[1,]

    # save results
    save_performance_train[[method_i]] <- performance
    save_best_model_fit[[method_i]] <- model_fit


    # WITH VALIDATION SET
    # It is necessary to check the performance of the best model without cv on validation set
    # after that, it is necessary to train again WITH ALL data and the best hyperparameters
    if (!is.null(hold_out)) {

      # check performance on valid data
      y_hat <- predict(model_fit, newdata = valid_data, type = "prob")[["efficient"]]
      y_obs <- valid_data$class_efficiency

      if (any(is.na(y_hat))) {

        out_names <- c("Accuracy", "Kappa",
                       "Recall", "Specificity", "Precision", "F1", "Balanced_Accuracy", "G_mean",
                       "ROC-AUC", "PR-AUC",
                       "Cross_Entropy", "Cross_Entropy_Efficient")

        performance <- setNames(rep(NA, length(out_names)), out_names)

      } else {

        # levels
        levls <- c("efficient", "not_efficient")

        # predictions
        y_hat_prob <- y_hat
        y_hat <- ifelse(y_hat > 0.5, "efficient", "not_efficient")
        y_hat <- factor(y_hat, levels = levls)

        # reference
        y_obs <- factor(y_obs, levels = levls)

        # confusion matrix
        cm <- confusionMatrix(
          data = y_hat,
          reference = y_obs,
          mode = "everything",
          positive = "efficient"
        )

        # G-mean
        gmean <- sqrt(cm$byClass["Recall"]*cm$byClass["Specificity"])

        # ROC-AUC
        roc_obj <- roc(
          response = y_obs,
          predictor = y_hat_prob,
          levels = rev(levls),
          direction = "<",
          quiet = TRUE)

        # PR-AUC
        pr_obj <- pr.curve(
          scores.class0 = y_hat_prob[y_obs == "not_efficient"],
          scores.class1 = y_hat_prob[y_obs == "efficient"],
          curve = TRUE
        )

        # cross_entropy
        y_bin <- ifelse(y_obs == "efficient", 1, 0)
        p <- y_hat_prob

        # no log(0)
        eps <- 1e-15
        p_clipped <- pmin(pmax(p, eps), 1 - eps)

        cross_entropy <- -mean(
          y_bin * log(p_clipped) + (1 - y_bin) * log(1 - p_clipped)
        )

        # cross entropy efficient class
        p <- y_hat_prob[y_bin == 1]
        y_bin <- y_bin[y_bin == 1]

        # no log(0)
        eps <- 1e-15
        p_clipped <- pmin(pmax(p, eps), 1 - eps)

        cross_entropy_efficient <- -mean(
          y_bin * log(p_clipped) + (1 - y_bin) * log(1 - p_clipped)
        )

        performance <- c(
          cm$overall[c("Accuracy", "Kappa")],
          cm$byClass[c("Recall", "Specificity",
                       "Precision", "F1")],
          "Balanced_Accuracy" = cm$byClass[["Balanced Accuracy"]],
          "G_mean" = unname(gmean),
          "ROC_AUC" = roc_obj$auc,
          "PR_AUC" = unname(pr_obj$auc.integral),
          "Cross_Entropy" = cross_entropy,
          "Cross_Entropy_Efficient" = cross_entropy_efficient
        )

        # save a copy
        same_performance <- performance

        best_imbalance_rate <- as.data.frame(best_balance_i)
        names(best_imbalance_rate) <- "imbalance_rate"

        best_hyperparameter <- best_results_train[1, names_tuneGrid]

        if (method_i == "glm") {

          parameters <- methods[["glm"]]

          if (parameters[["weights"]][1] == "dinamic") {
            w0 <- nrow(datasets_to_train[[as.character(best_balance_i)]]) / (2 * length(which(datasets_to_train[[as.character(best_balance_i)]]$class_efficiency == "not_efficient")))
            w1 <- nrow(datasets_to_train[[as.character(best_balance_i)]]) / (2 * length(which(datasets_to_train[[as.character(best_balance_i)]]$class_efficiency == "efficient")))
          } else if (is.data.frame(parameters[["weights"]])) {
            w0 <- parameters[["weights"]][["w0"]]
            w1 <- parameters[["weights"]][["w1"]]
          } else {
            w0 <- 1
            w1 <- 1
          }

          best_hyperparameter <- data.frame(w0, w1)

        }

        performance <- t(performance)

        performance <- cbind(
          best_imbalance_rate,
          best_hyperparameter,
          performance)

      }

      save_performance_validation[[method_i]] <- performance

    } # end hold out

  } # end method

  # show results without hold out
  if(is.null(hold_out)) {

    if (is.null(hold_out)) {

      output_PEAXAI_models <- list(
        best_model_fit = save_best_model_fit,
        performance_train = save_performance_train,
        performance_train_all = performance_train_all_by_dataset
      )

      return(output_PEAXAI_models)

    }

  }

  # It is necessary to determine the complete facets with ALL data
  get_imbalance <- function(df) as.numeric(sub("\\*$", "", as.character(df["imbalance_rate"])))

  best_imbalance <- sort(unique(unlist(lapply(save_performance_validation, get_imbalance))), na.last = NA)

  if (length(best_imbalance) == 1 & best_imbalance[1] == real_balance) {

    datasets_all <- list(all_data)
    real_balance <- round(prop.table(table(all_data$class_efficiency)), 4)
    names(datasets_all)[1] <- paste0(as.character(round(real_balance, 4)),"*")

  } else {

    if(round(real_balance, 4) %in% best_imbalance) {

      best_imbalance_SMOTE <- best_imbalance[which(round(real_balance, 4) != best_imbalance)]

      datasets_all <- list(all_data)
      real_balance <- round(prop.table(table(all_data$class_efficiency)), 4)[1]
      names(datasets_all)[1] <- paste0(as.character(round(real_balance, 4)),"*")

      if (as.character(RTS) %in% RTS_available) {

        # determine complete facets
        train_data_SMOTE <- SMOTE_data(
          data = all_data,
          x = x,
          y = y,
          RTS = RTS,
          balance_data = best_imbalance_SMOTE,
          seed = seed
        )

        datasets_all <- append(train_data_SMOTE, datasets_all, after = 0)

      }

    } else {

      datasets_all <- list()
      best_imbalance_SMOTE <- best_imbalance

      if (as.character(RTS) %in% RTS_available) {

        # determine complete facets
        train_data_SMOTE <- SMOTE_data(
          data = all_data,
          x = x,
          y = y,
          RTS = RTS,
          balance_data = best_imbalance_SMOTE,
          seed = seed
        )

        datasets_all <- append(train_data_SMOTE, datasets_all)

      }

    }

  }

  # train with ALL data: best balance + hypeparameters
  for (method_i in names(methods)) {

    imbalance_rate_i <- save_performance_validation[[method_i]]$imbalance_rate

    if (grepl("\\*$", imbalance_rate_i)) {
      imbalance_rate_i <- paste0(round(real_balance, 4), "*")
    }

    names_tuneGrid <- names(best_methods[[method_i]][["tuneGrid"]])

    if (method_i != "glm") {
      best_methods[[method_i]][["tuneGrid"]] <- save_performance_validation[[method_i]][names_tuneGrid]
    }

    model_fit <- train_PEAXAI(
      data = datasets_all[[as.character(imbalance_rate_i)]],
      method = method_i,
      parameters = best_methods[[method_i]],
      trControl = trainControl(method = "none", classProbs = TRUE),
      metric_priority = metric_priority,
      seed = seed
    )

    # actualize weights
    if (method_i == "glm") {

      parameters <- methods[["glm"]]

      if (parameters[["weights"]][1] == "dinamic") {
        w0 <- nrow(datasets_all[[as.character(imbalance_rate_i)]]) / (2 * length(which(datasets_all[[as.character(imbalance_rate_i)]]$class_efficiency == "not_efficient")))
        w1 <- nrow(datasets_all[[as.character(imbalance_rate_i)]]) / (2 * length(which(datasets_all[[as.character(imbalance_rate_i)]]$class_efficiency == "efficient")))
      } else if (is.data.frame(parameters[["weights"]])) {
        w0 <- parameters[["weights"]][["w0"]]
        w1 <- parameters[["weights"]][["w1"]]
      } else {
        w0 <- 1
        w1 <- 1
      }

      save_performance_validation[[method_i]]$w0 <- w0
      save_performance_validation[[method_i]]$w1 <- w1

    }

    save_best_model_fit[[method_i]] <- model_fit
  }

  output_PEAXAI_models <- list(
    best_model_fit = save_best_model_fit,
    performance_validation = save_performance_validation,
    performance_train = save_performance_train,
    performance_train_all = performance_train_all_by_dataset
    )

  return(output_PEAXAI_models)

}
