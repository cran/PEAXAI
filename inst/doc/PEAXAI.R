## ----setup, include=FALSE-----------------------------------------------------
# Enable a "fast mode" during R CMD check to keep CRAN timings under control
is_check <- ("CheckExEnv" %in% search()) || any(c("_R_CHECK_TIMINGS_","_R_CHECK_LICENSE_") %in% names(Sys.getenv()))

# Global knitr options for compact output and consistent figures
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE, message = FALSE,
  fig.width = 7, fig.height = 4.8, dpi = 96
)

# Seed for reproducibility
set.seed(1)

## -----------------------------------------------------------------------------
library(PEAXAI)
data("firms", package = "PEAXAI")

## -----------------------------------------------------------------------------
data <- subset(firms, autonomous_community == "Comunidad Valenciana")[, -ncol(firms)]
rm(firms)

## -----------------------------------------------------------------------------
str(data)
summary(data)

## -----------------------------------------------------------------------------
x <- c(1:4)
y <- c(5)

## -----------------------------------------------------------------------------
RTS <- "vrs"

## -----------------------------------------------------------------------------
imbalance_rate <- seq(0.05, 0.4, 0.05) 

## -----------------------------------------------------------------------------
methods <- list(
  "nnet" = list(
    tuneGrid = expand.grid(
      size = c(1, 5, 10, 20),
      decay = 10^seq(-5, -1, by = 1)
    ),
    maxit = 100,
    preProcess = c("center", "scale"),
    # # --- arguments nnet ---
    entropy = TRUE,
    skip = TRUE,
    maxit = 1000,
    MaxNWts = 100000,
    trace = FALSE,
    weights = NULL
  )
)

## -----------------------------------------------------------------------------
trControl <- list(
  method = "cv",
  number = 5
)

## -----------------------------------------------------------------------------
metric_priority <- c("Balanced_Accuracy", "F1", "ROC_AUC")

## -----------------------------------------------------------------------------
hold_out <- NULL

## -----------------------------------------------------------------------------
seed <- 1

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
models <- PEAXAI_fitting(
  data = data,
  x = x,
  y = y,
  RTS = RTS,
  imbalance_rate = imbalance_rate,
  methods = methods,
  trControl = trControl,
  metric_priority = metric_priority,
  hold_out = hold_out,
  verbose = TRUE,
  seed = seed
)

## ----echo=FALSE---------------------------------------------------------------
models$best_model_fit

## ----echo=FALSE---------------------------------------------------------------
models$performance_train

## -----------------------------------------------------------------------------
importance_method <- list(
    name = "SHAP",
    nsim = 200
)

## -----------------------------------------------------------------------------
relative_importance <- PEAXAI_global_importance(
  final_model = models[["best_model_fit"]][["nnet"]],
  x = x,
  y = y,
  explain_data = data,
  reference_data = data,
  importance_method = importance_method
)

## ----echo=FALSE---------------------------------------------------------------
relative_importance

## -----------------------------------------------------------------------------
efficiency_thresholds <- seq(0.75, 0.95, 0.1) 

## -----------------------------------------------------------------------------
relative_importance_custom <- t(matrix(
  data = c(0.2, 0, 0.2, 0.2, 0.4),
))
relative_importance_custom <- as.data.frame(relative_importance_custom)
names(relative_importance_custom) <- names(data)[c(x,y)]

## ----echo=FALSE---------------------------------------------------------------
relative_importance_custom

## -----------------------------------------------------------------------------
directional_vector <- list(
  relative_importance = relative_importance,
  baseline  = "mean"        
)

## -----------------------------------------------------------------------------
targets <- PEAXAI_counterfactuals(
  data = data,
  x = x,
  y = y,
  final_model = models[["best_model_fit"]][["nnet"]],
  efficiency_thresholds = efficiency_thresholds,
  directional_vector = directional_vector,
  n_expand = 0.25,
  n_grid = 50,
  max_y = 1,
  min_x = 1
)

## -----------------------------------------------------------------------------
head(targets[["0.85"]][["counterfactual_dataset"]], 10)

## -----------------------------------------------------------------------------
head(targets[["0.85"]][["inefficiencies"]], 10)

## -----------------------------------------------------------------------------
ranking <- PEAXAI_ranking(
  data = data,
  x = x,
  y = y,
  final_model = models[["best_model_fit"]][["nnet"]],
  rank_basis = "predicted"
)

ranking2 <- PEAXAI_ranking(
    data = data,
    x = x,
    y = y,
    final_model = models[["best_model_fit"]][["nnet"]],
    efficiency_thresholds = efficiency_thresholds,
    targets = targets,
    rank_basis = "attainable"
  )

## -----------------------------------------------------------------------------
head(round(ranking, 4), 50)

## -----------------------------------------------------------------------------
head(round(ranking2[["0.85"]], 4), 50)

## -----------------------------------------------------------------------------
peers <- PEAXAI_peer(
  data = data,
  x = x,
  y = y,
  final_model = models[["best_model_fit"]][["nnet"]],
  targets = targets,
  efficiency_thresholds = efficiency_thresholds,
  weighted = FALSE,
  relative_importance = relative_importance
)

## -----------------------------------------------------------------------------
head(peers, 50)

