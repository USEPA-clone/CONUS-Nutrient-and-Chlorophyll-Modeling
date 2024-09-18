#K-fold cross validation on variations of each original "optimal" random forest model---------
#Fit multiple models to the data: regular random forest, spatial random forest (on the predicted regular random forest model residuals),
#spatial logistic regression (generalized linear model), and a regular generalized linear model (glm)

##Required libraries--------
library(tidyverse)
library(spmodel)
library(ranger)
library(pROC)
library(parallel)

###Chl-a model---------
# set seed (some randomness in ranger/big data spmodel)
set.seed(0)

# store functions
cv <- function(unique_index, data, local) {
  
  print(unique_index)
  
  # find index
  which_index <- which(data$index == unique_index)
  # # training/test data
  data_train <- data[-which_index, ]
  data_test <- data[which_index, ]
  
  # fit random forest models
  # "regular" random forest ranger model
  rf <- ranger(Eutrophic_num ~ WOODY_VEG_PCT_HUC12 +
                 ATEMP_SUMMER_HIST_HUC12 + N_RED_DRYDEP_HUC12, 
               data = data_train, probability = TRUE)
  
  
  rf_preds <- rf$predictions
  success_column <- which(colnames(rf_preds) == "1")  
  data_train$ranger_resid <- as.numeric(as.character(data_train$Eutrophic_num)) - rf_preds[, success_column]
  
  # model spatial residuals
  # spatial random forest model
  spmod <- splm(ranger_resid ~ 1, data = data_train, spcov_type = "exponential",
                xcoord = Long, ycoord = Lat, local = local)
  
  # predict for ranger models
  preds_ranger <- predict(rf, data = data_test)$predictions[, success_column]
  preds_spatial_resid <- predict(spmod, newdata = data_test, local = local)
  preds_ranger_spatial <- preds_ranger + preds_spatial_resid
  preds_ranger_spatial <- pmin(1, preds_ranger_spatial)
  preds_ranger_spatial <- pmax(0, preds_ranger_spatial)
  
  # fit glm model
  # binomial family, so essentially logistic regression
  gmod <- glm(
    formula = Eutrophic_num ~ WOODY_VEG_PCT_HUC12 + ATEMP_SUMMER_HIST_HUC12 + N_RED_DRYDEP_HUC12, 
    family = "binomial",
    data = data_train
  ) 
  
  # fir spatial glm model
  spgmod <- spglm(
    formula = Eutrophic_num ~ WOODY_VEG_PCT_HUC12 + ATEMP_SUMMER_HIST_HUC12 + N_RED_DRYDEP_HUC12, 
    family = "binomial",
    data = data_train,
    spcov_type = "exponential",
    xcoord = Long,
    ycoord = Lat,
    local = local
  )
  
  # predictions for glm models
  preds_gmod <- predict(gmod, newdata = data_test, type = "response")
  preds_spgmod <- predict(spgmod, newdata = data_test, type = "response", local = local)
  
  # return output
  list(
    obs = as.numeric(as.character(data_test$Eutrophic_num)),
    pred_ranger = preds_ranger,
    pred = preds_ranger_spatial,
    pred_gmod = preds_gmod,
    pred_spgmod = preds_spgmod
  )
}

#Function for generating model performance statistics
get_statistics <- function(out) {
  
  #observed values from models and predicted values from models
  obs <- as.numeric(as.character(unlist(lapply(out, function(x) x$obs))))
  pred_ranger <- unlist(lapply(out, function(x) x$pred_ranger))
  pred_ranger_spatial <- unlist(lapply(out, function(x) x$pred))
  pred_gmod <- unlist(lapply(out, function(x) x$pred_gmod))
  pred_spgmod <- unlist(lapply(out, function(x) x$pred_spgmod))
  
  #prediction error for each model
  ranger_error <- obs - pred_ranger
  spatial_error <- obs - pred_ranger_spatial
  gmod_error <- obs - pred_gmod
  spgmod_error <- obs - pred_spgmod
  
  # bias
  bias <- c(
    mean(ranger_error),
    mean(spatial_error),
    mean(gmod_error),
    mean(spgmod_error)
  )
  
  # mse
  mse <- c(
    # mse
    mean(ranger_error^2),
    mean(spatial_error^2),
    mean(gmod_error^2),
    mean(spgmod_error^2)
  )
  
  # rmse
  rmse <- sqrt(mse)
  
  # mae
  mae <- c(
    mean(abs(ranger_error)),
    mean(abs(spatial_error)),
    mean(abs(gmod_error)),
    mean(abs(spgmod_error))
  )
  
  # medae
  medae <- c(
    median(abs(ranger_error)),
    median(abs(spatial_error)),
    median(abs(gmod_error)),
    median(abs(spgmod_error))
  )
  
  #auc
  auc <- c(
    unlist(auc(roc(response = obs, predictor = pred_ranger, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_ranger_spatial, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_gmod, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_spgmod, direction = "<", quiet = TRUE)))
  )
  
  #classification accuracy based on 0.5 probability threshold
  accuracy <- c(
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger_spatial >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger_spatial >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_gmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_gmod >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_spgmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_spgmod >= 0.5, 1, 0))))
  )
  
  type <- c("rf", "rf", "glm", "glm")
  spatial <- c("nonspatial", "spatial", "nonspatial", "spatial")
  
  data.frame(
    type, spatial, bias, mse, rmse, mae, medae, auc, accuracy
  )
  
}


#Load in data
Chla <- read.csv('https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/Chla.csv')

#Convert dependent variable to type factor
Chla$Eutrophic_num <- factor(Chla$Eutrophic_num, levels = c("0", "1"))

# 10-fold cross validation
# k=10 here
unique_index <- 1:10
Chla$index <- sample(rep(unique_index, length.out = NROW(Chla)))

# if not parallel
# out <- lapply(unique_index, cv, data = Chla, local = TRUE)
# get_statistics(out)

# if parallel
nc <- 5
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(ranger)
  library(spmodel)
})
out_safely <- parLapply(cl, unique_index, safely(cv), data = Chla, local = FALSE)
stopCluster(cl)

# any errors?
sum(sapply(out_safely, function(x) is.null(x$result)))

# find statistics
out <- lapply(out_safely, function(x) x$result)
get_statistics(out)

###TN50 model---------
# set seed (some randomness in ranger/big data spmodel)
set.seed(0)

# store functions
cv <- function(unique_index, data, local) {
  
  print(unique_index)
  
  # find index
  which_index <- which(data$index == unique_index)
  # # training/test data
  data_train <- data[-which_index, ]
  data_test <- data[which_index, ]
  
  # fit random forest models
  # "regular" random forest ranger model
  rf <- ranger(Cat_TN ~ WOODY_VEG_PCT_HUC12 + N_RED_DRYDEP_HUC12, 
               data = data_train, probability = TRUE)
  
  
  rf_preds <- rf$predictions
  success_column <- which(colnames(rf_preds) == "1")  
  data_train$ranger_resid <- as.numeric(as.character(data_train$Cat_TN)) - rf_preds[, success_column]
  
  # model spatial residuals
  # spatial random forest model
  spmod <- splm(ranger_resid ~ 1, data = data_train, spcov_type = "exponential",
                xcoord = Long, ycoord = Lat, local = local)
  
  # predict for ranger models
  preds_ranger <- predict(rf, data = data_test)$predictions[, success_column]
  preds_spatial_resid <- predict(spmod, newdata = data_test, local = local)
  preds_ranger_spatial <- preds_ranger + preds_spatial_resid
  preds_ranger_spatial <- pmin(1, preds_ranger_spatial)
  preds_ranger_spatial <- pmax(0, preds_ranger_spatial)
  
  # fit glm model
  # binomial family, so essentially logistic regression
  gmod <- glm(
    formula = Cat_TN ~ WOODY_VEG_PCT_HUC12 + N_RED_DRYDEP_HUC12, 
    family = "binomial",
    data = data_train
  ) 
  
  # fir spatial glm model
  spgmod <- spglm(
    formula = Cat_TN ~ WOODY_VEG_PCT_HUC12 + N_RED_DRYDEP_HUC12, 
    family = "binomial",
    data = data_train,
    spcov_type = "exponential",
    xcoord = Long,
    ycoord = Lat,
    local = local
  )
  
  # predictions for glm models
  preds_gmod <- predict(gmod, newdata = data_test, type = "response")
  preds_spgmod <- predict(spgmod, newdata = data_test, type = "response", local = local)
  
  # return output
  list(
    obs = as.numeric(as.character(data_test$Cat_TN)),
    pred_ranger = preds_ranger,
    pred = preds_ranger_spatial,
    pred_gmod = preds_gmod,
    pred_spgmod = preds_spgmod
  )
}

#Function for generating model performance statistics
get_statistics <- function(out) {
  
  #observed values from models and predicted values from models
  obs <- as.numeric(as.character(unlist(lapply(out, function(x) x$obs))))
  pred_ranger <- unlist(lapply(out, function(x) x$pred_ranger))
  pred_ranger_spatial <- unlist(lapply(out, function(x) x$pred))
  pred_gmod <- unlist(lapply(out, function(x) x$pred_gmod))
  pred_spgmod <- unlist(lapply(out, function(x) x$pred_spgmod))
  
  #prediction error for each model
  ranger_error <- obs - pred_ranger
  spatial_error <- obs - pred_ranger_spatial
  gmod_error <- obs - pred_gmod
  spgmod_error <- obs - pred_spgmod
  
  # bias
  bias <- c(
    mean(ranger_error),
    mean(spatial_error),
    mean(gmod_error),
    mean(spgmod_error)
  )
  
  # mse
  mse <- c(
    # mse
    mean(ranger_error^2),
    mean(spatial_error^2),
    mean(gmod_error^2),
    mean(spgmod_error^2)
  )
  
  # rmse
  rmse <- sqrt(mse)
  
  # mae
  mae <- c(
    mean(abs(ranger_error)),
    mean(abs(spatial_error)),
    mean(abs(gmod_error)),
    mean(abs(spgmod_error))
  )
  
  # medae
  medae <- c(
    median(abs(ranger_error)),
    median(abs(spatial_error)),
    median(abs(gmod_error)),
    median(abs(spgmod_error))
  )
  
  #auc
  auc <- c(
    unlist(auc(roc(response = obs, predictor = pred_ranger, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_ranger_spatial, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_gmod, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_spgmod, direction = "<", quiet = TRUE)))
  )
  
  #classification accuracy based on 0.5 probability threshold
  accuracy <- c(
    as.vector(sum(diag(table(obs, 
      ifelse(pred_ranger >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
      ifelse(pred_ranger_spatial >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger_spatial >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
    ifelse(pred_gmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_gmod >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
   ifelse(pred_spgmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_spgmod >= 0.5, 1, 0))))
  )
  
  type <- c("rf", "rf", "glm", "glm")
  spatial <- c("nonspatial", "spatial", "nonspatial", "spatial")
  
  data.frame(
    type, spatial, bias, mse, rmse, mae, medae, auc, accuracy
  )
  
}

#Load in data
TN50 <- read.csv('https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/TN50.csv')

#Convert dependent variable to type factor
TN50$Cat_TN <- if_else(TN50$Cat_TN == "Low", 0, 1)
TN50$Cat_TN <- factor(TN50$Cat_TN, levels = c("0", "1"))

# 10-fold cross validation
# k=10 here
unique_index <- 1:10
TN50$index <- sample(rep(unique_index, length.out = NROW(TN50)))

# if not parallel
# out <- lapply(unique_index, cv, data = TN50, local = TRUE)
# get_statistics(out)

# if parallel
nc <- 5
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(ranger)
  library(spmodel)
})
out_safely <- parLapply(cl, unique_index, safely(cv), data = TN50, local = FALSE)
stopCluster(cl)

# any errors?
sum(sapply(out_safely, function(x) is.null(x$result)))

# find statistics
out <- lapply(out_safely, function(x) x$result)
get_statistics(out)

###TN75 model----------
# set seed (some randomness in ranger/big data spmodel)
set.seed(0)

# store functions
cv <- function(unique_index, data, local) {
  
  print(unique_index)
  
  # find index
  which_index <- which(data$index == unique_index)
  # # training/test data
  data_train <- data[-which_index, ]
  data_test <- data[which_index, ]
  
  # fit random forest models
  # "regular" random forest ranger model
  rf <- ranger(Cat_TN ~ TILE_DRAINAGE_PCT_HUC12 + N_DEP_TOTAL_HUC12, 
               data = data_train, probability = TRUE)
  
  
  rf_preds <- rf$predictions
  success_column <- which(colnames(rf_preds) == "1")  
  data_train$ranger_resid <- as.numeric(as.character(data_train$Cat_TN)) - rf_preds[, success_column]
  
  # model spatial residuals
  # spatial random forest model
  spmod <- splm(ranger_resid ~ 1, data = data_train, spcov_type = "exponential",
                xcoord = Long, ycoord = Lat, local = local)
  
  # predict for ranger models
  preds_ranger <- predict(rf, data = data_test)$predictions[, success_column]
  preds_spatial_resid <- predict(spmod, newdata = data_test, local = local)
  preds_ranger_spatial <- preds_ranger + preds_spatial_resid
  preds_ranger_spatial <- pmin(1, preds_ranger_spatial)
  preds_ranger_spatial <- pmax(0, preds_ranger_spatial)
  
  # fit glm model
  # binomial family, so essentially logistic regression
  gmod <- glm(
    formula = Cat_TN ~ TILE_DRAINAGE_PCT_HUC12 + N_DEP_TOTAL_HUC12, 
    family = "binomial",
    data = data_train
  ) 
  
  # fir spatial glm model
  spgmod <- spglm(
    formula = Cat_TN ~ TILE_DRAINAGE_PCT_HUC12 + N_DEP_TOTAL_HUC12, 
    family = "binomial",
    data = data_train,
    spcov_type = "exponential",
    xcoord = Long,
    ycoord = Lat,
    local = local
  )
  
  # predictions for glm models
  preds_gmod <- predict(gmod, newdata = data_test, type = "response")
  preds_spgmod <- predict(spgmod, newdata = data_test, type = "response", local = local)
  
  # return output
  list(
    obs = as.numeric(as.character(data_test$Cat_TN)),
    pred_ranger = preds_ranger,
    pred = preds_ranger_spatial,
    pred_gmod = preds_gmod,
    pred_spgmod = preds_spgmod
  )
}

#Function for generating model performance statistics
get_statistics <- function(out) {
  
  #observed values from models and predicted values from models
  obs <- as.numeric(as.character(unlist(lapply(out, function(x) x$obs))))
  pred_ranger <- unlist(lapply(out, function(x) x$pred_ranger))
  pred_ranger_spatial <- unlist(lapply(out, function(x) x$pred))
  pred_gmod <- unlist(lapply(out, function(x) x$pred_gmod))
  pred_spgmod <- unlist(lapply(out, function(x) x$pred_spgmod))
  
  #prediction error for each model
  ranger_error <- obs - pred_ranger
  spatial_error <- obs - pred_ranger_spatial
  gmod_error <- obs - pred_gmod
  spgmod_error <- obs - pred_spgmod
  
  # bias
  bias <- c(
    mean(ranger_error),
    mean(spatial_error),
    mean(gmod_error),
    mean(spgmod_error)
  )
  
  # mse
  mse <- c(
    # mse
    mean(ranger_error^2),
    mean(spatial_error^2),
    mean(gmod_error^2),
    mean(spgmod_error^2)
  )
  
  # rmse
  rmse <- sqrt(mse)
  
  # mae
  mae <- c(
    mean(abs(ranger_error)),
    mean(abs(spatial_error)),
    mean(abs(gmod_error)),
    mean(abs(spgmod_error))
  )
  
  # medae
  medae <- c(
    median(abs(ranger_error)),
    median(abs(spatial_error)),
    median(abs(gmod_error)),
    median(abs(spgmod_error))
  )
  
  #auc
  auc <- c(
    unlist(auc(roc(response = obs, predictor = pred_ranger, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_ranger_spatial, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_gmod, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_spgmod, direction = "<", quiet = TRUE)))
  )
  
  #classification accuracy based on 0.5 probability threshold
  accuracy <- c(
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger_spatial >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger_spatial >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_gmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_gmod >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_spgmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_spgmod >= 0.5, 1, 0))))
  )
  
  type <- c("rf", "rf", "glm", "glm")
  spatial <- c("nonspatial", "spatial", "nonspatial", "spatial")
  
  data.frame(
    type, spatial, bias, mse, rmse, mae, medae, auc, accuracy
  )
  
}

#Load in data
TN75 <- read.csv('https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/TN75.csv')

#Convert dependent variable to type factor
TN75$Cat_TN <- if_else(TN50$Cat_TN == "Low", 0, 1)
TN75$Cat_TN <- factor(TN50$Cat_TN, levels = c("0", "1"))

# 10-fold cross validation
# k=10 here
unique_index <- 1:10
TN75$index <- sample(rep(unique_index, length.out = NROW(TN75)))

# if not parallel
# out <- lapply(unique_index, cv, data = TN75, local = TRUE)
# get_statistics(out)

# if parallel
nc <- 5
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(ranger)
  library(spmodel)
})
out_safely <- parLapply(cl, unique_index, safely(cv), data = TN75, local = FALSE)
stopCluster(cl)

# any errors?
sum(sapply(out_safely, function(x) is.null(x$result)))

# find statistics
out <- lapply(out_safely, function(x) x$result)
get_statistics(out)

###TP50 model----------
# set seed (some randomness in ranger/big data spmodel)
set.seed(0)

# store functions
cv <- function(unique_index, data, local) {
  
  print(unique_index)
  
  # find index
  which_index <- which(data$index == unique_index)
  # # training/test data
  data_train <- data[-which_index, ]
  data_test <- data[which_index, ]
  
  # fit random forest models
  # "regular" random forest ranger model
  rf <- ranger(Cat_TP ~ WOODY_VEG_PCT_HUC12 + Average.Annual.Surface.Runoff, 
               data = data_train, probability = TRUE)
  
  
  rf_preds <- rf$predictions
  success_column <- which(colnames(rf_preds) == "1")  
  data_train$ranger_resid <- as.numeric(as.character(data_train$Cat_TP)) - rf_preds[, success_column]
  
  # model spatial residuals
  # spatial random forest model
  spmod <- splm(ranger_resid ~ 1, data = data_train, spcov_type = "exponential",
                xcoord = Long, ycoord = Lat, local = local)
  
  # predict for ranger models
  preds_ranger <- predict(rf, data = data_test)$predictions[, success_column]
  preds_spatial_resid <- predict(spmod, newdata = data_test, local = local)
  preds_ranger_spatial <- preds_ranger + preds_spatial_resid
  preds_ranger_spatial <- pmin(1, preds_ranger_spatial)
  preds_ranger_spatial <- pmax(0, preds_ranger_spatial)
  
  # fit glm model
  # binomial family, so essentially logistic regression
  gmod <- glm(
    formula = Cat_TP ~ WOODY_VEG_PCT_HUC12 + Average.Annual.Surface.Runoff, 
    family = "binomial",
    data = data_train
  ) 
  
  # fir spatial glm model
  spgmod <- spglm(
    formula = Cat_TP ~ WOODY_VEG_PCT_HUC12 + Average.Annual.Surface.Runoff, 
    family = "binomial",
    data = data_train,
    spcov_type = "exponential",
    xcoord = Long,
    ycoord = Lat,
    local = local
  )
  
  # predictions for glm models
  preds_gmod <- predict(gmod, newdata = data_test, type = "response")
  preds_spgmod <- predict(spgmod, newdata = data_test, type = "response", local = local)
  
  # return output
  list(
    obs = as.numeric(as.character(data_test$Cat_TP)),
    pred_ranger = preds_ranger,
    pred = preds_ranger_spatial,
    pred_gmod = preds_gmod,
    pred_spgmod = preds_spgmod
  )
}

#Function for generating model performance statistics
get_statistics <- function(out) {
  
  #observed values from models and predicted values from models
  obs <- as.numeric(as.character(unlist(lapply(out, function(x) x$obs))))
  pred_ranger <- unlist(lapply(out, function(x) x$pred_ranger))
  pred_ranger_spatial <- unlist(lapply(out, function(x) x$pred))
  pred_gmod <- unlist(lapply(out, function(x) x$pred_gmod))
  pred_spgmod <- unlist(lapply(out, function(x) x$pred_spgmod))
  
  #prediction error for each model
  ranger_error <- obs - pred_ranger
  spatial_error <- obs - pred_ranger_spatial
  gmod_error <- obs - pred_gmod
  spgmod_error <- obs - pred_spgmod
  
  # bias
  bias <- c(
    mean(ranger_error),
    mean(spatial_error),
    mean(gmod_error),
    mean(spgmod_error)
  )
  
  # mse
  mse <- c(
    # mse
    mean(ranger_error^2),
    mean(spatial_error^2),
    mean(gmod_error^2),
    mean(spgmod_error^2)
  )
  
  # rmse
  rmse <- sqrt(mse)
  
  # mae
  mae <- c(
    mean(abs(ranger_error)),
    mean(abs(spatial_error)),
    mean(abs(gmod_error)),
    mean(abs(spgmod_error))
  )
  
  # medae
  medae <- c(
    median(abs(ranger_error)),
    median(abs(spatial_error)),
    median(abs(gmod_error)),
    median(abs(spgmod_error))
  )
  
  #auc
  auc <- c(
    unlist(auc(roc(response = obs, predictor = pred_ranger, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_ranger_spatial, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_gmod, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_spgmod, direction = "<", quiet = TRUE)))
  )
  
  #classification accuracy based on 0.5 probability threshold
  accuracy <- c(
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger_spatial >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger_spatial >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_gmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_gmod >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_spgmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_spgmod >= 0.5, 1, 0))))
  )
  
  type <- c("rf", "rf", "glm", "glm")
  spatial <- c("nonspatial", "spatial", "nonspatial", "spatial")
  
  data.frame(
    type, spatial, bias, mse, rmse, mae, medae, auc, accuracy
  )
  
}

#Load in data
TP50 <- read.csv('https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/TP50.csv')

#Convert dependent variable to type factor
TP50$Cat_TP <- if_else(TP50$Cat_TP == "Low", 0, 1)
TP50$Cat_TP <- factor(TP50$Cat_TP, levels = c("0", "1"))

# 10-fold cross validation
# k=10 here
unique_index <- 1:10
TP50$index <- sample(rep(unique_index, length.out = NROW(TP50)))

# if not parallel
# out <- lapply(unique_index, cv, data = TP50, local = TRUE)
# get_statistics(out)

# if parallel
nc <- 5
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(ranger)
  library(spmodel)
})
out_safely <- parLapply(cl, unique_index, safely(cv), data = TP50, local = FALSE)
stopCluster(cl)

# any errors?
sum(sapply(out_safely, function(x) is.null(x$result)))

# find statistics
out <- lapply(out_safely, function(x) x$result)
get_statistics(out)

###TP75 model----------
# set seed (some randomness in ranger/big data spmodel)
set.seed(0)

# store functions
cv <- function(unique_index, data, local) {
  
  print(unique_index)
  
  # find index
  which_index <- which(data$index == unique_index)
  # # training/test data
  data_train <- data[-which_index, ]
  data_test <- data[which_index, ]
  
  # fit random forest models
  # "regular" random forest ranger model
  rf <- ranger(Cat_TP ~ WOODY_VEG_PCT_HUC12 + RUNOFF_SPRING_HIST_HUC12, 
               data = data_train, probability = TRUE)
  
  
  rf_preds <- rf$predictions
  success_column <- which(colnames(rf_preds) == "1")  
  data_train$ranger_resid <- as.numeric(as.character(data_train$Cat_TP)) - rf_preds[, success_column]
  
  # model spatial residuals
  # spatial random forest model
  spmod <- splm(ranger_resid ~ 1, data = data_train, spcov_type = "exponential",
                xcoord = Long, ycoord = Lat, local = local)
  
  # predict for ranger models
  preds_ranger <- predict(rf, data = data_test)$predictions[, success_column]
  preds_spatial_resid <- predict(spmod, newdata = data_test, local = local)
  preds_ranger_spatial <- preds_ranger + preds_spatial_resid
  preds_ranger_spatial <- pmin(1, preds_ranger_spatial)
  preds_ranger_spatial <- pmax(0, preds_ranger_spatial)
  
  # fit glm model
  # binomial family, so essentially logistic regression
  gmod <- glm(
    formula = Cat_TP ~ WOODY_VEG_PCT_HUC12 + RUNOFF_SPRING_HIST_HUC12, 
    family = "binomial",
    data = data_train
  ) 
  
  # fir spatial glm model
  spgmod <- spglm(
    formula = Cat_TP ~ WOODY_VEG_PCT_HUC12 + RUNOFF_SPRING_HIST_HUC12, 
    family = "binomial",
    data = data_train,
    spcov_type = "exponential",
    xcoord = Long,
    ycoord = Lat,
    local = local
  )
  
  # predictions for glm models
  preds_gmod <- predict(gmod, newdata = data_test, type = "response")
  preds_spgmod <- predict(spgmod, newdata = data_test, type = "response", local = local)
  
  # return output
  list(
    obs = as.numeric(as.character(data_test$Cat_TP)),
    pred_ranger = preds_ranger,
    pred = preds_ranger_spatial,
    pred_gmod = preds_gmod,
    pred_spgmod = preds_spgmod
  )
}

#Function for generating model performance statistics
get_statistics <- function(out) {
  
  #observed values from models and predicted values from models
  obs <- as.numeric(as.character(unlist(lapply(out, function(x) x$obs))))
  pred_ranger <- unlist(lapply(out, function(x) x$pred_ranger))
  pred_ranger_spatial <- unlist(lapply(out, function(x) x$pred))
  pred_gmod <- unlist(lapply(out, function(x) x$pred_gmod))
  pred_spgmod <- unlist(lapply(out, function(x) x$pred_spgmod))
  
  #prediction error for each model
  ranger_error <- obs - pred_ranger
  spatial_error <- obs - pred_ranger_spatial
  gmod_error <- obs - pred_gmod
  spgmod_error <- obs - pred_spgmod
  
  # bias
  bias <- c(
    mean(ranger_error),
    mean(spatial_error),
    mean(gmod_error),
    mean(spgmod_error)
  )
  
  # mse
  mse <- c(
    # mse
    mean(ranger_error^2),
    mean(spatial_error^2),
    mean(gmod_error^2),
    mean(spgmod_error^2)
  )
  
  # rmse
  rmse <- sqrt(mse)
  
  # mae
  mae <- c(
    mean(abs(ranger_error)),
    mean(abs(spatial_error)),
    mean(abs(gmod_error)),
    mean(abs(spgmod_error))
  )
  
  # medae
  medae <- c(
    median(abs(ranger_error)),
    median(abs(spatial_error)),
    median(abs(gmod_error)),
    median(abs(spgmod_error))
  )
  
  #auc
  auc <- c(
    unlist(auc(roc(response = obs, predictor = pred_ranger, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_ranger_spatial, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_gmod, direction = "<", quiet = TRUE))),
    as.vector(auc(roc(response = obs, predictor = pred_spgmod, direction = "<", quiet = TRUE)))
  )
  
  #classification accuracy based on 0.5 probability threshold
  accuracy <- c(
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_ranger_spatial >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_ranger_spatial >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_gmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_gmod >= 0.5, 1, 0)))),
    as.vector(sum(diag(table(obs, 
                             ifelse(pred_spgmod >= 0.5, 1, 0)))) / sum(table(obs, ifelse(pred_spgmod >= 0.5, 1, 0))))
  )
  
  type <- c("rf", "rf", "glm", "glm")
  spatial <- c("nonspatial", "spatial", "nonspatial", "spatial")
  
  data.frame(
    type, spatial, bias, mse, rmse, mae, medae, auc, accuracy
  )
  
}

#Load in data
TP75 <- read.csv('https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/TP75.csv')

#Convert dependent variable to type factor
TP75$Cat_TP <- if_else(TP75$Cat_TP == "Low", 0, 1)
TP75$Cat_TP <- factor(TP75$Cat_TP, levels = c("0", "1"))

# 10-fold cross validation
# k=10 here
unique_index <- 1:10
TP75$index <- sample(rep(unique_index, length.out = NROW(TP75)))

# if not parallel
# out <- lapply(unique_index, cv, data = TP75, local = TRUE)
# get_statistics(out)

# if parallel
nc <- 5
cl <- makeCluster(nc)
clusterEvalQ(cl, {
  library(ranger)
  library(spmodel)
})
out_safely <- parLapply(cl, unique_index, safely(cv), data = TP75, local = FALSE)
stopCluster(cl)

# any errors?
sum(sapply(out_safely, function(x) is.null(x$result)))

# find statistics
out <- lapply(out_safely, function(x) x$result)
get_statistics(out)


#Predictions from optimal models for remaining HUC12s-------------
#turn off scientific notation for proper upload of HUC12 IDs
options(scipen=123)

##Chl-a---------
# fit spatial glm model on all data (not just training data)
spgmod_Chla <- spglm(
 formula = Eutrophic_num ~ WOODY_VEG_PCT_HUC12 + ATEMP_SUMMER_HIST_HUC12 + N_RED_DRYDEP_HUC12,
  family = "binomial",
  data = Chla, local = FALSE,
  spcov_type = "exponential",
  xcoord = Long,
  ycoord = Lat)

# make predictions on remaining HUC12s using newdata from original Chl-a model code
# which contains HUC12 centroid coordinates
HUC12_predict_corrected <- read.csv("https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/HUC12_predict_corrected.csv")

HUC12_predict_corrected$preds <- 
  predict(spgmod_Chla, newdata = HUC12_predict_corrected, type = "response",
          local = FALSE)

#Export the data
write.csv(HUC12_predict_corrected, "Chla_HUC12_predictions.csv")

##TN50-------
# fit spatial glm model on all data (not just training data)
spgmod_TN50 <- spglm(
  formula = Cat_TN ~ WOODY_VEG_PCT_HUC12 + N_RED_DRYDEP_HUC12,
  family = "binomial", local = FALSE,
  data = TN50,
  spcov_type = "exponential",
  xcoord = Long,
  ycoord = Lat)

# make predictions on remaining HUC12s using newdata from original TN/TP model code
# which contains HUC12 centroid coordinates
HUC12_predict_TN50 <- read.csv("https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/HUC12_predict_TN_1.csv")

HUC12_predict_TN50$preds <- 
  predict(spgmod_TN50, newdata = HUC12_predict_TN50, type = "response",
          local = FALSE)

#Export the data
write.csv(HUC12_predict_TN50, "TN50_HUC12_predictions.csv")


##TN75-------
# fit spatial glm model on all data (not just training data)
spgmod_TN75 <- spglm(
  formula = Cat_TN ~ TILE_DRAINAGE_PCT_HUC12 + N_DEP_TOTAL_HUC12,
  family = "binomial",
  data = TN75, local = FALSE,
  spcov_type = "exponential",
  xcoord = Long,
  ycoord = Lat)

# make predictions on remaining HUC12s using newdata from original TN/TP model code
# which contains HUC12 centroid coordinates
HUC12_predict_TN75 <- read.csv("https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/HUC12_predict_TN_2.csv")

HUC12_predict_TN75$preds <- 
  predict(spgmod_TN75, newdata = HUC12_predict_TN75, 
          local = FALSE, type = "response")

#Export the data
write.csv(HUC12_predict_TN75, "TN75_HUC12_predictions.csv")


##TP50--------
# fit spatial glm model on all data (not just training data)
spgmod_TP50 <- spglm(
  formula = Cat_TP ~ WOODY_VEG_PCT_HUC12 + Average.Annual.Surface.Runoff,
  family = "binomial",
  data = TP50,  local = FALSE,
  spcov_type = "exponential",
  xcoord = Long,
  ycoord = Lat)

# make predictions on remaining HUC12s using newdata from original TN/TP model code
# which contains HUC12 centroid coordinates
HUC12_predict_TP50 <- read.csv("https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/HUC12_predict_TP_1.csv")

HUC12_predict_TP50$preds <- 
  predict(spgmod_TP50, newdata = HUC12_predict_TP50, 
          local = FALSE, type = "response")

#Export the data
write.csv(HUC12_predict_TP50, "TP50_HUC12_predictions.csv")

##TP75-------
# fit spatial glm model on all data (not just training data)
spgmod_TP75 <- spglm(
  formula = Cat_TP ~ WOODY_VEG_PCT_HUC12 + RUNOFF_SPRING_HIST_HUC12,
  family = "binomial",
  data = TP75,  local = FALSE,
  spcov_type = "exponential",
  xcoord = Long,
  ycoord = Lat)

# make predictions on remaining HUC12s using newdata from original TN/TP model code
# which contains HUC12 centroid coordinates
HUC12_predict_TP75 <- read.csv("https://raw.githubusercontent.com/USEPA/CONUS-Nutrient-and-Chlorophyll-Modeling/main/HUC12_predict_TP_2.csv")

HUC12_predict_TP75$preds <- 
  predict(spgmod_TP75, newdata = HUC12_predict_TP75, 
          local = FALSE, type = "response")

#Export the data
write.csv(HUC12_predict_TP75, "TP75_HUC12_predictions.csv")
