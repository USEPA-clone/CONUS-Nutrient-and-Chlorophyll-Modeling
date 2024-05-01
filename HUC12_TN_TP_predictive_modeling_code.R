#Predictive modeling of TP and TN within CONUS at HUC12 scale based on input SPARROW TN and TP data
#Representative of streams throughout CONUS
#From: Saad, D. A., Argue, D.M., Schwarz, G.E., Anning, D.W., Ator, S.W., Hoos, A.B., Preston, S.D., Robertson, D.M., and Wise, D.R., 2019. Water-quality and streamflow datasets used for estimating long-term mean daily streamflow and annual loads to be considered for use in regional streamflow, nutrient and sediment SPARROW models, United States, 1999-2014.  (2019). https://doi.org:https://doi.org/10.5066/F7DN436B

#Read in joined TP and TN HUC12 data from Excel--------------
library(readxl)
#Make sure the TN and TP raw data is sorted to the top rows in Excel
#Otherwise R will read in as class "logical" and not as numeric data
HUC12_P_N_wetland <- 
  read_excel("Research/U.S. Datasets/WSIO Data/HUC12_N_P_filtered_WSIO_variables_final.xlsx", 
                                                        sheet = "WSIO_HUC12_Variables")

#Explore the proportion of missing data for this best fit variable dataset
library(naniar)  # for proportion of missing variables
library(UpSetR)
gg_miss_upset(HUC12_P_N_wetland,
              nsets = 15, nintersects = 20)


##Data cleaning and transform-----------------
library(dplyr)

#Filter out for only HUC12s with N:P ratio data
HUC12_P_N_wetland_1 <- HUC12_P_N_wetland %>% 
  dplyr::filter(mean_Mass_ratio_N_P != "")

#Transform Mean_TN and Mean_TP data into categorical data based on greater or less than the median (50th percentile) or 75th percentile of data
  HUC12_P_N_wetland_1$Cat_TN <- ifelse(HUC12_P_N_wetland_1$mean_Mean_TN >= quantile(HUC12_P_N_wetland_1$mean_Mean_TN,
                                                                                    prob=c(.50), type=1), 
                                       "High","Low")
  
  HUC12_P_N_wetland_1$Cat_TP <- ifelse(HUC12_P_N_wetland_1$mean_Mean_TP >= quantile(HUC12_P_N_wetland_1$mean_Mean_TP,
                                                                                    prob=c(.50),type=1), 
                                       "High","Low")

#Convert TN and TP category columns to type factor
HUC12_P_N_wetland_1$Cat_TN <- as.factor(HUC12_P_N_wetland_1$Cat_TN)
HUC12_P_N_wetland_1$Cat_TP <- as.factor(HUC12_P_N_wetland_1$Cat_TP)

#Remove data if more than 10% is missing
HUC12_P_N_wetland_2 <-
  HUC12_P_N_wetland_1[lapply(HUC12_P_N_wetland_1, 
                               function(x) sum(is.na(x)) / length(x) ) < 0.1 ]

#Prep the data in case a log-transform is applied...
#convert columns that contain negative values to all positive numbers with a constant added
HUC12_P_N_wetland_2$WETLANDS_CHG_PCT_HUC12 <- 
  9+HUC12_P_N_wetland_2$WETLANDS_CHG_PCT_HUC12

HUC12_P_N_wetland_2$P_AG_BALANCE_HUC12 <-
  21+HUC12_P_N_wetland_2$P_AG_BALANCE_HUC12

#Add a small value to every numeric dataframe (except HUC12 ID) 
#so the log10 of 0 isn't -inf if you log10 transform the variables at all
HUC12_P_N_wetland_3 <- 
  mutate_if(HUC12_P_N_wetland_2[,-c(63,66)], is.numeric, ~.+0.000001)

#Optional log10 transform, not needed
#HUC12_P_N_wetland_3 <- 
  #mutate_if(HUC12_P_N_wetland_3, is.numeric, log10)

#Remove TN and TP raw data, number of samples, and N:P molar ratio data
HUC12_P_N_wetland_3 <- HUC12_P_N_wetland_3[,-c(63,65:70)]

##Random forest modeling---------------
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(tibble)
library(ggplot2)      # for plotting the data
library(tidyr)
library(broom)

#Default basic random forest model to figure out most important predictor variables
#On entire dataset
#With progress bar because it will be computationally expensive
#Default TN model
#set seed for reproducibility
set.seed(123)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  default_TN <- ranger(
    dependent.variable.name = "Cat_TN", 
    data      = na.omit(HUC12_P_N_wetland_3[,-c(65)]), #remove either TP or TN category column depending on whether running default TP or TN model
    probability = FALSE, 
    importance  = 'impurity_corrected',
    num.trees = 500,
    mtry      = 24)
  setTxtProgressBar(pb, i)
}
close(pb)

default_TN

#Default TP model
#set seed for reproducibility
set.seed(123)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  default_TP <- ranger(
    dependent.variable.name = "Cat_TP", 
    data      = na.omit(HUC12_P_N_wetland_3[,-c(64)]), #remove either TP or TN category column depending on whether running default TP or TN model
    probability = FALSE, 
    importance  = 'impurity_corrected',
    num.trees = 500,
    mtry      = 24)
  setTxtProgressBar(pb, i)
}
close(pb)

default_TP

#p-values of variables
#For TN
TN_default_importance <- importance_pvalues(default_TN, method = "altmann", formula = Cat_TN ~ ., 
                   data = na.omit(HUC12_P_N_wetland_3[,-c(65)]))
#For TP
TP_default_importance <- importance_pvalues(default_TP, method = "altmann", formula = Cat_TP ~ ., 
                   data = na.omit(HUC12_P_N_wetland_3[,-c(64)]))

#Variable importance graphs-----
#Importance of predictor variables (run for default TN and TP ranger model)

#TN
#make a dataframe out of importance measures of variables
im <- data.frame(importance(default_TN)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.default_TN.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.default_TN.),y=importance.default_TN.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables (TN)")+ylab("Importance (corrected impurity)")+
  xlab("Input Variables")

#TP
#make a dataframe out of importance measures of variables
im <- data.frame(importance(default_TP)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.default_TP.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.default_TP.),y=importance.default_TP.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables (TP)")+ylab("Importance (corrected impurity)")+
  xlab("Input Variables")

#Now, figure out optimal random forest models---------

# subset the data with top predictors based on importance in default models
Best_var_train_P <- HUC12_P_N_wetland_3 %>%
dplyr::select(Average.Annual.Surface.Runoff,
       WOODY_VEG_PCT_HUC12, Cat_TP) #For phosphorus
Best_var_train_N <- HUC12_P_N_wetland_3 %>%
  dplyr::select(WOODY_VEG_PCT_HUC12, 
          N_RED_DRYDEP_HUC12,
                          Cat_TN) #For nitrogen

#Test for collinearity of the top predictor variables
#Only use variables that pass the conditions set for collinearity
#TN model
library(collinear)
selected_variables <- collinear(
  df= Best_var_train_N, #your data frame
  response=NULL, #name of your response variable
  predictors=NULL, #names of your predictors,
  preference_order=NULL, #your predictors in order of interest
  max_cor=0.7, #maximum bivariate correlation
  max_vif=5, #maximum variance inflation factor
  encoding_method="mean", #method to convert categorical predictors into numerics
)

selected_predictors_vif <- vif_df(
  df = Best_var_train_N,
  response = NULL,
  predictors = selected_variables
)

#TP Model
library(collinear)
selected_variables <- collinear(
  df= Best_var_train_P, #your data frame
  response=NULL, #name of your response variable
  predictors=NULL, #names of your predictors,
  preference_order=NULL, #your predictors in order of interest
  max_cor=0.7, #maximum bivariate correlation
  max_vif=5, #maximum variance inflation factor
  encoding_method="mean", #method to convert categorical predictors into numerics
)

selected_predictors_vif <- vif_df(
  df = Best_var_train_P,
  response = NULL,
  predictors = selected_variables
)


#Using only top predictor variables from default random forest models
# hyperparameter grid search
# reset for TP and TN
hyper_grid <- expand.grid(
  mtry       = seq(1, 2, by = 1),
  node_size  = seq(2, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0)

# total number of combinations
nrow(hyper_grid)
## [1] 32

# for reproduciblity
set.seed(123)
#for-loop through grid of parameters
#For either TP or TN
for(i in 1:nrow(hyper_grid)) {
  # train model
  optimal_model <- ranger(
    dependent.variable.name = "Cat_TN", 
    data            = Best_var_train_N, 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 123, 
    probability     = FALSE,
    importance      = 'impurity_corrected')
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(optimal_model$prediction.error)
}    

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

#Run optimal model
#For TN
#for reproduciblity
set.seed(123)
#With progress bar because it will be computationally expensive
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  optimal_ranger_N <- ranger(
    dependent.variable.name = "Cat_TN", 
    data            = Best_var_train_N, 
    num.trees       = 500,
    mtry            = 1,
    min.node.size   = 6, 
    sample.fraction = 0.632, probability=FALSE, #Run with probability=TRUE for partial dependency plots
    importance      = 'impurity')
  setTxtProgressBar(pb, i)
}
close(pb)

optimal_ranger_N

#For TP
# for reproduciblity
set.seed(123)
#With progress bar because it will be computationally expensive
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  optimal_ranger_P <- ranger(
    dependent.variable.name = "Cat_TP", 
    data            = Best_var_train_P, 
    num.trees       = 500,
    mtry            = 1,
    min.node.size   = 8, 
    sample.fraction = 0.55, probability=FALSE, #Run with probability=TRUE for partial dependency plots
    importance      = 'impurity')
  setTxtProgressBar(pb, i)
}
close(pb)

optimal_ranger_P

#Importance of predictor variables TP
#make a dataframe out of importance measures of variables
im <- data.frame(importance(optimal_ranger_P)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.optimal_ranger_P.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.optimal_ranger_P.),y=importance.optimal_ranger_P.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables")+ylab("Importance (impurity)")+
  xlab("Input Variables")

#Importance of predictor variables TN
#make a dataframe out of importance measures of variables
im <- data.frame(importance(optimal_ranger_N)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.optimal_ranger_N.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.optimal_ranger_N.),y=importance.optimal_ranger_N.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables")+ylab("Importance (impurity)")+
  xlab("Input Variables")

#Partial dependency plots
#Must have probability = TRUE in ranger model call for optimal models
library(pdp)
library(gridExtra)

#TP
grid.arrange(
  partial(optimal_ranger_P, "WOODY_VEG_PCT_HUC12", plot = TRUE),
  partial(optimal_ranger_P, "RUNOFF_SPRING_HIST_HUC12",  plot = TRUE),
  ncol = 2
)

#TN
grid.arrange(
  partial(optimal_ranger_N, "TILE_DRAINAGE_PCT_HUC12", plot = TRUE),
  partial(optimal_ranger_N, "N_DEP_TOTAL_HUC12",  plot = TRUE),
  ncol = 2
)

##Predict TN and TP in other HUC12s-------------
#Read in Excel version of original HUC12 data without N:P molar ratios
#Renamed the predictor variables from random forest model to match the optimal_ranger model
library(readxl)
HUC12_CONUS_for_N_P_ratio_modeling_prediction <- 
  read_excel("Research/U.S. Datasets/SPARROW/HUC12_CONUS_for_N_P_ratio_modeling_prediction.xlsx", 
                                                            sheet = "HUC12_CONUS_for_N_P_ratio_model")
#Convert HUC12 join column to type character
HUC12_CONUS_for_N_P_ratio_modeling_prediction$HUC12_join <-
  as.character(HUC12_CONUS_for_N_P_ratio_modeling_prediction$HUC12_join)

##Extrapolate model to estimate high/low TN and TP across all of CONUS at HUC12 scale
#Filter for just HUC12s of unknown TN or TP data
 HUC12_predict <- HUC12_CONUS_for_N_P_ratio_modeling_prediction %>%
   dplyr::filter(Prediction == "Y")

#Make predictions based on optimal ranger model on data without TN and TP

 #If used log10 transformed data...
 #Make sure to log10 transform the data in order to fit the model formula appropriately!
 #HUC12_predict <-
  # mutate_if(HUC12_predict, is.numeric, log10)
 
 #For TN 50th percentile
 #Subset data to account for missing data
 HUC12_predict_TN_1 <- HUC12_predict %>%
   dplyr::select(N_RED_DRYDEP_HUC12,
                 WOODY_VEG_PCT_HUC12, HUC12_join)
 
 HUC12_predict_TN_1 <- na.omit(HUC12_predict_TN_1)
 
 pred_randomForest1 <-  
   predict(optimal_ranger_N, HUC12_predict_TN_1)
 
 #For TN 75th percentile
 #Subset data to account for missing data
 HUC12_predict_TN_2 <- HUC12_predict %>%
   dplyr::select(N_DEP_TOTAL_HUC12,
                 TILE_DRAINAGE_PCT_HUC12, HUC12_join)
 
 HUC12_predict_TN_2 <- na.omit(HUC12_predict_TN_2)
 
 pred_randomForest2 <-  
   predict(optimal_ranger_N, HUC12_predict_TN_2)
 
 #For TP 50th percentile
 #Subset data to account for missing data
 HUC12_predict_TP_1 <- HUC12_predict %>%
   dplyr::select(Average.Annual.Surface.Runoff,
                 WOODY_VEG_PCT_HUC12, HUC12_join)
 
 HUC12_predict_TP_1 <- na.omit(HUC12_predict_TP_1)
 
 pred_randomForest3 <-  
   predict(optimal_ranger_P, HUC12_predict_TP_1)
 
 #For TP 75th percentile
 #Subset data to account for missing data
 HUC12_predict_TP_2 <- HUC12_predict %>%
   dplyr::select(RUNOFF_SPRING_HIST_HUC12,
                 WOODY_VEG_PCT_HUC12, HUC12_join)
 
 HUC12_predict_TP_2 <- na.omit(HUC12_predict_TP_2)
 
 pred_randomForest4 <-  
   predict(optimal_ranger_P, HUC12_predict_TP_2)

# convert to a data frame
#For TN 50th percentile
pred_randomForest1 <- data.frame(pred_randomForest1)
#For TN 75th percentile
pred_randomForest2 <- data.frame(pred_randomForest2)

#For TP 50th percentile
pred_randomForest3 <- data.frame(pred_randomForest3)
#For TP 75th percentile
pred_randomForest4 <- data.frame(pred_randomForest4)

#Merge with HUC12s
#For TN 50th percentile
HUC12s_predict_TN_50th <- 
  cbind(HUC12_predict_TN_1[,c(3)], pred_randomForest1)
#For TN 75th percentile
HUC12s_predict_TN_75th <- 
  cbind(HUC12_predict_TN_2[,c(3)], pred_randomForest2)

#For TP 50th percentile
HUC12s_predict_TP_50th <- 
  cbind(HUC12_predict_TP_1[,c(3)], pred_randomForest3)
#For TP 75th percentile
HUC12s_predict_TP_75th <- 
  cbind(HUC12_predict_TP_2[,c(3)], pred_randomForest4)

#Export to csv file for plotting in ArcGIS Pro with actual data
#TN
write.csv(HUC12s_predict_TN_50th, "HUC12_predict_TN_50th_final.csv")
write.csv(HUC12s_predict_TN_75th, "HUC12_predict_TN_75th_final.csv")
#TP
write.csv(HUC12s_predict_TP_50th, "HUC12_predict_TP_50th_final.csv")
write.csv(HUC12s_predict_TP_75th, "HUC12_predict_TP_75th_final.csv")

