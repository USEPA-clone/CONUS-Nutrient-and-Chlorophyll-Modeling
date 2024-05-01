#Predictive modeling of Chl a in lakes/reservoirs within CONUS at HUC12 scale
#Chlorophyll data from: 
#Platt, L. R., Koenig, L., Padilla, J., Spaulding, S. A., Covert, A., & Murphy, J. C. (2023). Source code: A national harmonized dataset of discrete chlorophyll from lakes and streams (2005-2022) (v1.0). Zenodo. https://doi.org/10.5281/zenodo.7879199

#Read in data----------
#Read in joined TP and TN HUC12 predictor variable data from Excel
#Make sure the TN and TP raw data is sorted to the top rows in Excel
#Otherwise R will read in as class "logical" and not as numeric data
library(readxl)
HUC12_P_N_wetland <- 
  read_excel("Research/U.S. Datasets/WSIO Data/HUC12_N_P_filtered_WSIO_variables_final.xlsx", 
             sheet = "WSIO_HUC12_Variables")

#Read in chlorophyll data
#Make sure first column is numeric HUC12s and not another type of ID
#Uncorrected Chl a 
HUC12_chlorophyll_a_uncorrected <- 
  read.csv("~/Research/U.S. Datasets/US_chlorophyll/USGS_chlorophyll/HUC12_chlorophyll_a_uncorrected.csv")
#Corrected Chl a
HUC12_chlorophyll_a_corrected <- 
  read.csv("~/Research/U.S. Datasets/US_chlorophyll/USGS_chlorophyll/HUC12_chlorophyll_a_corrected.csv")

#Convert the trophic state category columns to data type factor
HUC12_chlorophyll_a_uncorrected$Trophic.State <- as.factor(HUC12_chlorophyll_a_uncorrected$Trophic.State)
HUC12_chlorophyll_a_corrected$Trophic.State <- as.factor(HUC12_chlorophyll_a_corrected$Trophic.State)

#Change name of HUC12 column
library(dplyr)
HUC12_chlorophyll_a_uncorrected <-
  HUC12_chlorophyll_a_uncorrected %>% rename(
  HUC12_join = HUC_12)

HUC12_chlorophyll_a_corrected <-
  HUC12_chlorophyll_a_corrected %>% rename(
    HUC12_join = HUC_12)

#Merge with HUC12 predictor variables
HUC12_chlorophyll_a_uncorrected_join <-
  merge(HUC12_chlorophyll_a_uncorrected, HUC12_P_N_wetland, by="HUC12_join")

HUC12_chlorophyll_a_corrected_join <-
  merge(HUC12_chlorophyll_a_corrected, HUC12_P_N_wetland, by="HUC12_join")

#Remove data if more than 10% is missing
HUC12_chlorophyll_a_uncorrected_join <-
  HUC12_chlorophyll_a_uncorrected_join[lapply(HUC12_chlorophyll_a_uncorrected_join, 
                             function(x) sum(is.na(x)) / length(x) ) < 0.1 ]

HUC12_chlorophyll_a_corrected_join <-
  HUC12_chlorophyll_a_corrected_join[lapply(HUC12_chlorophyll_a_corrected_join, 
                                              function(x) sum(is.na(x)) / length(x) ) < 0.1 ]

#convert columns that contain negative values to all positive numbers with a constant added 
#in case you log-transform the data
HUC12_chlorophyll_a_uncorrected_join$WETLANDS_CHG_PCT_HUC12 <- 
  9+HUC12_chlorophyll_a_uncorrected_join$WETLANDS_CHG_PCT_HUC12

HUC12_chlorophyll_a_corrected_join$WETLANDS_CHG_PCT_HUC12 <- 
  9+HUC12_chlorophyll_a_corrected_join$WETLANDS_CHG_PCT_HUC12

HUC12_chlorophyll_a_uncorrected_join$P_AG_BALANCE_HUC12 <-
  21+HUC12_chlorophyll_a_uncorrected_join$P_AG_BALANCE_HUC12

HUC12_chlorophyll_a_corrected_join$P_AG_BALANCE_HUC12 <-
  21+HUC12_chlorophyll_a_corrected_join$P_AG_BALANCE_HUC12

#Add a small value to every numeric dataframe (except HUC12 ID) 
#so the log10 of 0 isn't -inf if you choose to log10 transform the data
HUC12_chlorophyll_a_uncorrected_join <- 
  mutate_if(HUC12_chlorophyll_a_uncorrected_join[,-c(1)], is.numeric, ~.+0.000001)

HUC12_chlorophyll_a_corrected_join <- 
  mutate_if(HUC12_chlorophyll_a_corrected_join[,-c(1)], is.numeric, ~.+0.000001)

#Optional log10 transform, not needed
#HUC12_P_N_wetland_3 <- 
#mutate_if(HUC12_P_N_wetland_3, is.numeric, log10)

#Remove any TN and TP raw data (shouldn't be any), chlor-a raw data, geographic shape data, 
#number of samples, and N:P molar ratio data (shouldn't be any)
HUC12_chlorophyll_a_uncorrected_join <- 
  HUC12_chlorophyll_a_uncorrected_join[,-c(2:6,70)]

HUC12_chlorophyll_a_corrected_join <- 
  HUC12_chlorophyll_a_corrected_join[,-c(2:6,70)]

##Applying models to just two trophic classes of chlorophyll-a, eutrophic or not-------------

#Make new columns for Eutrophic or not
HUC12_chlorophyll_a_corrected_join$Eutrophic <-
  ifelse(HUC12_chlorophyll_a_corrected_join$Trophic.State=="Oligo-mesotrophic",
         "No","Yes")

HUC12_chlorophyll_a_uncorrected_join$Eutrophic <-
  ifelse(HUC12_chlorophyll_a_uncorrected_join$Trophic.State=="Oligo-mesotrophic",
         "No","Yes")

#Change to type factor
HUC12_chlorophyll_a_corrected_join$Eutrophic <-
  as.factor(HUC12_chlorophyll_a_corrected_join$Eutrophic)

HUC12_chlorophyll_a_uncorrected_join$Eutrophic <-
  as.factor(HUC12_chlorophyll_a_uncorrected_join$Eutrophic)

##Random forest modeling---------------
library(rsample)      # data splitting 
library(randomForest) # basic implementation
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many machine learning models
library(tibble)
library(ggplot2)      # for plotting the data
library(tidyr)
library(broom)

#Running models to figure out optimal variables----------
#Corrected and uncorrected chlor-a data 

#Default basic random forest models, not fine-tuned
#On entire dataset
#With progress bar because it will be computationally expensive

#uncorrected chlor-a
#set seed for reproducibility
set.seed(123)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  default_TP_chlor_un <- ranger(
    dependent.variable.name = "Eutrophic", 
    data      = na.omit(HUC12_chlorophyll_a_uncorrected_join[,-c(1:3,65)]), 
    probability = FALSE, 
    importance  = 'impurity_corrected',
    num.trees = 500,
    mtry      = 4)
  setTxtProgressBar(pb, i)
}
close(pb)

default_TP_chlor_un

#Importance of predictor variables
#make a dataframe out of importance measures of variables
im <- data.frame(importance(default_TP_chlor_un)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.default_TP_chlor_un.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.default_TP_chlor_un.),y=importance.default_TP_chlor_un.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables (Chlor-a uncorrected)")+ylab("Importance (impurity corrected)")+
  xlab("Input Variables")

#corrected chlor-a
#set seed for reproducibility
set.seed(123)
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  default_TP_chlor_cor <- ranger(
    dependent.variable.name = "Eutrophic", 
    data      = na.omit(HUC12_chlorophyll_a_corrected_join[,-c(1:3,65)]), 
    probability = FALSE, 
    importance  = 'impurity_corrected',
    num.trees = 500,
    mtry      = 4)
  setTxtProgressBar(pb, i)
}
close(pb)

default_TP_chlor_cor

#Importance of predictor variables
#make a dataframe out of importance measures of variables
im <- data.frame(importance(default_TP_chlor_cor)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.default_TP_chlor_cor.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn,importance.default_TP_chlor_cor.),y=importance.default_TP_chlor_cor.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables (Chlor-a corrected)")+ylab("Importance (impurity corrected)")+
  xlab("Input Variables")

##Optimal chlor-a corrected and uncorrected predictive models---------
#for eutrophic or non-eutrophic classification

##Check for collinearity first with top predictor variables
library(collinear)

#Uncorrected Chlor-a
Uncor_chlor_best <- HUC12_chlorophyll_a_uncorrected_join[,c(19,48)]

selected_variables <- collinear(
  df= Uncor_chlor_best, #your data frame
  response=NULL, #name of your response variable
  predictors=NULL, #names of your predictors,
  preference_order=NULL, #your predictors in order of interest
  max_cor=0.7, #maximum bivariate correlation
  max_vif=5, #maximum variance inflation factor
  encoding_method="mean", #method to convert categorical predictors into numerics
)

selected_predictors_vif <- vif_df(
  df = Uncor_chlor_best,
  response = NULL,
  predictors = selected_variables
)

#Corrected Chlor-a
Cor_chlor_best <- HUC12_chlorophyll_a_corrected_join[,c(14,30)]

selected_variables <- collinear(
  df= Cor_chlor_best, #your data frame
  response=NULL, #name of your response variable
  predictors=NULL, #names of your predictors,
  preference_order=NULL, #your predictors in order of interest
  max_cor=0.7, #maximum bivariate correlation
  max_vif=5, #maximum variance inflation factor
  encoding_method="mean", #method to convert categorical predictors into numerics
)

selected_predictors_vif <- vif_df(
  df = Cor_chlor_best,
  response = NULL,
  predictors = selected_variables
)

#Using only top predictor variables from default random forest model
# hyperparameter grid search
# reset for each optimal model for loop 
hyper_grid <- expand.grid(
  mtry       = seq(1, 2, by = 1),
  node_size  = seq(2, 9, by = 2),
  sample_size = c(.55, .632, .70, .80),
  OOB_RMSE   = 0)

# total number of combinations
nrow(hyper_grid)
## [1] 32

#Subset the data for best variables (as done above), but make sure to include Eutrophic binary classification
#Uncorrected Chlor-a
Uncor_chlor_best_vars <- HUC12_chlorophyll_a_uncorrected_join[,c(19,48,67)]
#Corrected Chlor-a
Cor_chlor_best_vars <- HUC12_chlorophyll_a_corrected_join[,c(14,30,67)]

#for-loop through grid of parameters to find best combination for lowest OOB RMSE
#perform for both subsets of data (corrected and uncorrected)--need to refresh the hyper_grid each time (rerun above code to reset the hyper_grid dataframe)
# for reproduciblity
set.seed(123)
#Uncorrected Chlor-a
for(i in 1:nrow(hyper_grid)) {
  uncorrected_optimal_model <- ranger(
    dependent.variable.name = "Eutrophic", 
    data            = na.omit(Uncor_chlor_best_vars), 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 123, 
    probability     = FALSE,
    importance      = 'impurity_corrected')
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(uncorrected_optimal_model$prediction.error)
}    

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

#Corrected Chlor-a
# for reproduciblity
set.seed(123)
for(i in 1:nrow(hyper_grid)) {
  corrected_optimal_model <- ranger(
    dependent.variable.name = "Eutrophic", 
    data            = na.omit(Cor_chlor_best_vars), 
    num.trees       = 500,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    seed            = 123, 
    probability     = FALSE,
    importance      = 'impurity_corrected')
  # add OOB error to grid
  hyper_grid$OOB_RMSE[i] <- sqrt(corrected_optimal_model$prediction.error)
}    

hyper_grid %>% 
  dplyr::arrange(OOB_RMSE) %>%
  head(10)

#Run optimal model
#For corrected chlor-a
# for reproduciblity
set.seed(123)
#With progress bar because it will be computationally expensive
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  optimal_ranger_cor <- ranger(
    dependent.variable.name = "Eutrophic", 
    data            = na.omit(Cor_chlor_best_vars), 
    num.trees       = 500,
    mtry            = 1,
    min.node.size   = 8, 
    sample.fraction = 0.550, probability=FALSE, #Run with probability=TRUE for partial dependency plots
    importance      = 'impurity')
  setTxtProgressBar(pb, i)
}
close(pb)

optimal_ranger_cor

#For uncorrected chlor-a
# for reproduciblity
set.seed(123)
#With progress bar because it will be computationally expensive
pb <- txtProgressBar(min = 0, max = 100, style = 3)
for(i in 1:100) {
  optimal_ranger_un <- ranger(
    dependent.variable.name = "Eutrophic", 
    data            = na.omit(Uncor_chlor_best_vars), 
    num.trees       = 500,
    mtry            = 1,
    min.node.size   = 2, 
    sample.fraction = 0.550, probability=FALSE, #Run with probability=TRUE for partial dependency plots
    importance      = 'impurity')
  setTxtProgressBar(pb, i)
}
close(pb)

optimal_ranger_un

##Importance of predictor variables

#Corrected Chlor-a
#make a dataframe out of importance measures of variables
im <- data.frame(importance(optimal_ranger_cor)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.optimal_ranger_cor.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.optimal_ranger_cor.),y=importance.optimal_ranger_cor.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables (Corrected Chlor-a)")+ylab("Importance (impurity)")+
  xlab("Input Variables")

#Uncorrected Chlor-a
#make a dataframe out of importance measures of variables
im <- data.frame(importance(optimal_ranger_un)) 

#convert rownames to a column
# Apply rownames_to_column on the copy of 
# DataFrame and put name of function rn
im <- tibble::rownames_to_column(im, "rn") 

#Plot the top indicator variables
library(dplyr)
im1 <- im %>% dplyr::arrange(desc(importance.optimal_ranger_un.)) %>% dplyr::top_n(50) 

ggplot(data=im1, aes(x=reorder(rn, importance.optimal_ranger_un.),y=importance.optimal_ranger_un.)) +
  geom_bar(stat="identity") +
  coord_flip()+
  ggtitle("Order of important variables (Uncorrected Chlor-a)")+ylab("Importance (impurity)")+
  xlab("Input Variables")


#Partial dependency plots
#Must have probability = TRUE in ranger model call
library(pdp)
library(gridExtra)

#Corrected chlor-a
grid.arrange(
  partial(optimal_ranger_cor, "WOODY_VEG_PCT_HUC12", plot = TRUE),
  partial(optimal_ranger_cor, "Average.Annual.Temperature", plot=TRUE),
  ncol = 2
)

#Uncorrected chlor-a
grid.arrange(
  partial(optimal_ranger_un, "RUNOFF_SPRING_HIST_HUC12", plot = TRUE),
  partial(optimal_ranger_un, "N_SYNTH_FERT_APP_HUC12",  plot = TRUE),
  ncol = 2
)

##Predict chlor-a in other HUC12s-------------

#Merge chlorophyll-a corrected and uncorrected data with HUC12 variable data, but keeping all HUC12 variable data
#Left join
HUC12_chlorophyll_a_uncorrected_join_full <-
  left_join(HUC12_P_N_wetland, HUC12_chlorophyll_a_uncorrected, by="HUC12_join")

HUC12_chlorophyll_a_corrected_join_full <-
  left_join(HUC12_P_N_wetland, HUC12_chlorophyll_a_corrected, by="HUC12_join")

#Export datasets for external use
library(openxlsx)
list_of_datasets <- list("Uncorrected chlorophyll-a" = HUC12_chlorophyll_a_uncorrected_join_full,
                         "Corrected chlorophyll-a" = HUC12_chlorophyll_a_corrected_join_full)
write.xlsx(list_of_datasets, file = "chlor_a_predictor_variables_combined.xlsx")

##Extrapolate model to estimate high/low chlor-a across all of CONUS at HUC12 scale
#Filter for just HUC12s of unknown chlor-a corrected or uncorrected
HUC12_predict_chlor_un <- HUC12_chlorophyll_a_uncorrected_join_full %>%
  dplyr::filter(is.na(Trophic.State))
HUC12_predict_chlor_cor <- HUC12_chlorophyll_a_corrected_join_full %>%
  dplyr::filter(is.na(Trophic.State))

#Make predictions based on optimal ranger model on data without chlor-a

#If used log10 transformed data...
#Make sure to log10 transform the data in order to fit the model formula appropriately!
#HUC12_predict <-
# mutate_if(HUC12_predict, is.numeric, log10)

#For uncorrected chlor-a
#Subset data to account for missing data
HUC12_predict_uncorrected <- 
  na.omit(HUC12_predict_chlor_un[,c(16,45,63)])
pred_randomForest1 <-  
  predict(optimal_ranger_un, HUC12_predict_uncorrected)

#For corrected chlor-a
#Subset data to account for missing data
HUC12_predict_corrected <- 
  na.omit(HUC12_predict_chlor_cor[,c(11,27,63)])
pred_randomForest2 <-  
  predict(optimal_ranger_cor, HUC12_predict_corrected)

# convert to a data frame
#For uncorrected chlor-a
pred_randomForest1 <- data.frame(pred_randomForest1)

#For corrected chlor-a
pred_randomForest2 <- data.frame(pred_randomForest2)

#Merge with HUC12s
#For uncorrected chlor-a
HUC12s_predict_un_chlora <- 
  cbind(HUC12_predict_uncorrected[,c(3)], pred_randomForest1)

#For corrected chlor-a
HUC12s_predict_cor_chlora <- 
  cbind(HUC12_predict_corrected[,c(3)], pred_randomForest2)

#Export to csv file for plotting in ArcGIS Pro with actual data
#uncorrected chlor-a
write.csv(HUC12s_predict_un_chlora, "HUC12_predict_uncorrected_chlora_final.csv")
#corrected chlor-a
write.csv(HUC12s_predict_cor_chlora, "HUC12_predict_corrected_chlora_final.csv")
