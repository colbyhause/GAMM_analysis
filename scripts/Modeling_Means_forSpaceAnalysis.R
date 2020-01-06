
library(tidyverse)
library(mgcv)
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(ClusterR)
library(clValid)
library(NbClust)
library(Distance)

#FOR SPATIAL DATA:
library(naniar)
library(tidyverse)
library(udunits2)
library(fpc)
library(vegan)
library(corrplot)
library(mclust)
library(tmap)
library(viridis)
library(sf)

#1.VERSION1----
#  This means model is using the RAW data, not the predicted data from the GAMM
# The ISSUE with this method is that some of the bins means are calculated from just one number, as some of the bin only have data from one of the transects. I would not trust to use this raw version as the final analysis. 

# read in csv with all flame data----
flame_all <- read_csv("data/Data_alltransects_withRKM_FIXED_120919.csv") 
flame_all <- flame_all[, -c(1:2)] # remove first 2 cols, done need them
# Remove the one weird temperature outlier:
flame_all <- flame_all %>% 
  filter(temp >=2)

# Divide each reach into bins .05 km bins (50 m)
# Have to have your column for distance be named "distance" for code to work
flame_all<- data.frame(flame_all)
colnames(flame_all)[1] <- "distance"

# create dataframe with the bins
binned_flame <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.05)) # the +1 is to prevent the code from cutting of the distance prematurely 

# also binning the data another way (I've seen Gabe use this) where you determine how the data is cut up based on the size bin you want:
max(binned_flame$distance)
146/.05
dat <- data.frame(binned_flame, cut(binned_flame$distance, breaks = 146/.05, labels = F))
colnames(dat)[18] <- "cuts"

# One day turn this into ONE for loop:

# CHL loop----
chl_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("CHL", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_chl_transect = mean(CHL)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_chl_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  chl_df  <- rbind(chl_df, dat_Tcounts)
}


write_csv(chl_df, "data_output/CHL_raw_means_for_clustering.csv")
chl_df <- read_csv("data_output/CHL_raw_means_for_clustering.csv")

# Turb loop----
turb_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("turb", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_turb_transect = mean(turb)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_turb_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  turb_df  <- rbind(turb_df, dat_Tcounts)
}

write_csv(turb_df, "data_output/turb_raw_means_for_clustering.csv")
turb_df<- read_csv("data_output/turb_raw_means_for_clustering.csv")

# NO3 loop-----
NO3_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("NO3", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_NO3_transect = mean(NO3)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_NO3_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  NO3_df  <- rbind(NO3_df, dat_Tcounts)
}

write_csv(NO3_df, "data_output/NO3_raw_means_for_clustering.csv")
NO3_df <- read_csv("data_output/NO3_raw_means_for_clustering.csv")

#pH loop----
pH_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("pH", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_pH_transect = mean(pH)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_pH_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  pH_df  <- rbind(pH_df, dat_Tcounts)
}

write_csv(pH_df, "data_output/pH_raw_means_for_clustering.csv")
pH_df<- read_csv("data_output/pH_raw_means_for_clustering.csv")
#DO loop----
DO_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("ODO", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_DO_transect = mean(ODO)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_DO_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  DO_df  <- rbind(DO_df, dat_Tcounts)
}

write_csv(DO_df, "data_output/DO_raw_means_for_clustering.csv")
DO_df<- read_csv("data_output/DO_raw_means_for_clustering.csv")

# spCond loop----
spCond_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("spCond", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_spCond_transect = mean(spCond)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_spCond_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  spCond_df  <- rbind(spCond_df, dat_Tcounts)
}

write_csv(spCond_df, "data_output/spCond_raw_means_for_clustering.csv")
spCond_df <- read_csv("data_output/spCond_raw_means_for_clustering.csv")


#fDOM loop----
fDOM_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("fDOM", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_fDOM_transect = mean(fDOM)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_fDOM_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  fDOM_df  <- rbind(fDOM_df, dat_Tcounts)
}

write_csv(fDOM_df, "data_output/fDOM_raw_means_for_clustering.csv")
fDOM_df <- read_csv( "data_output/fDOM_raw_means_for_clustering.csv")

# temp loop----
temp_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("temp", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_temp_transect = mean(temp)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_temp_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  temp_df  <- rbind(temp_df, dat_Tcounts)
}

write_csv(temp_df, "data_output/temp_raw_means_for_clustering.csv")
temp_df<- read_csv("data_output/temp_raw_means_for_clustering.csv")

# Make one dataframe of the grand mean of the parameters:
grandMean_allParams <- data.frame(cbind(temp.grand.mean = temp_df$grand_mean,spCond.grand.mean = spCond_df$grand_mean, pH.grand.mean = pH_df$grand_mean, ODO.grand.mean = DO_df$grand_mean, turb.grand.mean = turb_df$grand_mean, fDOM.grand.mean = fDOM_df$grand_mean, CHL.grand.mean = chl_df$grand_mean, NO3.grand.mean = NO3_df$grand_mean), chl_df[, c(2:10, 13)]) # the last part is just binding in the rest of the data that is the same for all params, so I just pulled it from the chl dataframe 

# write to a csv:----
write_csv(grandMean_allParams, "data_output/grandMeans.fromRaw.allParams.forClustering.csv")

# figure out where there is only 1 transect in the calculation:
NAmeans <- temp_df %>% 
  filter(tran_counts == 1)
NAmeans<- NAmeans[, c(1, 6:13)]
  

histo<- hist(NAmeans$cuts)
histo$counts
length(unique(NAmeans$cuts))

# VERSION 2 AND 3:----
# This means model uses the values predicted by the GAMM in Version2 in the script Modeling_VariationforTimeAnalysis.R
# The same code is used for version 2 and 3, the only difference is the dataframe that is used at the start,

# Read in df of predicted values estimated by GAMM:

#This is the data for version3, made in #2a. of script predict.values.from.GAMMS.R. This dataframe has one predicted measurement for each 0.05 increment rkm:

dat_predV3 <- read_csv("data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams.csv") 
# Dataset for version 2: # this is the dataframe is made in #2b of script predict.values.from.GAMMS.R. It has one predicted measurement for each actual rkm
#Note: the n in the dataframe is the number of coords that went into the calculate of the average coordinate ( see predict.values.from.GAMMS.R for more detailed explaination):

dat_predV2 <- read_csv("data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams_across_actualRKMs_withGeomeans.csv")


# choose which version you are doing:
# dat_pred <- dat_predV3 # assign either the dataset from V3 or V2 here
dat_pred <- dat_predV2 # assign either the dataset from V3 or V2 here 
rkm_list <- unique(dat_pred$rkm)

# CHL loop----
chl_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("CHL.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(CHL.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  chl_pred_df  <- rbind(chl_pred_df, dat_Tcounts)
}

# Turb loop----
turb_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("turb.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(turb.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  turb_pred_df  <- rbind(turb_pred_df, dat_Tcounts)
}

#NO3 loop----
NO3_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("NO3.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(NO3.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  NO3_pred_df  <- rbind(NO3_pred_df, dat_Tcounts)
}

# pH loop----
pH_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("pH.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(pH.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  pH_pred_df  <- rbind(pH_pred_df, dat_Tcounts)
}

# DO loop----
DO_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("DO.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(DO.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  DO_pred_df  <- rbind(DO_pred_df, dat_Tcounts)
}

#SpCond loop----
spCond_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("spCond.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(spCond.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  spCond_pred_df  <- rbind(spCond_pred_df, dat_Tcounts)
}

#fDOM loop----
fDOM_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("fDOM.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(fDOM.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  fDOM_pred_df  <- rbind(fDOM_pred_df, dat_Tcounts)
}

#Temp loop----
temp_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("temp.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(mean = mean(temp.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  temp_pred_df  <- rbind(temp_pred_df, dat_Tcounts)
}

# WRITE TO FILES:
# Combine to make one file:----
# If used Version 3 data:
#CV_allParams_predicted_0.05km <- data.frame(cbind(CHL.CV.predicted = chl_pred_df$cv,turb.CV.predicted =  turb_pred_df$cv, NO3.CV.predicted = NO3_pred_df$cv, spCond.CV.predicted = spCond_pred_df$cv, temp.CV.predicted = temp_pred_df$cv, pH.CV.predicted = pH_pred_df$cv, fDOM.CV.predicted = fDOM_pred_df$cv, DO.CV.predicted = DO_pred_df$cv, chl_pred_df[, 2:5])) 

#write to a csv:
#write_csv(CV_allParams_predicted_0.05km, "data_output/CV/CV_fromPredicted.05km_AllParams.csv") ### this is for version 3

# If used version2 data:
mean_allParams_predicted_V2<- data.frame(cbind(CHL.mean.predicted = chl_pred_df$mean,turb.mean.predicted =  turb_pred_df$mean, NO3.mean.predicted = NO3_pred_df$mean, spCond.mean.predicted = spCond_pred_df$mean, temp.mean.predicted = temp_pred_df$mean, pH.mean.predicted = pH_pred_df$mean, fDOM.mean.predicted = fDOM_pred_df$mean, DO.mean.predicted = DO_pred_df$mean, chl_pred_df[,5], dat_predV2[, c(1, 11:13)]))


#write to a csv:
#write_csv(CV_allParams_predicted_0.05km, "data_output/CV/CV_fromPredicted.05km_AllParams.csv") ### csv for version 3

write_csv(mean_allParams_predicted_V2, "data_output/Mean/Mean_fromPredicted_AllParams_across_actualRKMS.csv") ### csv for version 2
# checking trancounts:
which(mean_allParams_predicted_V2$tran_counts == 1) # should not have any 1s, only 2s and 3s

# WORKFLOW: now feed this data output into the Mean section of the ClusterCode_SpaceAndTimeModeling.R script 

