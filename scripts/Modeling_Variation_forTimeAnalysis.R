# The code below make the csvs for each version of the CV (coefficient of variation) analysis. 
#1.---- The first is calculating the CV using raw data: it takes the mean value of each parameter in each transect for each bin, and then takes those means and calculates the cv.The issue I ran into with this is the CV can't be calculated from just one number, and there are some bins that only had one transect in them. To prevent this, I had to increase the bin size to .35, which really loses our resolution.See ISSUES section below.
#2.---- The second version takes values from the GAMM that are predicted for each EXISTING rkm. The GAMM predicts a parameter values for every existing rkm value (this is what I set the new data "newd" to predcit across) so i get one predicted value per transect for each rkm. I then don't have to get a mean value first, I just get the cv on the 3 transect values for each rkm.
#3.---- The third version take the predicted values from the GAMM. The GAMM predicts a parameter values for every 50 km (this is what Iset the "new data" for it to predcit across to) so i get one predicted value per 50 km per transect. I then don't have to get a mean value first, I just get the cv on the 3 transect values per 50 km. Issue: then i also lose the coords associated witht these points, and have to estimate those ( see script "CalculatingCoords_forPredictedrkms)

library(sjstats)
library(Distance)
library(tidyverse)

# Version 1: CV on the raw data:----
# Set up data----

# read in csv with all flame data
flame_all <- read_csv("data/Data_alltransects_withRKM_FIXED_120919.csv") 
flame_all <- flame_all[, -c(1:2)] # remove first 2 cols, dont need them
# Remove the one weird temperature outlier:
flame_all <- flame_all %>% 
  filter(temp >=2)


# Divide each reach into bins .05 km bins (50 m)
# Have to have your column for distance be named "distance" for code to work
flame_all<- data.frame(flame_all)
colnames(flame_all)[1] <- "distance"


# create dataframe with the bins
binned_flame <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.35)) # the +1 is to prevent the code from cutting of the distance prematurely 

# also binning the data another way (I've seen Gabe use this) where you determine how the data is cut up based on the size bin you want:
max(binned_flame$distance)
146/.05
dat <- data.frame(binned_flame, cut(binned_flame$distance, breaks = 146/.05, labels = F))
colnames(dat)[18] <- "cuts"

cut_list<- unique(dat$cuts)

#### ISSUE:----
#### this bin is too small because it leaves some bins with only one transect in them, and the coeff of variation on one number is returned as an NA 
#### need to determine the bin size so there is no bin with only 1 transect in it 
#### Using the .05 bin size that leave 690 observations with only one transect in it. I determined this was the issue because for each parameter there would be 690 observations with NAs in them after running the CV code. Also online it says the cv function returns NA for only 1 number.
# This is how i determined number of obs with NAs:
test <- which(is.na(turb_df))
test2<- turb_df[test, ]

# Now need to figure out what is the bin size to use to prevent just 1 transect per bin:
# using the turb code to test:
# 0.1 drops it down to 458
# 0.15 drops it down to 159
# 0.2 drops it down to 117
# .25 drops it to 36
# 0.3 drops it down to 52
# 0.35 drops it to 0 #### so would need to set bin to .35 km
# 0.4 drops it down to 4
146/.35
dat <- data.frame(binned_flame, cut(binned_flame$distance, breaks = 146/.35, labels = F))
colnames(dat)[18] <- "cuts"

cut_list<- unique(dat$cuts)

turb_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("turb", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_turb_transect = mean(turb)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_turb_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  turb_df  <- rbind(turb_df, dat_Tcounts)
}

which(turb_df$tran_counts == 1)
test <- which(is.na(turb_df))
test2<- turb_df[test, ]

#LOOPS:----

# CHL loop----
chl_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("CHL", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_chl_transect = mean(CHL)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_chl_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  chl_df  <- rbind(chl_df, dat_Tcounts)
}

write_csv(chl_df, "data_output/CV/CHL_raw_cv_for_clustering.csv")


# Turb loop----
turb_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("turb", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_turb_transect = mean(turb)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_turb_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  turb_df  <- rbind(turb_df, dat_Tcounts)
}

write_csv(turb_df, "data_output/CV/turb_raw_cv_for_clustering.csv")
# NO3 loop-----
NO3_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("NO3", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_NO3_transect = mean(NO3)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_NO3_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  NO3_df  <- rbind(NO3_df, dat_Tcounts)
}

write_csv(NO3_df, "data_output/CV/NO3_raw_cv_for_clustering.csv")
#pH loop----
pH_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("pH", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_pH_transect = mean(pH)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_pH_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  pH_df  <- rbind(pH_df, dat_Tcounts)
}

write_csv(pH_df, "data_output/CV/pH_raw_cv_for_clustering.csv")
pHtest<- dat[ , c("pH", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")]

which(is.na(pH_df) == T)
#DO loop----
DO_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("ODO", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_DO_transect = mean(ODO)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_DO_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  DO_df  <- rbind(DO_df, dat_Tcounts)
}

write_csv(DO_df, "data_output/CV/DO_raw_cv_for_clustering.csv")
# spCond loop----
spCond_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("spCond", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_spCond_transect = mean(spCond)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_spCond_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  spCond_df  <- rbind(spCond_df, dat_Tcounts)
}

write_csv(spCond_df, "data_output/CV/spCond_raw_cv_for_clustering.csv")
#fDOM loop----
fDOM_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("fDOM", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_fDOM_transect = mean(fDOM)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_fDOM_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  fDOM_df  <- rbind(fDOM_df, dat_Tcounts)
}
write_csv(fDOM_df, "data_output/CV/fDOM_raw_cv_for_clustering.csv")

# temp loop----
temp_df <- NULL
for (i in cut_list) {
  dat_df<- dat[ , c("temp", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_temp_transect = mean(temp)) %>% 
    group_by(cuts) %>% 
    mutate(grand_cv = cv(unique(mean_temp_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  temp_df  <- rbind(temp_df, dat_Tcounts)
}

write_csv(temp_df, "data_output/CV/temp_raw_cv_for_clustering.csv")

# Make one dataframe of the grand mean of the parameters:
CVraw_allParams <- data.frame(cbind(temp.grand.cv = temp_df$grand_cv,spCond.grand.cv = spCond_df$grand_cv, pH.grand.cv = pH_df$grand_cv, ODO.grand.cv = DO_df$grand_cv, turb.grand.cv = turb_df$grand_cv, fDOM.grand.cv = fDOM_df$grand_cv, CHL.grand.cv = chl_df$grand_cv, NO3.grand.cv = NO3_df$grand_cv), chl_df[, c(2:10, 13)]) # the last part is just binding in the rest of the data that is the same for all params, so I just pulled it from the chl dataframe 

# write to a csv:----
write_csv(CVraw_allParams, "data_output/CV/CV.fromRaw.allParams.forClustering.csv")

# Version 2 AND 3: CV on predicted GAMM values ----
#Read in df of predicted values estimated by GAMM:
dat_predV3 <- read_csv("data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams.csv") #this is the data for version3, made in #2a. of script predict.values.from.GAMMS.R. This dataframe has one predicted measurement for each 0.05 increment rkm

dat_predV2 <- read_csv("data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams_across_actualRKMs_withGeomeans.csv")
# this is the data for version 2, made in #2b of script predict.values.from.GAMMS.R
# This dataframe has one predicted measurement for each actual rkm
#Note: the n in the dataframe is the number of coords that went into the calculate of the average coordinate ( see predict.values.from.GAMMS.R for more detailed explaination)

# choose which version you are doing:
dat_pred <- dat_predV2
rkm_list <- unique(dat_pred$rkm)

# CHL loop----
chl_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- dat_pred[ , c("CHL.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(cv = cv(CHL.predicted))
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
    mutate(cv = cv(turb.predicted))
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
    mutate(cv = cv(NO3.predicted))
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
    mutate(cv = cv(pH.predicted))
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
    mutate(cv = cv(DO.predicted))
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
    mutate(cv = cv(spCond.predicted))
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
    mutate(cv = cv(fDOM.predicted))
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
    mutate(cv = cv(temp.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  temp_pred_df  <- rbind(temp_pred_df, dat_Tcounts)
}

# Combine to make one file:----
#CV_allParams_predicted_0.05km <- data.frame(cbind(CHL.CV.predicted = chl_pred_df$cv,turb.CV.predicted =  turb_pred_df$cv, NO3.CV.predicted = NO3_pred_df$cv, spCond.CV.predicted = spCond_pred_df$cv, temp.CV.predicted = temp_pred_df$cv, pH.CV.predicted = pH_pred_df$cv, fDOM.CV.predicted = fDOM_pred_df$cv, DO.CV.predicted = DO_pred_df$cv, chl_pred_df[, 2:5])) # thisnis for version 3

#write to a csv:
#write_csv(CV_allParams_predicted_0.05km, "data_output/CV/CV_fromPredicted.05km_AllParams.csv") ### this is for version 3

CV_allParams_predicted_0.05km <- data.frame(cbind(CHL.CV.predicted = chl_pred_df$cv,turb.CV.predicted =  turb_pred_df$cv, NO3.CV.predicted = NO3_pred_df$cv, spCond.CV.predicted = spCond_pred_df$cv, temp.CV.predicted = temp_pred_df$cv, pH.CV.predicted = pH_pred_df$cv, fDOM.CV.predicted = fDOM_pred_df$cv, DO.CV.predicted = DO_pred_df$cv, chl_pred_df[,5], dat_predV2[, c(1, 11:13)]))
# use this one for version 2

#write to a csv:
#write_csv(CV_allParams_predicted_0.05km, "data_output/CV/CV_fromPredicted.05km_AllParams.csv") ### csv for version 3

write_csv(CV_allParams_predicted_0.05km, "data_output/CV/CV_fromPredicted_AllParams_across_actualRKMS.csv") ### csv for version 2
####CHECK:
dat <- read_csv("data_output/CV/CV_fromPredicted_AllParams_across_actualRKMS.csv")
# checking trancounts:
which(dat$tran_counts == 2) # should not have any 1s, only 2s and 3s (2s because the first tansect did not go as far as 2 and 3, so some rkms only have T2 and T3 in them)

# WORKFLOW: now feed this data output into "TIME (VARIATION) MODEL" section of the cluster code in script: ClusterCode_SpaceandTimeModels.R

############################### just testing code below:
# test to make sure this is working:----
test_df <- subset(dat_pred, dat_pred$rkm == 0 |dat_pred$rkm == .05 | dat_pred$rkm == .1)
rkm_list <- unique(test_df$rkm)
turb_pred_df <- NULL
for (i in rkm_list) {
  dat_df<- test_df[ , c("turb.predicted", "transect", "rkm")] %>% 
    filter(rkm == i) %>% 
    group_by(rkm) %>% 
    mutate(cv = cv(turb.predicted))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  turb_pred_df  <- rbind(turb_pred_df, dat_Tcounts)
}

#Trying to turn this into a full loop:----------------------
param_list <- c("CHL.predicted", "turb.predicted","NO3.predicted", "pH.predicted", "DO.predicted", "spCond.predicted", "fDOM.predicted", "temp.predicted")

param_pred_df <- NULL
for (param in param_list) {
  dat_df<- dat_pred[ , c(param,"transect", "rkm")]
for (i in rkm_list) { 
  param_df<- dat_df %>% 
    filter(rkm == i) %>% 
    #group_by(transect) %>% 
    group_by(rkm) %>% 
    mutate(cv = cv(CHL.predicted))
  dat_Tcounts <- param_df %>% 
    mutate(tran_counts = length(unique(transect)))
  param_pred_df  <- rbind(param_df, dat_Tcounts)
  #return(param_pred_df)
}
  return(param_pred_df)
  write_csv(param_pred_df, paste0("data_output/CV/CV_", param, ".csv"))
}

 