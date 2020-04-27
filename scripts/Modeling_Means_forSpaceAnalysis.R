
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
library(geosphere)
#1.VERSION1----
#  This means model is using the RAW data, not the predicted data from the GAMM
# The ISSUE with this method is that some of the bins means are calculated from just one number, as some of the bin only have data from one of the transects. I would not trust to use this raw version as the final analysis. 

# read in csv with all flame data----
flame_all <- read_csv("data/Data_alltransects_withRKM_FIXED_120919.csv") 
#flame_all <- flame_all[, -c(1:2)] # remove first 2 cols, done need them
# Remove the one weird temperature outlier:
flame_all <- flame_all %>% 
  filter(temp >=2)

# Divide each reach into bins .05 km bins (50 m)
# Have to have your column for distance be named "distance" for code to work
flame_all<- data.frame(flame_all)
colnames(flame_all)[1] <- "distance"

# create dataframe with the bins
#binned_flame <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.05)) # the +1 is to prevent the code from cutting of the distance prematurely 

# Binning to see how big the bin have to be so there is no bin with <2 transects
binned_flame<- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.1)) 

#binned_flame.15 <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.15))
#binned_flame.2 <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.2)) 
#binned_flame.3 <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.3)) 
#binned_flame.25 <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.25)) 


# also binning the data another way (I've seen Gabe use this) where you determine how the data is cut up based on the size bin you want:
max(binned_flame$distance)
146/.05

dat <- data.frame(binned_flame, cut(binned_flame$distance, breaks = 146/.1, labels = F))

#dat <- data.frame(binned_flame.1, cut(binned_flame.25$distance, breaks = 146/.1, labels = F))

colnames(dat)[18] <- "cuts"

# Deciding to break river into 100 m bins, and just drop all the points that only have 1 transect in that bin
# Remove records with bins that only have 1 transect in them. use turb loop to do this:
test_df <- NULL
for (i in 1:max(dat$cuts)) {
  dat_df<- dat[ , c("turb", "distance", "lat", "lon", "reach", "transect", "transect_ID", "distbegin", "distend", "cuts")] %>% 
    filter(cuts == i) %>% 
    group_by(transect) %>% 
    mutate(mean_turb_transect = mean(turb)) %>% 
    group_by(cuts) %>% 
    mutate(grand_mean = mean(unique(mean_turb_transect)))
  dat_Tcounts <- dat_df %>% 
    mutate(tran_counts = length(unique(transect)))
  test_df  <- rbind(test_df, dat_Tcounts)
}

counts_w_1_at.1 <- which(test_df$tran_counts ==1)
length(counts_w_1_at.1)

flame_all_truncated <- flame_all[-counts_w_1_at.1, ]
nrow(flame_all_truncated) #use this df now
nrow(flame_all)

write_csv(flame_all_truncated, "data_output/Data_alltransects_withRKM_FIXED_120919_NoBinsWith1_obs.csv")

#counts_w_1_at.25 <- which(chl_df$tran_counts <=1) # this is the smallest distance with the least number of data lost
#counts_w_1_at.2 <- which(chl_df$tran_counts <=1)
#counts_w_1_at.15 <- which(chl_df$tran_counts <=1)
#counts_w_1_at.05 <- which(chl_df$tran_counts <=1)
#counts_w_1_at.1<- which(chl_df$tran_counts <=1)
#length(counts_w_1_at.1)

#test <- chl_df[counts_w_1_at.05, ]
#test <- chl_df[counts_w_1_at.1, ]

# make test spatial to see where these point are
library(sf)

# first: # Run spatialize clusters function:
spatialize_clusters <- function(data, coords, filename) {
  spatial_data <- st_as_sf(data, coords = coords, #specify lat lon
                           remove = F, # don't remove these lat/lon cols from df
                           crs = 4326) # add projection
  # write spatial data to creat shp file----
  st_write(spatial_data, paste0("data_output/Spatial_Data/",filename,".shp"), delete_dsn = T)
}

spatialize_clusters(test, c("lon", "lat"), "chl_test_05bin")
spatialize_clusters(binned_flame_F1, c("lon", "lat"), "flame1_05bin")
spatialize_clusters(binned_flame_F, c("lon", "lat"), "flame1_05bin")

#.1 = 458
#.2 = 113
#.25 = 36 records

# OK LEFT OFF HERE: I NEED TO EITHER:
# SUBSAMPLE THE DATA, CLUSTER IT, AND SEE IF THAT IMPACTS CLUSTER RESULTS
# OR INTERPOLATE THE DATA OUT FOR EACH TRANSECT FOR EACH VARIABLE AT A SCALE WHERE THERE IS NO AUTOCORRELATION: THEN DO MEAN CALCULATIONS FROM THE KRIGGING RESULTS?- ASK ANDREW
# OR JUST DECIDE ON HOW MUCH TO BIN BY, THEN LOOK AT VARIOGRAM , THEN THIN DATA TO THAT SCALE, THEN CLUSTER ON THAT AND SEE IF RESULTS CHANGE....

# V1 start here----
# read in csv with all flame data----
flame_all <- read_csv("data/Data_alltransects_withRKM_FIXED_120919.csv") 

flame_all <- flame_all %>% 
  filter(temp >=2)

# Divide each reach into bins .05 km bins (50 m)
# Have to have your column for distance be named "distance" for code to work
flame_all<- data.frame(flame_all)
colnames(flame_all)[1] <- "distance"

# Binning to see how big the bin have to be so there is no bin with <2 transects
binned_flame<- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.1)) 

max(binned_flame$distance)
146/.1

dat <- data.frame(binned_flame, cut(binned_flame$distance, breaks = 146/.1, labels = F))
colnames(dat)[18] <- "cuts"

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

length(which(chl_df$tran_counts ==1)) #458
chl_1_record<- which(chl_df$tran_counts ==1)
 
chl_df_trunc <- chl_df[-chl_1_record, ]

write_csv(chl_df_trunc, "data_output/CHL_raw_means_for_clustering.1bins_NoSingleRecords.csv")
chl_df_trunc <- read_csv("data_output/CHL_raw_means_for_clustering.1bins_NoSingleRecords.csv")

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

nrow(turb_df)
length(which(turb_df$tran_counts ==1)) #458
turb_1_record<- which(turb_df$tran_counts ==1)

turb_df_trunc <- turb_df[-turb_1_record, ]
nrow(turb_df_trunc)

write_csv(turb_df_trunc, "data_output/turb_raw_means_for_clustering.1bins_NoSingleRecords.csv")
turb_df_trunc<- read_csv("data_output/turb_raw_means_for_clustering.1bins_NoSingleRecords.csv")

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

nrow(NO3_df)
length(which(NO3_df$tran_counts ==1)) #458
NO3_1_record<- which(NO3_df$tran_counts ==1)

NO3_df_trunc <- NO3_df[-NO3_1_record, ]
nrow(NO3_df_trunc)

write_csv(NO3_df_trunc, "data_output/NO3_raw_means_for_clustering.1bins_NoSingleRecords.csv")
NO3_df_trunc <- read_csv("data_output/NO3_raw_means_for_clustering.1bins_NoSingleRecords.csv")

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

nrow(pH_df)
length(which(pH_df$tran_counts ==1)) #458
pH_1_record<- which(pH_df$tran_counts ==1)

pH_df_trunc <- pH_df[-pH_1_record, ]
nrow(pH_df_trunc)


write_csv(pH_df_trunc, "data_output/pH_raw_means_for_clustering.1bins_NoSingleRecords.csv")
pH_df_trunc<- read_csv("data_output/pH_raw_means_for_clustering.1bins_NoSingleRecords.csv")

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


nrow(DO_df)
length(which(DO_df$tran_counts ==1)) #458
DO_1_record<- which(DO_df$tran_counts ==1)

DO_df_trunc <- DO_df[-DO_1_record, ]
nrow(DO_df_trunc)


write_csv(DO_df_trunc, "data_output/DO_raw_means_for_clustering.1bins_NoSingleRecords.csv")
DO_df_trunc<- read_csv("data_output/DO_raw_means_for_clustering.1bins_NoSingleRecords.csv")

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


nrow(spCond_df)
length(which(spCond_df$tran_counts ==1)) #458
spCond_1_record<- which(spCond_df$tran_counts ==1)

spCond_df_trunc <- spCond_df[-spCond_1_record, ]
nrow(spCond_df_trunc)

write_csv(spCond_df_trunc, "data_output/spCond_raw_means_for_clustering.1bins_NoSingleRecords.csv")
spCond_df_trunc <- read_csv("data_output/spCond_raw_means_for_clustering.1bins_NoSingleRecords.csv")


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


nrow(fDOM_df)
length(which(fDOM_df$tran_counts ==1)) #458
fDOM_1_record<- which(fDOM_df$tran_counts ==1)

fDOM_df_trunc <- fDOM_df[-fDOM_1_record, ]
nrow(fDOM_df_trunc)

write_csv(fDOM_df_trunc, "data_output/fDOM_raw_means_for_clustering.1bins_NoSingleRecords.csv")
fDOM_df_trunc <- read_csv( "data_output/fDOM_raw_means_for_clustering.1bins_NoSingleRecords.csv")

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


nrow(temp_df)
length(which(temp_df$tran_counts ==1)) #458
temp_1_record<- which(temp_df$tran_counts ==1)

temp_df_trunc <- temp_df[-temp_1_record, ]
nrow(temp_df_trunc)


write_csv(temp_df_trunc, "data_output/temp_raw_means_for_clustering.1bins_NoSingleRecords.csv")
temp_df_trunc<- read_csv("data_output/temp_raw_means_for_clustering.1bins_NoSingleRecords.csv")

# Make one dataframe of the grand mean of the parameters:
#grandMean_allParams <- data.frame(cbind(temp.grand.mean = temp_df$grand_mean,spCond.grand.mean = spCond_df$grand_mean, pH.grand.mean = pH_df$grand_mean, ODO.grand.mean = DO_df$grand_mean, turb.grand.mean = turb_df$grand_mean, fDOM.grand.mean = fDOM_df$grand_mean, CHL.grand.mean = chl_df$grand_mean, NO3.grand.mean = NO3_df$grand_mean), chl_df[, c(2:10, 13)]) # the last part is just binding in the rest of the data that is the same for all params, so I just pulled it from the chl dataframe 

grandMean_allParams_trunc <- data.frame(cbind(temp.grand.mean = temp_df_trunc$grand_mean,spCond.grand.mean = spCond_df_trunc$grand_mean, pH.grand.mean = pH_df_trunc$grand_mean, ODO.grand.mean = DO_df_trunc$grand_mean, turb.grand.mean = turb_df_trunc$grand_mean, fDOM.grand.mean = fDOM_df_trunc$grand_mean, CHL.grand.mean = chl_df_trunc$grand_mean, NO3.grand.mean = NO3_df_trunc$grand_mean), chl_df_trunc[, c(2:10, 13)])

# write to a csv:----
write_csv(grandMean_allParams_trunc, "data_output/grandMeans.fromRaw.allParams.forClustering.1bins_NoSingleRecords.csv")

# Now need to use geomean to get one coord for each bin, use geomean loop:
# Pulled this code from predict.values.from.Gamms.R code, and adapted for this need:
# Need to get the geomean for each bin, and then keep only the one mean value for each variable. 
#Geomean Loop----
#The geomean function find the mean in a set of coordinates.
# Get unique bins <- unique(grandMean_allParams_trunc$cuts)
#
cuts <- unique(grandMean_allParams_trunc$cuts)
length(cuts)
cut_geomean_df <- NULL # make sure to run this before starting the loop 
cuts_coords_trunc <- data.frame(grandMean_allParams_trunc[, c("cuts", "lon", "lat")])

dat_cuts <- cuts_coords_trunc

for (i in cuts) {
  print(i)
  cut <- subset(dat_cuts, dat_cuts$cuts == i)
  coords <- data.frame(cut[ , c("lon", "lat")])
  coords$lon <- as.numeric(as.character(coords$lon))
  coords$lat <- as.numeric(as.character(coords$lat))
  n <- nrow(cut)
  if (nrow(cut) == 1) {
    geomean_i$x <- data.frame(mean(coords$lon))
    geomean_i$y <- data.frame(mean(coords$lat))
    rkm_geomean_df <- rbind(cut_geomean_df, data.frame(unique(cut$cuts), geomean_i, n))
  } else {
    geomean_i <- data.frame(geomean(coords))
    cut_geomean_df <- rbind(cut_geomean_df, data.frame(unique(cut$cuts), geomean_i, n))
  }
}

colnames(cut_geomean_df) <- c("cut", "geomean_lon", "geomean_lat", "cut_n")
cut_geomean_df$geomean_lon <- as.numeric(cut_geomean_df$geomean_lon) # change from list to numeric
cut_geomean_df$geomean_lat <- as.numeric(cut_geomean_df$geomean_lat) # change from list to numeric

write_csv(cut_geomean_df, "data_output/geomean_of_bins_.1_NoSingleRecords")
read_csv("data_output/geomean_of_bins_.1_NoSingleRecords")

# NOw need to incorporate geomeans into data


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

