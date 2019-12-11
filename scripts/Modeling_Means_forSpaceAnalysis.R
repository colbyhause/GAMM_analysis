# Trying out the "means" model in 2 ways:
# 1. running a GAM on each transect for each parameter, then using those values to bin and get mean 
# 2. just getting means form raw data

library(tidyverse)
#install.packages("Distance")
library(Distance)



#1. Means model using coeff of variation:-----
# read in csv with all flame data
flame_all <- read_csv("data/Data_alltransects_withRKM_FIXED_120919.csv") 
flame_all <- flame_all[, -c(1:2)] # remove first 2 cols, done need them

# Divide each reach into bins .05 km bins (50 m)
# Have to have your column for distance be named "distance" for code to work
flame_all<- data.frame(flame_all)
colnames(flame_all)[1] <- "distance"

# create dataframe with the bins
binned_flame <- create.bins(data = flame_all, cutpoints = seq(from = 0, to = max(flame_all$distance +1), by = 0.05)) # the +1 is to prevent the code from cutting of the distance prematurely 

# Code for checking to see which ones the code was dropping before I added the +1 to it:
#determine_dropped_data <- flame_all$transect_ID %in% binned_flame$transect_ID
#removed_in_binning <- flame_all[determine_dropped_data == F, ]

# look to see if there are enough detections in all the bins, was working and now isnt for some reason:
#ggplot(data = binned_flame) +
  #hist(aes(x = binned_flame$distbegin)) #, breaks = seq(from = 0, to = max(rkm), by = 0.1))
#binned_flame$distbegin <- as.numeric(binned_flame$distbegin)

# calculate the mean for each transect within a bin, and then the global mean for that bin 
test<- binned_flame %>% 
filter(distbegin == 0.0) %>% 
  group_by(transect) %>% 
  mutate(mean_chl_transect = mean(CHL)) %>% 
  group_by(distbegin)  %>% 
  mutate(grand_mean = mean(unique(mean_chl_transect)))

test2<-binned_flame %>% 
  filter(distbegin == 0.15)
# need to turn this into a loop, make an example dataframe:
# 
#### WRITING LOOP FOR THIS:--------------------------------------------------
ex_df <- binned_flame %>% 
  filter(distbegin == 0.0 | distbegin == 0.05 | distbegin == 0.10 | distbegin == 0.15 | distbegin == 0.2| distbegin == 0.25| distbegin == 0.3 | distbegin == 0.65) # | distbegin == 0.35| distbegin == 0.4)# |distbegin == 0.45 | distbegin == 0.5| distbegin == 0.55| distbegin ==0.6)#| distbegin ==0.65)

data <- ex_df
param_list <- c("temp", "spCond", "pH","ODO", "turb", "fDOM", "CHL", "NO3")
param_mean_transect_df <- NULL

# Run first Loop:
 for (i in param_list) {
  print(i)
  param_df <- data[ , c(i, "distbegin", "transect")]
  param_df$newcolumn1 <- 0
  param_df$newcolumn2 <- 0
  start_bins <-  unique(param_df$distbegin)
  for (j in start_bins) {
    print(j)
    param_mean_transect <- param_df
    bin_df <- subset(param_mean_transect, param_mean_transect$distbegin == j)
    colnames(bin_df)[4] <- paste0("mean", colnames(param_df)[1], "byTransect")
    transect_list <- unique(param_mean_transect$transect)
    for (k in transect_list) {
      print(k) 
      param_mean_transect_temporary <- subset(bin_df, bin_df$transect == k)
      print("here")
      param_mean_transect_temporary[4] <- mean(param_mean_transect_temporary[, 1])
      print("here2")
      param_mean_transect_df <- rbind(param_mean_transect_df, param_mean_transect_temporary)
    }
    return(param_mean_transect_df)
  }
 }

# Use param_mean_transect_df from first loop in second loop:
for (j in start_bins ) {
    bin_df_globalmean_temporary <- subset(param_mean_transect_df, 
                                          param_mean_transect_df$distbegin == j)
    bin_df_globalmean_temporary[5] <- mean(unique(bin_df_globalmean_temporary[, 4]))
    colnames(bin_df_globalmean_temporary)[5] <- paste0(colnames(param_df)[1], "_global_mean")
    globalmean_df <- rbind(globalmean_df, data.frame(bin_df_globalmean_temporary))
}
globalmean_df <-NULL

##################


length(param_mean_transect_df$distbegin == 0.0)
rm(globalmean_df)
    #return((param_mean_transect_df))
    print(param_mean_transect_df)
    for (j in start_bins) {
      bin_df_globalmean_temporary <- subset(param_mean_transect_df, 
                                            param_mean_transect_df$distbegin == j)
      bin_df_globalmean <- bin_df_globalmean_temporary
      bin_df_globalmean[5] <-"test"
      #bin_df_globalmean[5] <- mean(bin_df_globalmean_temporary[, 4])
      bin_df_globalmean_df <- rbind(bin_df_globalmean_df, data.frame(bin_df_globalmean))
    }
  }
  return((bin_df_globalmean))
}
  #}
  return((bin_df_globalmean))
 }
    
param_mean_transect_df <- NULL
globalmean_df <- NULL
rm(bin_df_globalmean_temporary)
rm(bin_df_globalmean)
rm(bin_df)

test<- ex_df %>% 
  filter(distbegin == 0.2)

  param_mean_transect_df <- rbind(param_mean_transect_df, data.frame(param_mean_transect_temporary))
      #write_csv(x = param_mean_transect_df, path = paste0(param, "_means_for_clustering.csv"))
    }
  }
 }
  for (j in start_bins) {
    bin_df_globalmean <- subset(param_mean_transect_df, param_mean_transect_df$distbegin == j)
    param_mean_transect_df[5] <- mean(unique(param_mean_transect_df[, 4]))
    param_mean_transect_df_final <- rbind(param_mean_transect_df_final, 
                                          data.frame(param_mean_transect_df))
    
    }
  }
  return(param_mean_transect_df_final)
 }

  #write_csv(x = param_mean_transect_df, path = paste0(param, "_means_for_clustering.csv"))
 
mean(unique(param_mean_transect_df$newcolumn1))
param_mean_transect_df <- NULL
param_mean_transect_df_final <- NULL


param_mean_transect <- subset(param_mean_transect, param_mean_transect$transect == "T1_")

param_mean_transect$newcolumn1 <- mean(param_mean_transect[, 1])
test <- param_mean_transect 
test1<- test %>% 
  mutate(paste0("mean_", colnames(param_df)[1], "by_transect") == "test")

colnames(param_mean_transect)[4] <-  paste0("mean_", colnames(param_df)[1], "by_transect")
  mutate("mean_", param_df$1, "by_transect")

   df<-ex_df %>% 
     filter(distbegin == i) %>% 
     group_by(transect) 
 for (j in param_list) {
  print(j)
   df2 <- df

  
# Means model using GAMMs
# read in csv with all flame data
flame_all <- read_csv("data/Data_alltransects_withRKM.csv") 
flame_all <- flame_all[, -1]
# Remove the one weird temperature outlier:
flame_all <- flame_all %>% 
  filter(temp >=2)


# CHL data
chl <- flame_all[ , c("CHL", "rkm", "transect")]
chl$transect <- factor(chl$transect)

# Turb data
turb <- flame_all[ , c("turb", "rkm", "transect")]
turb$transect <- factor(turb$transect)

# DO data
DO <- flame_all[ , c("ODO", "rkm", "transect")]
DO$transect <- factor(DO$transect)

# fDOM data
fdom <- flame_all[ , c("fDOM", "rkm", "transect")]
fdom$transect <- factor(fdom$transect)

# sp Cond data
spCond <- flame_all[ , c("spCond", "rkm", "transect")]
spCond$transect <- factor(spCond$transect)

#NO3 data
NO3 <- flame_all[ , c("NO3", "rkm", "transect")]
NO3$transect <- factor(NO3$transect)

# Temp data
temp <- flame_all[ , c("temp", "rkm", "transect")]
temp$transect <- factor(temp$transect)


#pH data
pH <- flame_all[ , c("pH", "rkm", "transect")]
pH$transect <- factor(pH$transect)

data_for_each_param <- flame_all[, c(1, 12, 13, 15)]
