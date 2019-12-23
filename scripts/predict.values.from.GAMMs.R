# predicting out values for chlorophyll. This might be able to allow for all three transect to be in each bin (except for the last ~50 km)
library(tidyverse)
library(mgcv)
library(geosphere)

#1. read in csv with all flame data----
flame_all <- read_csv("data/Data_alltransects_withRKM_FIXED_120919.csv") 

# CHL data
chl <- flame_all[ , c("CHL", "lat", "lon", "reach", "rkm", "transect" )]
chl$transect <- factor(chl$transect)

# Turb data
turb <- flame_all[ , c("turb","lat", "lon", "reach", "rkm", "transect")]
turb$transect <- factor(turb$transect)

# DO data
DO <- flame_all[ , c("ODO", "rkm","lat", "lon", "reach", "transect")]
DO$transect <- factor(DO$transect)

# fDOM data
fdom <- flame_all[ , c("fDOM", "rkm","lat", "lon", "reach", "transect")]
fdom$transect <- factor(fdom$transect)

# sp Cond data
spCond <- flame_all[ , c("spCond", "rkm", "lat", "lon", "reach", "transect")]
spCond$transect <- factor(spCond$transect)

#NO3 data
NO3 <- flame_all[ , c("NO3", "rkm","lat", "lon", "reach", "transect")]
NO3$transect <- factor(NO3$transect)

# Temp data
temp <- flame_all[ , c("temp", "rkm","lat", "lon", "reach", "transect")]
temp$transect <- factor(temp$transect)


#pH data
pH <- flame_all[ , c("pH", "rkm","lat", "lon", "reach", "transect")]
pH$transect <- factor(pH$transect)

# Determine the rkm to predict out to for each transect:
T1 <- flame_all %>% 
  filter(transect == "T1_")
max(T1$rkm) # for T1 predict out to 111 rkm
min(T1$rkm)

T2 <- flame_all %>% 
  filter(transect == "T2_")
max(T2$rkm) #for T2 predict out to 145 rkm

T3 <- flame_all %>% 
  filter(transect == "T3_")
max(T3$rkm) #for T3 predict out to 145 rkm

data_for_each_param <- flame_all[, c(1, 10:15)]

#2a. New data option 1----
#Make dataframe to predict out to FOR PREDICTING MEASUREMENTS FOR EVERY 0.05 KM FROM 0 TO 145 KM. *****OR TO 2b below ******

new1 <- data.frame(transect = rep(x = "T1_", times = length(seq(0, 111, by = 0.05))), 
                   rkm = seq(0, 111, by = 0.05))
new2 <- data.frame(transect = rep(x = "T2_", times = length(seq(0, 145, by = 0.05))), 
                   rkm = seq(0, 145, by = 0.05))
new3 <- data.frame(transect = rep(x = "T3_", times = length(seq(0, 145, by = 0.05))),
                   rkm = seq(0, 145, by = 0.05))
newd <- rbind(new1, new2, new3)

write_csv(newd, "data_output/newdataRKM_predict_GAMM.csv")
newd<- read_csv("data_output/newdataRKM_predict_GAMM.csv")

# ^ Predicting the data this way, you will evenutally need coords to make the data spatial. This is where the script CalculatingCoords_forPredictedrkmsFINAL comes in

#2b. New data option 2----
# Make dataframe to predict across just the rkms that we CURRENTLY have in our dataset
# Need to get unique rkms to predict across for T1 (stopped earlier than T2 and T3):
T1_rkms_to_predict_across <- which(flame_all$rkm <= 111.0)
T1_rkms <- flame_all$rkm[T1_rkms_to_predict_across]
T1_rkms_u <- data.frame(unique(T1_rkms))
colnames(T1_rkms_u) <- "unique_rkms"
max(T1_rkms_u); min(T1_rkms_u)

T2_rkms_to_predict_across <- which(flame_all$rkm <= 145)
T2_rkms <- flame_all$rkm[T2_rkms_to_predict_across]
T2_rkms_u <- data.frame(unique(T2_rkms))
colnames(T2_rkms_u) <- "unique_rkms"
max(T2_rkms_u); min(T2_rkms_u)

T3_rkms_to_predict_across <- which(flame_all$rkm <= 145)
T3_rkms <- flame_all$rkm[T3_rkms_to_predict_across]
T3_rkms_u <- data.frame(unique(T3_rkms))
colnames(T3_rkms_u) <- "unique_rkms"
max(T3_rkms_u); min(T3_rkms_u)

all_transect_rkms <- rbind(T1_rkms_u,T2_rkms_u, T3_rkms_u )

# then get the unique rkms from the unique rkms for each transect :
all_unique_transect_rkms <- unique(all_transect_rkms)

# set rkms for each transect to predict across:
T1p <- all_unique_transect_rkms %>% 
  filter (all_unique_transect_rkms >= 0.4 & all_unique_transect_rkms <= 111) # bc T1 started a little later and ended much earlier

T2p <- all_unique_transect_rkms

T3p <- all_unique_transect_rkms

# Need to get unique rkms to predict across for T2 and T3:

new1 <- data.frame(transect = rep(x = "T1_", times = nrow(T1p)), 
                                  rkm = T1p$unique_rkms)
new2 <- data.frame(transect = rep(x = "T2_", times = nrow(T2p)), 
                                  rkm = T2p$unique_rkms)
new3 <- data.frame(transect = rep(x = "T3_", times = nrow(T3p)), 
                                  rkm = T3p$unique_rkms)
                   
newd <- rbind(new1, new2, new3)

#test <- merge(newd, flame_all, by = "rkm")

flame_ordered_by_rkm <- flame_all[order(flame_all$rkm), ]
unique_flame_ordered <- unique(flame_ordered_by_rkm[, c("rkm", "lat", "lon")])

rkm_coords_all <- data.frame(flame_all[ , c("rkm", "lon", "lat")]) # make sure class dataframe
class(rkm_coords_all)
rkm_coords_all$lon <- as.numeric(as.character(rkm_coords_all$lon)) # since the coords are factors right now, need to turn to a character and then a numeric 
rkm_coords_all$lat <-  as.numeric(as.character(rkm_coords_all$lat))# # since the coords are factors right now, need to turn to a character and then a numeric
#check:
class(rkm_coords_all$lon)
class(rkm_coords_all$lat)

#Geomean Loop----
rkms <- all_unique_transect_rkms$unique_rkms
rkms
rkm_geomean_df <- NULL # make sure to run this before starting the loop 
dat <- rkm_coords_all 

for (i in rkms) {
 print(i)
  rkm <- subset(dat, dat$rkm == i)
  coords <- data.frame(rkm[ , c("lon", "lat")])
  coords$lon <- as.numeric(as.character(coords$lon))
  coords$lat <- as.numeric(as.character(coords$lat))
  n <- nrow(rkm)
  if (nrow(rkm) == 1) {
    geomean_i$x <- data.frame(mean(coords$lon))
    geomean_i$y <- data.frame(mean(coords$lat))
    rkm_geomean_df <- rbind(rkm_geomean_df, data.frame(unique(rkm$rkm), geomean_i, n))
  } else {
  geomean_i <- data.frame(geomean(coords))
  rkm_geomean_df <- rbind(rkm_geomean_df, data.frame(unique(rkm$rkm), geomean_i, n))
  }
}

ones <- subset(rkm_geomean_df,  rkm_geomean_df$n == 1)
test <- subset(rkm_geomean_df, rkm_geomean_df$unique.rkm.rkm. == rkm_geomean_df[1, 1] )
               
length(which(rkm_geomean_df$y > 1000))
max(rkms)
test <- subset(dat, dat$rkm == rkms[length(rkms)])

test2<- rkm_geomean_df$unique.rkm.rkm. %in% rkms
test3 <- rkms %in% rkm_geomean_df$unique.rkm.rkm.
which(test3 == F)
which(test2 == F)
rkms[test3 == F, ]
rkm_geomean_df[test2==F, ]

length(nrow(test))
# troubleshooting code:

for (i in rkms) {
  print(i)
  i <- 49.4890489
  rkm <- subset(dat, dat$rkm == rkm[1] )
  coords <- data.frame(dat[ , c("lon", "lat")])
  coords$lon <- as.numeric(as.character(rkm_coords_all$lon))
  coords$lat <- as.numeric(as.character(rkm_coords_all$lat))
  geomean_i <- data.frame(geomean(coords))
  rkm_geomean_df <- rbind(rkm_geomean_df, data.frame(unique(rkm$rkm), geomean_i))
}

## just a though- if didnt want to use the geomean, only mybe the first coord, could throw that code in loop instead
rkm_geomean_df$geomean_lon <- as.numeric(rkm_geomean_df$geomean_lon) # change from list to numeric
rkm_geomean_df$geomean_lat <- as.numeric(rkm_geomean_df$geomean_lat) # change from list to numeric
colnames(rkm_geomean_df) <- c("rkm", "geomean_lon", "geomean_lat", "n")


test <- merge(rkm_geomean_df, newd, by = "rkm") # should be able to just use this line of code to put the coords back in after you get your predicted data 
# STOPPED HERE, NEED TO RUN CODE BELOW WITH NEW DATA


#3. GAMMS----
# CHL----
gam_all_chl <- mgcv::gam(CHL ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                         data= chl,
                         family = gaussian,
                         method = "REML")
#pdf(file = "figure_output/chl_gam.pdf")
chl_gam_plot <- plot(y = gam_all_chl$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "CHL-a fitted values", xlab = "Distance from Upper Release (km)")
chl_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_chl)
mgcv::summary.gam(gam_all_chl, re.test = TRUE)
gam.check(gam_all_chl)
mgcv::plot.gam(gam_all_chl)

# Make dataframe of predicted values 
chl_predict <- predict(gam_all_chl, newdata = newd,type="response",se.fit=T)
chl_predicted_df <- cbind(newd, chl_predict)
upper_se <- data.frame(upper_se = chl_predicted_df$fit + chl_predicted_df$se.fit)
lower_se <- data.frame(lower_se = chl_predicted_df$fit - chl_predicted_df$se.fit)
chl_predicted_df <- data.frame(cbind(chl_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(chl_predicted_df,"data_output/Predicted_from_GAMM/chl_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(chl_predicted_df,"data_output/Predicted_from_GAMM/chl_gam_predicted_across_actualrkms.csv")
chl_predicted_df <- read_csv("data_output/Predicted_from_GAMM/chl_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(chl$rkm, chl$CHL, cex = .1)  # raw data 
plot(x = chl_predicted_df$rkm, y = chl_predicted_df$fit, cex = .1)
points(chl_predicted_df$rkm, chl_predicted_df$upper_se, cex = .01, col = "red")
points(chl_predicted_df$rkm, chl_predicted_df$lower_se, cex = .01, col = "red")

# Turb----
gam_all_turb <- mgcv::gam(turb ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                         data= turb,
                         family = gaussian,
                         method = "REML")

#pdf(file = "figure_output/turb_gam.pdf")
turb_gam_plot <- plot(y = gam_all_turb$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "turb fitted values", xlab = "Distance from Upper Release (km)")
turb_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_turb)
mgcv::summary.gam(gam_all_turb, re.test = TRUE)
gam.check(gam_all_turb)
mgcv::plot.gam(gam_all_turb)

# Make dataframe of predicted values 
turb_predict <- predict(gam_all_turb, newdata = newd,type="response",se.fit=T)
turb_predicted_df <- cbind(newd, turb_predict)
upper_se <- data.frame(upper_se = turb_predicted_df$fit + turb_predicted_df$se.fit)
lower_se <- data.frame(lower_se = turb_predicted_df$fit - turb_predicted_df$se.fit)
turb_predicted_df <- data.frame(cbind(turb_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(turb_predicted_df,"data_output/Predicted_from_GAMM/turb_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(turb_predicted_df,"data_output/Predicted_from_GAMM/turb_gam_predicted_across_actualrkms.csv")
turb_predicted_df <- read_csv("data_output/Predicted_from_GAMM/turb_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(turb$rkm, turb$turb, cex = .1)  # raw data 
plot(x = turb_predicted_df$rkm, y = turb_predicted_df$fit, cex = .1)
points(turb_predicted_df$rkm, turb_predicted_df$upper_se, cex = .01, col = "red")
points(turb_predicted_df$rkm, turb_predicted_df$lower_se, cex = .01, col = "red")


# NO3----
gam_all_NO3 <- mgcv::gam(NO3 ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                          data= NO3,
                          family = gaussian,
                          method = "REML")
#pdf(file = "figure_output/NO3_gam.pdf")
NO3_gam_plot <- plot(y = gam_all_NO3$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "NO3 fitted values", xlab = "Distance from Upper Release (km)")
NO3_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_NO3)
mgcv::summary.gam(gam_all_NO3, re.test = TRUE)
gam.check(gam_all_NO3)
mgcv::plot.gam(gam_all_NO3)

# Make dataframe of predicted values 
NO3_predict <- predict(gam_all_NO3, newdata = newd,type="response",se.fit=T)
NO3_predicted_df <- cbind(newd, NO3_predict)
upper_se <- data.frame(upper_se = NO3_predicted_df$fit + NO3_predicted_df$se.fit)
lower_se <- data.frame(lower_se = NO3_predicted_df$fit - NO3_predicted_df$se.fit)
NO3_predicted_df <- data.frame(cbind(NO3_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(NO3_predicted_df,"data_output/Predicted_from_GAMM/NO3_gam_predicted.05km.csv")# made using newd data from option 1 (#2a)
write_csv(NO3_predicted_df,"data_output/Predicted_from_GAMM/NO3_gam_predicted_across_actualrkms.csv")
NO3_predicted_df <- read_csv("data_output/Predicted_from_GAMM/NO3_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(NO3$rkm, NO3$NO3, cex = .1)  # raw data 
plot(x = NO3_predicted_df$rkm, y = NO3_predicted_df$fit, cex = .1)
points(NO3_predicted_df$rkm, NO3_predicted_df$upper_se, cex = .01, col = "red")
points(NO3_predicted_df$rkm, NO3_predicted_df$lower_se, cex = .01, col = "red")

# spCond----
gam_all_spCond <- mgcv::gam(spCond ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                          data= spCond,
                          family = gaussian,
                          method = "REML")
#pdf(file = "figure_output/spCond_gam.pdf")
spCond_gam_plot <- plot(y = gam_all_spCond$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "spCond fitted values", xlab = "Distance from Upper Release (km)")
spCond_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_spCond)
mgcv::summary.gam(gam_all_spCond, re.test = TRUE)
gam.check(gam_all_spCond)
mgcv::plot.gam(gam_all_spCond)

# Make dataframe of predicted values 
spCond_predict <- predict(gam_all_spCond, newdata = newd,type="response",se.fit=T)
spCond_predicted_df <- cbind(newd, spCond_predict)
upper_se <- data.frame(upper_se = spCond_predicted_df$fit + spCond_predicted_df$se.fit)
lower_se <- data.frame(lower_se = spCond_predicted_df$fit - spCond_predicted_df$se.fit)
spCond_predicted_df <- data.frame(cbind(spCond_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(spCond_predicted_df,"data_output/Predicted_from_GAMM/spCond_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(spCond_predicted_df,"data_output/Predicted_from_GAMM/spCond_gam_predicted_across_actualrkms.csv")
spCond_predicted_df <-read_csv("data_output/Predicted_from_GAMM/spCond_gam_predicted_across_actualrkms.csv")
# plot for visual:
plot(spCond$rkm, spCond$spCond, cex = .1)  # raw data 
plot(x = spCond_predicted_df$rkm, y = spCond_predicted_df$fit, cex = .1)
points(spCond_predicted_df$rkm, spCond_predicted_df$upper_se, cex = .01, col = "red")
points(spCond_predicted_df$rkm, spCond_predicted_df$lower_se, cex = .01, col = "red")

# temp----
gam_all_temp <- mgcv::gam(temp ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                          data= temp,
                          family = gaussian,
                          method = "REML")
#pdf(file = "figure_output/temp_gam.pdf")
temp_gam_plot <- plot(y = gam_all_temp$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "temp fitted values", xlab = "Distance from Upper Release (km)")
temp_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_temp)
mgcv::summary.gam(gam_all_temp, re.test = TRUE)
gam.check(gam_all_temp)
mgcv::plot.gam(gam_all_temp)

# Make dataframe of predicted values 
temp_predict <- predict(gam_all_temp, newdata = newd,type="response",se.fit=T)
temp_predicted_df <- cbind(newd, temp_predict)
upper_se <- data.frame(upper_se = temp_predicted_df$fit + temp_predicted_df$se.fit)
lower_se <- data.frame(lower_se = temp_predicted_df$fit - temp_predicted_df$se.fit)
temp_predicted_df <- data.frame(cbind(temp_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(temp_predicted_df,"data_output/Predicted_from_GAMM/temp_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(temp_predicted_df,"data_output/Predicted_from_GAMM/temp_gam_predicted_across_actualrkms.csv")
temp_predicted_df <- read_csv("data_output/Predicted_from_GAMM/temp_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(temp$rkm, temp$temp, cex = .1)  # raw data 
plot(x = temp_predicted_df$rkm, y = temp_predicted_df$fit, cex = .1)
points(temp_predicted_df$rkm, temp_predicted_df$upper_se, cex = .01, col = "red")
points(temp_predicted_df$rkm, temp_predicted_df$lower_se, cex = .01, col = "red")


# pH----
gam_all_pH <- mgcv::gam(pH ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                          data= pH,
                          family = gaussian,
                          method = "REML")
#pdf(file = "figure_output/pH_gam.pdf")
pH_gam_plot <- plot(y = gam_all_pH$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "pH fitted values", xlab = "Distance from Upper Release (km)")
pH_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_pH)
mgcv::summary.gam(gam_all_pH, re.test = TRUE)
gam.check(gam_all_pH)
mgcv::plot.gam(gam_all_pH)

# Make dataframe of predicted values 
pH_predict <- predict(gam_all_pH, newdata = newd,type="response",se.fit=T)
pH_predicted_df <- cbind(newd, pH_predict)
upper_se <- data.frame(upper_se = pH_predicted_df$fit + pH_predicted_df$se.fit)
lower_se <- data.frame(lower_se = pH_predicted_df$fit - pH_predicted_df$se.fit)
pH_predicted_df <- data.frame(cbind(pH_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
# write_csv(pH_predicted_df,"data_output/Predicted_from_GAMM/pH_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(pH_predicted_df,"data_output/Predicted_from_GAMM/pH_gam_predicted_across_actualrkms.csv")
pH_predicted_df <- read_csv("data_output/Predicted_from_GAMM/pH_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(pH$rkm, pH$pH, cex = .1)  # raw data 
plot(x = pH_predicted_df$rkm, y = pH_predicted_df$fit, cex = .1)
points(pH_predicted_df$rkm, pH_predicted_df$upper_se, cex = .01, col = "red")
points(pH_predicted_df$rkm, pH_predicted_df$lower_se, cex = .01, col = "red")

# fDOM----
gam_all_fDOM <- mgcv::gam(fDOM ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                          data= fdom,
                          family = gaussian,
                          method = "REML")
#pdf(file = "figure_output/fDOM_gam.pdf")
fDOM_gam_plot <- plot(y = gam_all_fDOM$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "fDOM fitted values", xlab = "Distance from Upper Release (km)")
fDOM_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_fDOM)
mgcv::summary.gam(gam_all_fDOM, re.test = TRUE)
gam.check(gam_all_fDOM)
mgcv::plot.gam(gam_all_fDOM)

# Make dataframe of predicted values 
fDOM_predict <- predict(gam_all_fDOM, newdata = newd,type="response",se.fit=T)
fDOM_predicted_df <- cbind(newd, fDOM_predict)
upper_se <- data.frame(upper_se = fDOM_predicted_df$fit + fDOM_predicted_df$se.fit)
lower_se <- data.frame(lower_se = fDOM_predicted_df$fit - fDOM_predicted_df$se.fit)
fDOM_predicted_df <- data.frame(cbind(fDOM_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(fDOM_predicted_df,"data_output/Predicted_from_GAMM/fDOM_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(fDOM_predicted_df,"data_output/Predicted_from_GAMM/fDOM_gam_predicted_across_actualrkms.csv")
fDOM_predicted_df<- read_csv("data_output/Predicted_from_GAMM/fDOM_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(fdom$rkm, fdom$fDOM, cex = .1)  # raw data 
plot(x = fDOM_predicted_df$rkm, y = fDOM_predicted_df$fit, cex = .1)
points(fDOM_predicted_df$rkm, fDOM_predicted_df$upper_se, cex = .01, col = "red")
points(fDOM_predicted_df$rkm, fDOM_predicted_df$lower_se, cex = .01, col = "red")

# DO----
gam_all_DO <- mgcv::gam(ODO ~ s(rkm, bs="fs", by=transect, k=20) + s(transect, bs="re"),
                          data= DO,
                          family = gaussian,
                          method = "REML")
#pdf(file = "figure_output/DO_gam.pdf")
DO_gam_plot <- plot(y = gam_all_DO$fitted.values, x = data_for_each_param$rkm, pch = 16, cex = .7, ylab = "DO fitted values", xlab = "Distance from Upper Release (km)")
DO_gam_plot
#dev.off()
# Checking GAMM:
mgcv::gam.vcomp(gam_all_DO)
mgcv::summary.gam(gam_all_DO, re.test = TRUE)
gam.check(gam_all_DO)
mgcv::plot.gam(gam_all_DO)

# Make dataframe of predicted values 
DO_predict <- predict(gam_all_DO, newdata = newd,type="response",se.fit=T)
DO_predicted_df <- cbind(newd, DO_predict)
upper_se <- data.frame(upper_se = DO_predicted_df$fit + DO_predicted_df$se.fit)
lower_se <- data.frame(lower_se = DO_predicted_df$fit - DO_predicted_df$se.fit)
DO_predicted_df <- data.frame(cbind(DO_predicted_df, upper_se, lower_se))

# Write to a csv:
# This dataframe has ONE predicted value per rkm (binned by 0.05 rkm)
#write_csv(DO_predicted_df,"data_output/Predicted_from_GAMM/DO_gam_predicted.05km.csv") # made using newd data from option 1 (#2a)
write_csv(DO_predicted_df,"data_output/Predicted_from_GAMM/DO_gam_predicted_across_actualrkms.csv")
DO_predicted_df<- read_csv("data_output/Predicted_from_GAMM/DO_gam_predicted_across_actualrkms.csv")

# plot for visual:
plot(DO$rkm, DO$ODO, cex = .1)  # raw data 
plot(x = DO_predicted_df$rkm, y = DO_predicted_df$fit, cex = .1)
points(DO_predicted_df$rkm, DO_predicted_df$upper_se, cex = .01, col = "red")
points(DO_predicted_df$rkm, DO_predicted_df$lower_se, cex = .01, col = "red")

#Combine into 1 data frame:----
#predicted_gamm_AllParams <- data.frame(cbind(CHL.predicted = chl_predicted_df$fit,turb.predicted =  turb_predicted_df$fit, NO3.predicted = NO3_predicted_df$fit, spCond.predicted = spCond_predicted_df$fit, temp.predicted = temp_predicted_df$fit, pH.predicted = pH_predicted_df$fit, fDOM.predicted = fDOM_predicted_df$fit, DO.predicted = DO_predicted_df$fit, chl_predicted_df[, 1:2])) #### USE IF USED newd DATA FROM 2a

predicted_gamm_AllParams <- data.frame(cbind(CHL.predicted = chl_predicted_df$fit,turb.predicted =  turb_predicted_df$fit, NO3.predicted = NO3_predicted_df$fit, spCond.predicted = spCond_predicted_df$fit, temp.predicted = temp_predicted_df$fit, pH.predicted = pH_predicted_df$fit, fDOM.predicted = fDOM_predicted_df$fit, DO.predicted = DO_predicted_df$fit), chl_predicted_df[, 1:2]) ### USE IF USED newd DATA FROM 2b

write_csv(predicted_gamm_AllParams, "data_output/CV/CV_fromPredicted_AllParams_actualRKMS_NOcoordsYet.csv")

# MERGE coords into dataset ----
# merges the predicted dataframe with geomean dataframe made in step 2b:
predicted_gamm_AllParamsFinal <- merge(x = predicted_gamm_AllParams, y = rkm_geomean_df, by = "rkm")

#write to csv:
#write_csv(predicted_gamm_AllParams, "data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams.csv") ## made using newd data from option 1 (#2a)

write_csv(predicted_gamm_AllParamsFinal, "data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams_across_actualRKMs_withGeomeans.csv") ## made using newd data from option 2 (#2b)

# WORKFLOW: now read in this final csv output into the Modeling_Variation_forTimeAnalysis.R script to calculate the coefficients of variation 