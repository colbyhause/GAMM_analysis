# predicting out values for chlorophyll. This might be able to allow for all three transect to be in each bin (except for the last ~50 km)

library(mgcv)
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
newd <- data.frame(flame_all[ , c("transect", "rkm")])
class(newd)


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

# plot for visual:
plot(DO$rkm, DO$ODO, cex = .1)  # raw data 
plot(x = DO_predicted_df$rkm, y = DO_predicted_df$fit, cex = .1)
points(DO_predicted_df$rkm, DO_predicted_df$upper_se, cex = .01, col = "red")
points(DO_predicted_df$rkm, DO_predicted_df$lower_se, cex = .01, col = "red")

#Combine into 1 data frame:----
#predicted_gamm_AllParams <- data.frame(cbind(CHL.predicted = chl_predicted_df$fit,turb.predicted =  turb_predicted_df$fit, NO3.predicted = NO3_predicted_df$fit, spCond.predicted = spCond_predicted_df$fit, temp.predicted = temp_predicted_df$fit, pH.predicted = pH_predicted_df$fit, fDOM.predicted = fDOM_predicted_df$fit, DO.predicted = DO_predicted_df$fit, chl_predicted_df[, 1:2])) #### USE IF USED newd DATA FROM 2a

predicted_gamm_AllParams <- data.frame(cbind(CHL.predicted = chl_predicted_df$fit,turb.predicted =  turb_predicted_df$fit, NO3.predicted = NO3_predicted_df$fit, spCond.predicted = spCond_predicted_df$fit, temp.predicted = temp_predicted_df$fit, pH.predicted = pH_predicted_df$fit, fDOM.predicted = fDOM_predicted_df$fit, DO.predicted = DO_predicted_df$fit, data_for_each_param)) ### USE IF USED newd DATA FROM 2b

#write to csv:
#write_csv(predicted_gamm_AllParams, "data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams.csv") ## made using newd data from option 1 (#2a)

write_csv(predicted_gamm_AllParams, "data_output/Predicted_from_GAMM/GAMMPredictedValues_AllParams_across_actualRKMs.csv") ## made using newd data from option 2 (#2b)
