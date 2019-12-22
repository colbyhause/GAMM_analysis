# The purpose of this script is to solve this issue: When I predicted out values for each parameter for every 0.05 km interval I lost my lat lon coordinates in the process.
# ** Just to remind myself of the REASON why I decided to used predicted values to calculate CV:  because in the raw data some of the .05 km bins only had data from 1 transect in them, so I couldn't do a calculation without either dropping those data or making the bin .3 km. Predicting it out allowed meto have 1 estimate for each transect in each 0.05 km bin. Now I am realizing what I could have done instead was only predicted across the rkms that I actually had to avoid having to "estimate" coords for exatly every 0.05 km from 0 to 145 km, as the code below does. Oh well, Ill try it both ways i guess? ** 
# To make the cluster analysis on these data spatial, I need those coordinates. So there were 2 steps to this, first I had to figure out a way to assign coordinates to the predicted rkm values, which were every .05 km. To do this, I rounded each of the actual rkm values (calcluated by the SJ_River_Distabce script) to the nearest 0.05. Then I took those and found the "geomean" ie. the average of every bin of coords and that became the coord for that bin.
#  Below is the code for the "nearest" fucntion, which allows you to round a number to a specified nearest value. Again, I wanted to round the rkms from my raw data to the nearest 0.05 value so i can then pull the lat lon values for those


# This script includes code that I got off stackoverflow:
# https://stackoverflow.com/questions/38985186/how-to-round-to-0-25-or-0-75-but-not-0-00-or-0-50-in-r

FTOL <- 1e-8;
feq <- function(a,b,tol=FTOL) ifelse(is.finite(a) & is.finite(b),abs(a-b)<=max(abs(a),abs(b))*tol,a==b);
fne <- function(a,b,tol=FTOL) ifelse(is.finite(a) & is.finite(b),abs(a-b)>max(abs(a),abs(b))*tol,a!=b);
fgt <- function(a,b,tol=FTOL) ifelse(is.finite(a) & is.finite(b),a-b>max(abs(a),abs(b))*tol,a>b);
fge <- function(a,b,tol=FTOL) ifelse(is.finite(a) & is.finite(b),a-b>=-max(abs(a),abs(b))*tol,a>=b);
flt <- function(a,b,tol=FTOL) ifelse(is.finite(a) & is.finite(b),b-a>max(abs(a),abs(b))*tol,a<b);
fle <- function(a,b,tol=FTOL) ifelse(is.finite(a) & is.finite(b),b-a>=-max(abs(a),abs(b))*tol,a<=b);

HALFRULE_UP   <- 1L; ## round towards +Inf
HALFRULE_DOWN <- 2L; ## round towards -Inf
HALFRULE_IN   <- 3L; ## round towards zero
HALFRULE_OUT  <- 4L; ## round away from zero
HALFRULE_EVEN <- 5L; ## round to the even multiple of the two multiples of nearest that are tied
HALFRULE_ODD  <- 6L; ## round to the odd multiple of the two multiples of nearest that are tied
nearest <- function(x,nearest=1,offset=0,halfrule=HALFRULE_EVEN) {
  
  ## ensure correct types
  x <- as.double(x);
  nearest <- as.double(nearest);
  offset <- as.double(offset);
  halfrule <- as.integer(halfrule);
  
  ## validate
  v <- which(!halfrule%in%c(HALFRULE_UP,HALFRULE_DOWN,HALFRULE_IN,HALFRULE_OUT,HALFRULE_EVEN,HALFRULE_ODD)); if (length(v)>0L) stop(paste0('invalid halfrule: ',halfrule[v[1L]],'.'));
  
  ## normalize lengths
  len <- max(length(x),length(nearest),length(halfrule));
  x <- rep(x,len=len);
  nearest <- rep(nearest,len=len);
  offset <- rep(offset,len=len);
  halfrule <- rep(halfrule,len=len);
  
  ## apply offset
  x <- x-offset;
  
  ## must treat zero nearests different from non-zero
  nonzero <- fne(nearest,0);
  
  ## start building result
  res <- double(length(x));
  
  ## nearest zero doesn't really make any sense; but you can consider every possible number to be at the nearest zero
  res[!nonzero] <- x[!nonzero];
  
  ## simplify subsequent operations to only focus on non-zero nearest
  x <- x[nonzero];
  nearest <- nearest[nonzero];
  halfrule <- halfrule[nonzero];
  offset <- offset[nonzero];
  
  ## don't care about negativity of nearest
  nearest <- abs(nearest);
  
  ## get how many times nearest goes into x, truncated
  mult <- as.integer(x/nearest); ## note: can't use %/%, since that actually floors toward -Inf
  
  ## get round-truncated result
  r.trunc <- mult*nearest;
  
  ## get absolute excess over r.trunc
  excess <- abs(x - r.trunc);
  
  ## get half of nearest
  half.of.nearest <- nearest*0.5;
  
  ## add one to mult if necessary; whether we need to do this in the case of a tie depends on the specified tiebreaker rule and which side of the zero multiple x is on
  adds <- which(
    fgt(excess,half.of.nearest)
    | feq(excess,half.of.nearest) & (
      halfrule==HALFRULE_UP & fgt(x,0)
      | halfrule==HALFRULE_DOWN & flt(x,0)
      | halfrule==HALFRULE_OUT
      | halfrule==HALFRULE_EVEN & mult%%2L!=0L
      | halfrule==HALFRULE_ODD & mult%%2L==0L
    )
  );
  mult[adds] <- mult[adds] + ifelse(flt(x[adds],0),-1,1);
  
  ## recover target value from mult, and at the same time unshift offset
  res[nonzero] <- nearest*mult+offset;
  
  res;
  
}; ## end nearest()
nearest.halfup   <- function(x,nearest=1,offset=0) nearest(x,nearest,offset,HALFRULE_UP  );
nearest.halfdown <- function(x,nearest=1,offset=0) nearest(x,nearest,offset,HALFRULE_DOWN);
nearest.halfin   <- function(x,nearest=1,offset=0) nearest(x,nearest,offset,HALFRULE_IN  );
nearest.halfout  <- function(x,nearest=1,offset=0) nearest(x,nearest,offset,HALFRULE_OUT );
nearest.halfeven <- function(x,nearest=1,offset=0) nearest(x,nearest,offset,HALFRULE_EVEN);
nearest.halfodd  <- function(x,nearest=1,offset=0) nearest(x,nearest,offset,HALFRULE_ODD );

test <- nearest(flame_all$rkm,0.05) # this works

# Add the rounded rkms from the test df to the data :
data_nearest_rkm<- cbind(test, flame_all)

# take a subset of this to test geomean function
nearest_subset <- data_nearest_rkm %>% 
  filter(test == 0.55)

# Determine average coordinate of the subset using geomean function from geosphere package
#install.packages("geosphere")
library(geosphere)

# First I thought that yuou had to make the coords into spatial datapoints, but thats not true. The coordinates just have to be numeric. 
#spatpoints_2<- SpatialPoints(flame_all[,10:11], proj4string=CRS(as.character("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")), bbox = NULL)

# make into matrix:
nearest_subset<- as.matrix(nearest_subset) 
nearest_subset_coords<- nearest_subset[, c("lon", "lat")]
nearest_subset_coords <- data.frame(nearest_subset_coords) # make sure class dataframe
class(nearest_subset)
nearest_subset_coords$lon <- as.numeric(as.character(nearest_subset_coords$lon)) # since the coords are factors right now, need to turn to a character and then a numeric 
nearest_subset_coords$lat <-  as.numeric(as.character(nearest_subset_coords$lat))# # since the coords are factors right now, need to turn to a character and then a numeric
#check:
class(nearest_subset_coords)
class(nearest_subset_coords$lon)
class(nearest_subset_coords$lat)

mean.55km <- geomean(nearest_subset_coords) # it works!

testrkms <- seq(1, 145, 0.05)
i<- testrkms %in% data_nearest_rkm$test
length(testrkms[i == F])

# now make the subset spatial and then that one mean spatial and see what it looks like in QGIS

nearest_subset_spatial <- st_as_sf(nearest_subset_coords, 
                                 coords = c("lon", "lat"), # for point data
                                 remove = F, # don't remove these lat/lon cols from df
                                 crs = 4269) # add projection (this is WGS84)

st_write(nearest_subset_spatial, "data_output/Spatial_data/test_subset.55_to_calc_meanCoords.shp")

# make the man point spatial 
colnames(mean.55km ) <- c("lon", "lat")
mean.55km <- data.frame(mean.55km)
geomean_spatial <- st_as_sf(mean.55km, 
                                 coords = c("lon", "lat"), # for point data
                                 remove = F, # don't remove these lat/lon cols from df
                                 crs = 4269) # add projection (this is WGS84)

st_write(geomean_spatial, "data_output/Spatial_data/geomean.55_calc_coords.shp")


# Looks good in QGIS
# Now do this for the whole dataset:----
nearest.05 <- nearest(flame_all$rkm,0.05)

# Add the rounded rkms from the test df to the data :
data_nearest_rkm<- cbind(nearest.05, flame_all)

# Determine average coordinate of the subset using geomean function from geosphere package
#install.packages("geosphere")
library(geosphere)

# make into matrix:
data_nearest_rkm<- as.matrix(data_nearest_rkm) 
#nearest_subset_coords<- nearest_subset[, c("lon", "lat")]
data_nearest_rkm_df <- data.frame(data_nearest_rkm) # make sure class dataframe
class(data_nearest_rkm_df)
data_nearest_rkm_df$lon <- as.numeric(as.character(data_nearest_rkm_df$lon)) # since the coords are factors right now, need to turn to a character and then a numeric 
data_nearest_rkm_df$lat <-  as.numeric(as.character(data_nearest_rkm_df$lat))# # since the coords are factors right now, need to turn to a character and then a numeric
#check:
class(data_nearest_rkm_df)
class(data_nearest_rkm_df$lon)
class(data_nearest_rkm_df$lat)

### NOW NEED TO WRITE LOOP TO OUTPUT A DATAFRAME WITH THE GEOMEAN FOR EACH RKM IN UNIQUE(DATA_NEAREST_RKM_DF)
mean.05km <- geomean(nearest_subset_coords) # this needas to go in loop ^

# just seeing how many rkms I am going to lose when i assign coords back to the predicted dataframe :
testrkms <- seq(1, 145, 0.05)
i<- testrkms %in% data_nearest_rkm$test
length(testrkms[i == F])

# now make the subset spatial and then that one mean spatial and see what it looks like in QGIS

nearest_subset_spatial <- st_as_sf(nearest_subset_coords, 
                                   coords = c("lon", "lat"), # for point data
                                   remove = F, # don't remove these lat/lon cols from df
                                   crs = 4269) # add projection (this is WGS84)

st_write(nearest_subset_spatial, "data_output/Spatial_data/test_subset.55_to_calc_meanCoords.shp")

# make the man point spatial 
colnames(mean.55km ) <- c("lon", "lat")
mean.55km <- data.frame(mean.55km)
geomean_spatial <- st_as_sf(mean.55km, 
                            coords = c("lon", "lat"), # for point data
                            remove = F, # don't remove these lat/lon cols from df
                            crs = 4269) # add projection (this is WGS84)

st_write(geomean_spatial, "data_output/Spatial_data/geomean.55_calc_coords.shp")



