# This script includes code that I got off stackoverflow:
# https://stackoverflow.com/questions/38985186/how-to-round-to-0-25-or-0-75-but-not-0-00-or-0-50-in-r
# My issue was when I predicted out values for each parameter for every 0.05 km interval I lost my lat lon coordinates in the process. To make the cluster analysis on these data spatial, I need those coordinates. Below is the code for the "nearest" fucntion, which allows you to round a number to a specified nearest value. I want to round the rkms from my raw data to the nearest 0.05 value so i can then pull the lat lon values for those 
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

test <- nearest(flame_all$rkm,0.05)

# Make spatial points object:
spatpoints_2<-SpatialPoints(flame_all[,10:11], proj4string=CRS(as.character("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")), bbox = NULL)

geomean(spatpoints_2)
test2_subset<- test2_subset[, c("lon", "lat")]

mean.4kmtest <- geomean(test2_subset) # this worked!

class(test2_subset)
class(test2_subset$lon)
class(test2_subset$lat)

SpatialPointsDataFrame(coords, data, coords.nrs = numeric(0), 
                       proj4string = CRS(as.character(NA)), match.ID, bbox = NULL)


line2network(path = "data_output/nhdplus_flameLine.gpkg", layer="NHDFlowline_Network", reproject = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

#########THIS IS WHERE I LEFT OFF: ^^^^
#########issue: description above. I left off trying to figure out how to make coords into a spatial point dataframe. the cvoords i was putting into geomean are working, and i think they need to be spatial points. use the river dist code to getb the projections. after that try to get the mean for each subset of points, and then convert that back to coords. 

test2<- cbind(test, flame_all)
test2_subset <- test2 %>% 
  filter(test == 0.55)
mean(test2_subset$lat) 37.31174
mean(test2_subset$lon)  -120.935
# ^^^ this method might actually work, look at google earth page- i plotted the ream against 3 of the raw points and it does a pretty good job at hitting the center. Maybe do this for all and then plot all of the in qgis 




install.packages("geosphere")
library(geosphere)
# make into matrix:
test2_subset <- as.matrix(test2_subset)
mean.4km <- geomean(test2_subset[, c("lon", "lat")])

class(test2_subset$lon)
class(test2_subset$lat)
class(test2_subset)

