#' Analyses of ROI-level MODIS-MCD14ML fire risk factors 
#' Example ROI: Tanzania 
#' Source data courtesy UMD (ftp://fuoco.geog.umd.edu/modis/C5/mcd14ml, login=fire, pwd=burnt)
#' M. Walsh & J. Chen, May 2015

# Required packages
# install.packages(c("downloader","rgdal","maptools","plyr","raster")), dependencies=TRUE)
require(downloader)
require(rgdal)
require(maptools)
require(plyr)
require(raster)
require(gstat)

#+ Data downloads ---------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("MCD14ML", showWarnings=F)
dat_dir <- "./MCD14ML"

# download Africa-wide MCD14ML data (note large download: ~490 Mb)
download("https://www.dropbox.com/s/oi0nvimjawmz7xf/AF_MCD14ML.txt.zip?dl=0", "./MCG14ML/AF_MCD14ML.txt.zip", mode="wb")
unzip("./MCD14ML/AF_MCD14ML.txt.zip", exdir="./MCD14ML", overwrite=T)
fires <- read.table(paste(dat_dir, "/AF_MCD14ML.txt", sep=""), header=T)

# download ROI shape file
download("https://www.dropbox.com/s/bb5277wvuce0b2h/TZ_bound.zip?dl=0", "./MCD14ML/TZ_bound.zip", mode="wb")
unzip("./MCD14ML/TZ_bound.zip", exdir="./MCD14ML", overwrite=T)
bound <- readShapeSpatial("./MCD14ML/TZ_bound.shp")

# download GeoSurvey prediction grids and stack in raster
download("https://www.dropbox.com/s/e1hsnz7scc3qf9n/TZ_GS_preds.zip?dl=0", "./MCD14ML/TZ_GS_preds.zip", mode="wb")
unzip("./MCD14ML/TZ_GS_preds.zip", exdir="./MCD14ML", overwrite=T)
glist <- list.files(path="./MCD14ML", pattern="tif", full.names=T)
grids <- stack(glist)

#+ Data setup -------------------------------------------------------------
# Variable subset selection
fires <- fires[c(1,4:7,9:10)] ## variable (column) subset ... drop unnecessary variables
hc_fires <- fires[which(fires$conf > 80), ] ## high confidence (>80%) fire detection subset

# ROI bounding box subset selection
bbx <- bbox(bound) ## ROI bounding box
bbx_fires <- hc_fires[which(hc_fires$lon >= bbx[1,1] & hc_fires$lon <= bbx[1,2] & hc_fires$lat >= bbx[2,1] & hc_fires$lat <= bbx[2,2]), ]

# Project fire observation coords to Africa LAEA CRS
bbx_proj <- as.data.frame(project(cbind(bbx_fires$lon, bbx_fires$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(bbx_proj) <- c("x","y")
bbx_fires <- cbind(bbx_fires, bbx_proj)

# Generate AfSIS grid cell ID's (GID) & dates
res.pixel <- 1000
xgid <- ceiling(abs(bbx_fires$x)/res.pixel)
ygid <- ceiling(abs(bbx_fires$y)/res.pixel)
gidx <- ifelse(bbx_fires$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(bbx_fires$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
bbx_fires <- cbind(bbx_fires, GID)

# Aggregate bounding box-wide observations by date
DID_fires <- ddply(bbx_fires, .(YYYYMMDD), summarize,
                   N = length(YYYYMMDD),
                   FRP = mean(FRP))
DID_fires$Date <- strptime(DID_fires$YYYYMMDD, "%Y%m%d")

# Bounding box time series plots
plot(as.Date(DID_fires$Date), DID_fires$N, type="l", xlim=as.Date(c("2000-01-01","2015-01-01")), xlab="Date", ylab="Number of fires")
plot(as.Date(DID_fires$Date), DID_fires$FRP, type="l", xlim=as.Date(c("2000-01-01","2015-01-01")), xlab="Date", ylab="Mean fire intensities")

# Aggregate observations by AfSIS GID's
GID_fires <- ddply(bbx_fires, .(GID), summarize,
                   x = mean(x),
                   y = mean(y),
                   N = length(GID),
                   maxD = max(YYYYMMDD))
GID_fires$maxD <- strptime(GID_fires$maxD, "%Y%m%d") ## date of last fire
GID_fires$DSLF <- as.numeric(difftime(max(GID_fires$maxD), GID_fires$maxD, units="days")) ## number of days since last fire

# Extract gridded GeoSurvey predictions at fire locations
coordinates(GID_fires) = ~x+y
projection(GID_fires) <- projection(grids)
xgrid <- extract(grids, GID_fires)
ROI_fires <- cbind(GID_fires@coords, GID_fires@data, xgrid)
ROI_fires <- na.omit(ROI_fires)

# ROI ecdf plots of main variables
plot(ecdf(ROI_fires$DSLF), main="", xlab="No. days since last fire", ylab="Cum. proportion of observations", xlim=c(0, 5000), verticals=T, lty=1, lwd=1, col="black", do.points=F)
plot(ecdf(ROI_fires$N), main="", xlab="No. of fires per GID (2000-2015)", ylab="Cum. proportion of observations", xlim=c(0, 6), verticals=T, lty=1, lwd=1, col="black", do.points=F)

# Export fire locations
write.csv(ROI_fires[1:6], "./MCD14ML/TZ_fire_locs.csv", row.names=F)

#+ Exploratory models -----------------------------------------------------
coordinates(ROI_fires) = ~x+y
projection(ROI_fires) = projection(grids)

# GLM's: Days since last fire event (DSLF)
DSLF.glm <- glm(DSLF+1~CRP_ens*RSP_ens*WCP_ens, family=poisson(link="log"), data=ROI_fires)
summary(DSLF.glm)
DSLF.pred <- predict(DSLF.glm, ROI_fires) ## GID-level predictions
quantile(ROI_fires$DSLF, probs=c(0.05, 0.5, 0.95))
quantile(exp(DSLF.pred), probs=c(0.05, 0.5, 0.95))

# DSLF.glm spatial residuals
DSLF.var <- variogram(residuals(DSLF.glm) ~ 1, cutoff = 10000, width = 1000, ROI_fires)
DSLF.vgm <- vgm(model = "Sph", nugget = 700, range = 5000, psill = 900)
DSLF.fit <- fit.variogram(DSLF.var, model = DSLF.vgm)
plot(DSLF.var, DSLF.fit, pc = "+", cex = 2)

# Number of fire events on record (N)
N.glm <- glm(N~CRP_ens*RSP_ens*WCP_ens, family=poisson(link="log"), data=ROI_fires)
summary(N.glm)
N.pred <- predict(N.glm, ROI_fires) ## GID-level predictions
quantile(ROI_fires$N, probs=c(0.05, 0.5, 0.95))
quantile(exp(N.pred), probs=c(0.05, 0.5, 0.95))

# N.glm spatial residuals
N.var <- variogram(residuals(N.glm) ~ 1, cutoff = 10000, width = 1000, ROI_fires)
N.vgm <- vgm(model = "Sph", nugget = 0.4, range = 5000, psill = 0.5)
N.fit <- fit.variogram(N.var, model = N.vgm)
plot(N.var, N.fit, pc = "+", cex = 2)
