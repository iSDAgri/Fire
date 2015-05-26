#' Analyses of ROI-level MODIS-MCD14ML fire data 
#' Example ROI: Tanzania 
#' Source data courtesy UMD (ftp://fuoco.geog.umd.edu/modis/C5/mcd14ml, login=fire, pwd=burnt)
#' M. Walsh, May 2015

# Required packages
# install.packages(c("downloader","rgdal","maptools","raster")), dependencies=TRUE)
require(downloader)
require(rgdal)
require(maptools)
require(plyr)
require(raster)

#+ Data downloads ---------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("MCD14ML", showWarnings=F)
dat_dir <- "./MCD14ML"

# download Africa-wide MCD14ML data (large download: ~490 Mb)
download("https://www.dropbox.com/s/oi0nvimjawmz7xf/AF_MCD14ML.txt.zip?dl=0", "./MCG14ML/AF_MCD14ML.txt.zip", mode="wb")
unzip("./MCD14ML/AF_MCD14ML.txt.zip", exdir="./MCD14ML", overwrite=T)
fires <- read.table(paste(dat_dir, "/AF_MCD14ML.txt", sep=""), header=T)

# download ROI shape file
download("https://www.dropbox.com/s/bb5277wvuce0b2h/TZ_bound.zip?dl=0", "./MCD14ML/TZ_bound.zip", mode="wb")
unzip("./MCD14ML/TZ_bound.zip", exdir="./MCD14ML", overwrite=T)
bound <- readShapeSpatial("./MCD14ML/TZ_bound.shp")

# download GeoSurvey predictions and stack in raster
download("https://www.dropbox.com/s/0ucx0oqx8c4lpej/TZ_grids.zip?dl=0", "./MCD14ML/TZ_grids.zip", mode="wb")
unzip("./MCD14ML/TZ_grids.zip", exdir="./MCD14ML", overwrite=T)
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

# Generate AfSIS grid cell ID's (GID)
res.pixel <- 1000
xgid <- ceiling(abs(bbx_fires$x)/res.pixel)
ygid <- ceiling(abs(bbx_fires$y)/res.pixel)
gidx <- ifelse(bbx_fires$x<0, paste("W", xgid, sep=""), paste("E", xgid, sep=""))
gidy <- ifelse(bbx_fires$y<0, paste("S", ygid, sep=""), paste("N", ygid, sep=""))
GID <- paste(gidx, gidy, sep="-")
bbx_fires <- cbind(bbx_fires, GID)

# Aggregate observations by GID's
GID_fires <- ddply(bbx_fires, c("GID"), summarise,
                   X = mean(x),
                   Y = mean(y),
                   N = length(GID),
                   minD = min(YYYYMMDD),
                   maxD = max(YYYYMMDD),                                     
                   FRP = mean(FRP))

# Aggregate observations by date
DID_fires <- ddply(bbx_fires, c("YYYYMMDD"), summarise,
                   D = mean(YYYYMMDD),
                   N = length(YYYYMMDD),
                   FRP = mean(FRP))

