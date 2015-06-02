#' Prediction of recent ROI-level MODIS-MCD14ML fire domains 
#' Example ROI: Tanzania
#' Pre-processing with https://github.com/mgwalsh/Fire/blob/master/TZ_fires.R 
#' Source data courtesy UMD (ftp://fuoco.geog.umd.edu/modis/C5/mcd14ml, login=fire, pwd=burnt)
#' M. Walsh, June 2015
 
# Required packages
# install.packages(c("downloader","rgdal","raster","caret")), dependencies=TRUE)
require(downloader)
require(rgdal)
require(raster)
require(caret)

#+ Data downloads ---------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("TZ_MCD14ML", showWarnings=F)
dat_dir <- "./TZ_MCD14ML"

# download Tanzania MCD14ML data
download("https://www.dropbox.com/s/lz4lt31s0h94i9j/TZ_fire_locs.csv.zip?dl=0", "./TZ_MCG14ML/TZ_fire_locs.csv.zip", mode="wb")
unzip("./TZ_MCD14ML/TZ_fire_locs.csv.zip", exdir="./TZ_MCD14ML", overwrite=T)
fires <- read.table(paste(dat_dir, "/TZ_fire_locs.csv", sep=","), header=T)

# download TZ covariate grids and stack in raster
download("https://www.dropbox.com/s/otiqe78s0kf1z1s/TZ_grids.zip?dl=0", "./TZ_MCD14ML/TZ_grids.zip", mode="wb")
unzip("./TZ_MCD14ML/TZ_grids.zip", exdir="./TZ_MCD14ML", overwrite=T)
glist <- list.files(path="./TZ_MCD14ML", pattern="tif", full.names=T)
grids <- stack(glist)

