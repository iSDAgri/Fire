library(maptools)

# Create a "Data" folder in your current working directory
# dir.create("MCD14ML", showWarnings=F)

filelist <- list.files(path="./MCD14ML", pattern="asc", full.names=T)
bound <- readShapeSpatial("~/MCD14ML/AF_bound.shp")

MCD14ML_data <- NULL

for(file in filelist){
    data <- read.table(file, header=TRUE)
    data <- as.data.frame(data)
    data_loc <- cbind(lon=data$lon, lat=data$lat)
    data_loc <- as.data.frame(data_loc)
    coordinates(data_loc) = ~lon+lat
    temp <- over(data_loc, geometry(bound))
    MCD14ML_data <- rbind(MCD14ML_data, data[!is.na(temp), ])
    print(paste("extracting", file))
}

write.table(MCD14ML_data, "AF_MCD14ML.txt", quote=FALSE, row.names=FALSE)

