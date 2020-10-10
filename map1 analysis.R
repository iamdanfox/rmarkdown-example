
# https://rspatial.org/spatial/5-files.html - tutorials on spatial data in R
# install programmes
install.packages("dplyr")
install.packages("tidyverse") # dplyr for data manipulation
# ggplot2 for data visualization
# tidyr for data reshaping (which is a sub-task of data manipulation)
# lubridate to deal with dates
# forcats to deal with R's factor variable type
# stringr to deal with string data
install.packages("readxl")
install.packages("readr")
install.packages("testthat")
install.packages("sp")
install.packages("sf")
install.packages("raster")
install.packages("rgdal")
install.packages("maptools")
install.packages("rasterVis")
install.packages("tmap")
install.packages("tmaptools")
install.packages("shiny")
install.packages("shinyjs")
install.packages("datasets")
install.packages("forcats")
install.packages("ggplot2")
install.packages("grDevices")
install.packages("graphics")
install.packages("lattice")
install.packages("latticeExtra")
install.packages("methods")
install.packages("RColorBrewer")
install.packages("stats")
install.packages("stringr")
install.packages("testthat")
install.packages("tibble")
install.packages("tidyr")
install.packages("tidyverse")
install.packages("utils")
# install.packages("purr") # Error: package 'purr' is not available
# install.packages("GIS")
install.packages("rgeos")
# install.packages("velox")
install.packages("sf")
install.packages("future")
install.packages("future.apply")
install.packages("tictoc")

library(dplyr)
library(tidyverse)
library(readxl)
library(readr)
library(testthat)
library(sp)
library(sf)
library(raster)
library(rgdal)
library(maptools)
library(rasterVis)
library(tmap)
library(tmaptools)
library(shiny)
library(shinyjs)
library(datasets)
library(forcats)
library(ggplot2)
library(grDevices)
library(graphics)
library(lattice)
library(latticeExtra)
library(methods)
library(RColorBrewer)
library(stats)
library(stringr)
library(testthat)
library(tibble)
library(tidyr)
library(tidyverse)
library(utils)
library(rgeos)
library(velox)
library(sf)
library(future)
library(future.apply)
library(tictoc)


# test practice
#####################################################
# #testexample - making a map using point, line and polygon data (points, lines and polygons are the 3 classes of vector data - a class is a geometry only and needs to have attributes added)
# ##create some coordinate data and assign a CRS
# longitude <- c(-116.7, -120.4, -116.7, -113.5, -115.5, -120.8, -119.5, -113.7, -113.7, -110.7)
# latitude <- c(45.3, 42.6, 38.9, 42.1, 35.7, 38.9, 36.2, 39, 41.6, 36.9)
# lonlat <- cbind(longitude, latitude)
# pts <- SpatialPoints(lonlat)
# # class(pts)
# # showDefault(pts)
# coordinate_reference<-CRS("+proj=longlat +datum=WGS84")#this defines a CRS  # https://proj.org/usage/projections.html - details of syntax to change coordinate systems in R (using regal)
# pts<- SpatialPoints(lonlat, proj4string=coordinate_reference)# this assignes the CRS to the dataset
# ##create a dataframe with the coordinates and some attribute values
# precip_value<-runif(nrow(lonlat),min=0,max=100)
# df<-data.frame(ID=1:nrow(lonlat),precip=precip_value) #first create a dataframe with same number of rows as coordinate data
# pts_df<-SpatialPointsDataFrame(pts,df)#then add coordinate data
# showDefault(pts_df)#show what's inside
# ##use line and polygon data
# lns<-spLines(lonlat,crs = coordinate_reference)
# pols<-spPolygons(lonlat, crs= coordinate_reference)
# ##plot data
# plot(pols, axes=TRUE, las=1)
# plot(pols, border='blue', col='yellow', lwd=3, add=TRUE)
# points(pts, col='red', pch=20, cex=3)
#
# #testexample - making a map using raster data
# r<-raster(ncol=10,nrow=10, xmx=80,xmn=-80,ymx=50,ymn=-50)
# values(r)<-runif(ncell(r))
# values(r)<-1:ncell(r)
# plot(r)
# lon <- c(-68.8, 27.2, -35.9, -34.9, -46.2, -35.4, -15.7)
# lat <- c(41.3, 42.9, 42.4, 39.8, 37.6, 38.3, 37.6)
# lonlat <- cbind(lon, lat)
# pols <- spPolygons(lonlat, crs='+proj=longlat +datum=WGS84')
# points(lonlat, col='red', pch=20, cex=3)
# plot(pols, border='blue', lwd=2, add=TRUE)
#
# #n.b. Rasterstack and rasterbrick objects both contain multiple spatial layers, but all layers from a rasterstack can come from multiple files, whereas all layers from a rasterbrick must come from the same file.
# #e.g.
# r2 <- r * r
# r3  <- sqrt(r)
# s <- stack(r, r2, r3)
# b<-brick(s)
#
# #test example - raster algebra
# r1<-raster(ncol=10,nrow=10, xmx=20, xmn=-20, ymx=20,ymn=-20)
# values(r1)<-1:ncell(r1)
# r2<-raster(ncol=10,nrow=10, xmx=20, xmn=-20, ymx=20,ymn=-20)
# values(r2)<-1:ncell(r2)
# s<-stack(r1,r2) ##these two rasters cannot be stacked if they do not have same resolution and origin value (i.e. value closest to raster corners)
#
#
# a<-r1+10 # this changes every value of the raster
# b<-r1*10
# r1[]<- runif(ncell(r1))
# r1<-round(r1)
# s[r1]<-0.5
#
# c<-mean(r1,r2,10) #this gives mean of these values for every cell (not a summary)
# cellStats(r1,sum) # this gives summary statistics for cell as a whole
#
#
# #########################################################
# Question 8. Where are smallholder farmers?

# import data
setwd("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/Map 1")
# readTIFF("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/Worldpop population data/Africa_1km_Population/AFR_PPP_2020_adj_v2.tif", native = FALSE, all = FALSE, convert = FALSE, info = FALSE, indexed = FALSE, as.is = FALSE)
library(raster)
library(sf)
gadm_SSA_dissolved <- st_read("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/Map 1/gadm_36_0 SSA dissolved.shp")
# farming_systems<-shapefile("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/Dixons farming systems/ssa_fs_final.shp")
## n.b. this imports .shp file as a Spatial*DataFrame (SpatialPolygonsDataFrame in this case)
# africa_pop_raster<-raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/Worldpop population data/Africa_1km_Population/AFR_PPP_2020_adj_v2.tif")
# field_sizes_raster<-raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/IIASA global field sizes/Global Field Sizes 2/dominant_field_size_categories.tif")
hi_binary <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/Map 1/hi_av_2030_rcp85_clipped_to_2020_pop_binary.tif")
heat_crop_binary <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/Map 1//heat_crop_generic_2030_Africa_monthly_clipped_to_SPAM_crops_binary.tif")
rasterised_SSA_weird_values <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/gadm36_levels_shp_SSA/rasterised_weird_values.tif")

thi <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/Julian's data/thi/thi_av_2030_rcp85.asc")


Bf_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Bf_2010_Aw.tif")
Ch_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Ch_2010_Aw.tif")
Ct_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Ct_2010_Aw.tif")
Dk_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Dk_2010_Aw.tif")
Gt_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Gt_2010_Aw.tif")
Ho_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Ho_2010_Aw.tif")
Pg_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Pg_2010_Aw.tif")
Sh_2010_Aw <- raster("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/All raw data/FAO gridded livestock /FAO (2018) -- Global Livestock Distribution in 2010.geotiff/6_Sh_2010_Aw.tif")


# set palette
palette_1 <- c("#f2f268", "#add45b", "#52b35d", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff")
palette_2 <- c("#f2f268", "#add45b", "#52b35d", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff")
palette_3 <- brewer.pal(n = 8, name = "YlGnBu")
palette_4 <- brewer.pal(n = 8, name = "RdYlGn")
palette_5_sat_all <- c("#e30613", "#f39200", "#feeb18", "#009640", "#a90f09", "#b47014", "#c1b100", "#007330", "#6c0d00", "#724700", "#78701e", "#004b1d", "#1d1d1b")
palette_5_sat_selected_old <- c("#1d1d1b", "#004b1d", "#78701e", "#007330", "#c1b100", "#b47014", "#009640", "#feeb18", "#f39200", "#e30613")
palette_5_sat_selected <- c("#9d9c9c", "#d6e9d8", "#fff9c7", "#9bcea4", "#fff381", "#f8b77c", "#009640", "#feeb18", "#f39200", "#e30613")














theme1 <- custom.theme(
  symbol = palette_1,
  fill = c("#f2f268", "#add45b", "#52b35d", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff"),
  region = c("#f2f268", "#add45b", "#52b35d", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff"),
  bg = "#fcf6e8", fg = "#000000"
)
theme2 <- custom.theme(
  symbol = palette_2,
  fill = c("#f2f268", "#add45b", "#52b35d", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff"),
  region = c("#f2f268", "#add45b", "#52b35d", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff", "#00fbff"),
  bg = "#fcf6e8", fg = "#000000"
)
theme3 <- custom.theme(
  symbol = palette_4,
  fill = brewer.pal(n = 6, name = "YlGnBu"),
  region = brewer.pal(n = 6, name = "YlGnBu"),
  bg = "#fcf6e8", fg = "#000000"
)



# shapefile(s,"test.shp", overwrite=TRUE)
# writeRaster(africa_pop_raster,"/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/delete.tif")

# set coordinate reference system
## find some coordinate reference systems from the EPSG database which match the location
library(rgdal)
# espg<-make_EPSG()
# espg$note
# German_CRSs<-grep("DRC",espg$note, ignore.case = TRUE)
# espg[German_CRSs[1:2],] #n.b specific rows and columns can be selected by name, e.g. x[c("row1","row2),]
# ##check the coordinate reference system of the dataset
crs(hi_binary)
# crs(thi_binary)
crs(heat_crop_binary)
# change the coordinate reference system if required (fine for vector, not fine for raster unless absolutely necessary because values of the new cells must be based on values of the old cells, so loss of precision). Same process used for raster if required.
# library(rgdal)
# new_crs<-CRS("+proj=robin + datum=WGS84")
# countries_raster_robinson_crs<-spTransform(countries_shp,new_crs)

# manipulate vector and raster data
# geom(countries_shp)
# head(countries_shp)
# countries_shp$new <- sample(x=letters, size=length(countries_shp),replace=TRUE)
# countries_shp$new<-NULL

# display vector and raster data
# plot(countries_shp)

# display rasters

## plot and get information
# plot(thi_binary)
# res(thi_binary)# this is the size of each cell (x by y) - higher means bigger cells
# dim(thi_binary)# this is the number of cells on x and y - it cannot be changed without changing the values in the cells
# xmax(thi_binary)# this is the max x extent (bounding box) of the raster  - it can be changed without changing the values inside the cells
# inMemory(thi_binary)# tells us whether cell values are stored in R's memory
# plot(thi_binary)
# res(thi_binary)
# dim(thi_binary)
# xmax(thi_binary)
# inMemory(thi_binary)
## this is what happens if you change number of cells
# xmax(field_sizes_raster)<-0
# ncol(field_sizes_raster)<-6
# hasValues(field_sizes_raster)


# mask field sizes with population map
# s<-stack(field_sizes_raster, africa_pop_raster) # can't do this because different origin and resolution

# mask smallest field size with population
# thi_binary_cropped<-crop(x=thi_binary,y=gadm_SSA_dissolved) n.b. doesnt work - maybe because crop wrong function for vecotr and raster?

# field_size_africa_pop<-merge(field_size_cropped, africa_pop_raster)# not possible due to different resolution and origins

# res(thi_binary) #resolution is different
# res(gadm_SSA_dissolved)
# origin(thi_binary)# origin is different
# origin(gadm_SSA_dissolved)
# extent(thi_binary)
# extent(gadm_SSA_dissolved)
# getValuesBlock(gadm_SSA_dissolved,4000,4000)

## use the dis(aggregate) or resample functions to make resolution same (this also makes origin the same)
# field_size_cropped_resampled<-resample(field_size_cropped, africa_pop_raster, "ngb")

# mask population by smallest field size
# is_u0.64ha<- calc(field_size_cropped_resampled, fun=function(field_size_cropped_resampled){ field_size_cropped_resampled[field_size_cropped_resampled != 3506] <- NA; return(field_size_cropped_resampled)} )
# as.matrix(is_u0.64ha)
# pop_u0.64ha<-mask(x=africa_pop_raster,mask=is_u0.64ha, updatevalue=NA)


# reclassify values away from 0 and 1 to 1 and 10 and change NAs to 0
rcl <- matrix(c(0, 1, 1, 10), nrow = 2, ncol = 2)


hi_binary_reclassified <- reclassify(hi_binary, rcl)
heat_crop_binary_reclassified <- reclassify(heat_crop_binary, rcl)

rcl <- matrix(c(0, 1, 1, 1000000000000, 1, 10), nrow = 2, ncol = 3)
Bf_2010_Aw <- reclassify(Bf_2010_Aw, rcl)
Ch_2010_Aw <- reclassify(Ch_2010_Aw, rcl)
Ct_2010_Aw <- reclassify(Ct_2010_Aw, rcl)
Dk_2010_Aw <- reclassify(Dk_2010_Aw, rcl)
Gt_2010_Aw <- reclassify(Gt_2010_Aw, rcl)
Ho_2010_Aw <- reclassify(Ho_2010_Aw, rcl)
Pg_2010_Aw <- reclassify(Pg_2010_Aw, rcl)
Sh_2010_Aw <- reclassify(Sh_2010_Aw, rcl)


Bf_2010_Aw[is.na(Bf_2010_Aw)] <- 0
Ch_2010_Aw[is.na(Ch_2010_Aw)] <- 0
Ct_2010_Aw[is.na(Ct_2010_Aw)] <- 0
Dk_2010_Aw[is.na(Dk_2010_Aw)] <- 0
Gt_2010_Aw[is.na(Gt_2010_Aw)] <- 0
Ho_2010_Aw[is.na(Ho_2010_Aw)] <- 0
Pg_2010_Aw[is.na(Pg_2010_Aw)] <- 0
Sh_2010_Aw[is.na(Sh_2010_Aw)] <- 0


thi_binary <- Bf_2010_Aw + Ch_2010_Aw + Ct_2010_Aw + Dk_2010_Aw + Gt_2010_Aw + Ho_2010_Aw + Pg_2010_Aw + Sh_2010_Aw

rcl <- matrix(c(0, 8, 8, 81, 1, 10), nrow = 2, ncol = 3)
thi_binary_reclassified <- reclassify(thi_binary, rcl)
thi <-

  rasterised_SSA_reclass_matrix <- matrix(c(0, 10000000000000, 1), nrow = 1, ncol = 3)
rasterised_SSA <- reclassify(rasterised_SSA_weird_values, rasterised_SSA_reclass_matrix)


thi_binary_reclassified <- crop(thi_binary_reclassified, rasterised_SSA)
thi_binary_reclassified <- thi_binary_reclassified * rasterised_SSA

thi_binary_reclassified[is.na(thi_binary_reclassified)] <- 0

hi_binary_reclassified[is.na(hi_binary_reclassified)] <- 0
heat_crop_binary_reclassified[is.na(heat_crop_binary_reclassified)] <- 0

hi_binary_reclassified <- hi_binary_reclassified * rasterised_SSA
heat_crop_binary_reclassified <- heat_crop_binary_reclassified * rasterised_SSA






# NAvalue(thi_binary_reclassified)
# NAvalue(thi_binary_reclassified) <- 0
# NAvalue(hi_binary_reclassified)
# NAvalue(hi_binary_reclassified) <- 0
# NAvalue(heat_crop_binary_reclassified)
# NAvalue(heat_crop_binary_reclassified) <- 0




# Combine values into new map
combined <- thi_binary_reclassified + hi_binary_reclassified + heat_crop_binary_reclassified
# histogram(combined)

reclass_matrix <- matrix(c(3, 12, 21, 30, 2, 1, 0, 11, 20, 10, 30, 31, 32, 33, 20, 10, 0, 21, 22, 11), nrow = 10, ncol = 2)

combined_reclassified <- reclassify(combined, reclass_matrix)





##### stoppped here Firday 2nd october

### started here 5th october


combined <- mask(combined_reclassified, gadm_SSA_dissolved)


combined_plot <- tm_shape(combined) +
  tm_raster(
    col = NA,
    alpha = NA,
    palette = palette_5_sat_selected,
    n = 10,
    style = "cat",
    style.args = list(),
    as.count = NA,
    breaks = NULL,
    interval.closure = "left",
    labels = c(
      "0 datasets",
      "no heat stress (score 0), 1 dataset (incomplete data)",
      "low heat stress (score 1), 1 dataset (incomplete data)",
      "no heat stress (score 0), 2 datasets (incomplete data)",
      "low heat stress (score 1), 2 datasets (incomplete data)",
      "medium heat stress (score 2), 2 datasets (incomplete data)",
      "no heat stress (score 0), 3 datasets (complete data)",
      "low heat stress (score 1), 3 datasets (complete data)",
      "medium heat stress (score 2), 3 datasets (complete data)",
      "high heat stress (score 3), 3 datasets (complete data)"
    ),
    drop.levels = FALSE,
    midpoint = NULL,
    stretch.palette = TRUE,
    contrast = NA,
    saturation = 1,
    interpolate = NA,
    colorNA = NULL,
    textNA = "Missing",
    showNA = TRUE,
    colorNULL = NULL,
    title = "",
    legend.show = TRUE,
    legend.format = list(),
    legend.is.portrait = TRUE,
    legend.reverse = FALSE,
    legend.hist = FALSE,
    legend.hist.title = NA,
    legend.z = NA,
    legend.hist.z = NA,
    zindex = NA,
    group = NULL,
    auto.palette.mapping = NULL,
    max.categories = NULL,
    max.value = 255
  ) +
  tm_layout(
    legend.frame = FALSE,
    legend.outside = FALSE, legend.position = c("left", "centre"),
    title = "Heat stress in livestock, humans and crops across SSA",
    title.position = c("centre", "top"),
    title.size = 0.5,
    legend.bg.color = "#ffffff",
    legend.text.size = 0.4,
    legend.width = 0.4
  )





writeRaster(combined, "/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/Map 1/combined.tif", overwrite = TRUE)

combined_plot
dev.off()
save.image("/Users/katherinehickson/Library/Mobile Documents/com~apple~CloudDocs/Adaptation atlas/Investor report/Map 1/combined_plot.pdf")















## remove outlier values
hist(pop_u0.64ha,
  main = "Distribution of population values",
  xlab = "population in ~1km2", ylab = "Frequency",
  col = "springgreen"
)
maxValue(pop_u0.64ha)
unique(pop_u0.64ha)
vls <- rasterToPoints(pop_u0.64ha, spatial = FALSE)
vls <- as.data.frame(vls)
head(vls)
smm <- summary(vls$AFR_PPP_2020_adj_v2)
pop_u0.64ha_no_outliers <- pop_u0.64ha
pop_u0.64ha_no_outliers[which(pop_u0.64ha_no_outliers[] > 10000)] <- NA


# make map

# pop_field_size_map<-levelplot(pop_u0.64ha, margin=FALSE, par.settings = theme3,
#           main="Smallholder population per ~1km2 where dominent \nfield size is <0.64Ha, Africa", maxpixels=2e6)

#

countries_shp_cropped <- crop(x = countries_shp, y = pop_u0.64ha)


pop_field_size_map <- tm_shape(pop_u0.64ha) +
  tm_layout(
    legend.frame = TRUE,
    legend.outside = FALSE, legend.position = c("left", "top"),
    main.title = "Smallholder population (<0.64ha land plots)",
    main.title.size = 1,
    main.title.position = "centre",
    legend.bg.color = "#ffffff",
    legend.text.size = 0.7
  ) +
  tm_raster(
    title = "",
    col = "std", palette = palette_4,
    legend.hist = FALSE, legend.is.portrait = TRUE
  ) +
  tmap_style("white") +
  tm_shape(countries_shp_cropped) +
  tm_borders("black", lwd = .5)

# sum values
total_u0.64ha_pop <- cellStats(pop_u0.64ha, sum)
mean_u0.64ha_pop <- cellStats(pop_u0.64ha, mean)
total_africa_pop <- cellStats(africa_pop_raster, sum)
total_pop_u0.64ha_no_outliers <- cellStats(pop_u0.64ha_no_outliers, sum)

# find number of smallholder farmers in each farming system
# farming_systems()
# showDefault(farming_systems)#show what's inside
# unique(farming_systems$DESC)
# geom(farming_systems)
# farming_systems_raster<-rasterize(farming_systems, africa_pop_raster)

# extent(farming_systems)
# extent(farming_systems_raster)
# extent(africa_pop_raster)

#
# # farming_systems_dataframe<-as.data.frame(farming_systems_raster) - vector memory exhausted
# # merge(farming_system_pop,farming_systems_raster)
#
# ##Try to make 14 vector layers and extract for each . Simple.
# farming_systems_1<-farming_systems[which(farming_systems$DESC=="1. Irrigated"),]
# farming_systems_2<-farming_systems[which(farming_systems$DESC=="2. Tree crop" ),]
# farming_systems_3<-farming_systems[which(farming_systems$DESC=="3. Forest based"),]
# farming_systems_4<-farming_systems[which(farming_systems$DESC=="4. Rice-tree crop"),]
# farming_systems_5<-farming_systems[which(farming_systems$DESC=="5. Highland perennial"),]
# farming_systems_6<-farming_systems[which(farming_systems$DESC=="6. Highland temperate mixed"),]
# farming_systems_7<-farming_systems[which(farming_systems$DESC=="7. Root crop"),]
# farming_systems_8<-farming_systems[which(farming_systems$DESC=="8. Cereal-root crop mixed"),]
# farming_systems_9<-farming_systems[which(farming_systems$DESC=="9. Maize mixed"),]
# farming_systems_10<-farming_systems[which(farming_systems$DESC=="10. Large commercial & smallholder"),]
# farming_systems_11<-farming_systems[which(farming_systems$DESC=="11. Agro-pastoral millet/sorghum"),]
# farming_systems_12<-farming_systems[which(farming_systems$DESC=="12. Pastoral"),]
# farming_systems_13<-farming_systems[which(farming_systems$DESC=="13. Sparse (arid)"),]
# farming_systems_14<-farming_systems[which(farming_systems$DESC=="14. Coastal artisanal fishing"),]
# farming_systems_list<-c(farming_systems_1, farming_systems_2, farming_systems_3, farming_systems_4
#                            , farming_systems_5, farming_systems_6, farming_systems_7, farming_systems_8
#                            , farming_systems_9, farming_systems_10, farming_systems_11, farming_systems_12
#                            , farming_systems_13, farming_systems_14)
#
#
# ##method 1 - by farming system, extract the sum of each polygon
# library(raster)
# m1_pop_farming_systems_1_sum_by_feature<-raster::extract(africa_pop_raster,farming_systems_1, fun=sum)
# m1_pop_farming_systems_1_sum_by_feature_df<-as.data.frame(m1_pop_farming_systems_1_sum_by_feature)
# m1_pop_farming_systems_1_total_1<-sum(m1_pop_farming_systems_1_sum_by_feature_df$V1, na.rm=TRUE)
# m1_pop_farming_systems_1_total_2<-sum(m1_pop_farming_systems_1_sum_by_feature, na.rm=TRUE)
#
# m1_pop_farming_systems_2_sum_by_feature<-raster::extract(africa_pop_raster,farming_systems_2, fun=sum)
# m1_pop_farming_systems_2_sum_by_feature_df<-as.data.frame(m1_pop_farming_systems_2_sum_by_feature)
# m1_pop_farming_systems_2_total_1<-sum(m1_pop_farming_systems_2_sum_by_feature_df$V1, na.rm=TRUE)
# m1_pop_farming_systems_2_total_2<-sum(m1_pop_farming_systems_2_sum_by_feature, na.rm=TRUE)
#
# for (i in 1:length(farming_systems_list)){
#   extract<-extract(africa_pop_raster,farming_systems_1[i], fun=sum) #xx class(farming_systems_1[i]) is list (even though farming_systems_1[i] shows as SpatialPolygonsDataFrame). How get class to show as spaatial polygons dataframe so that extract can be performed in loop?
#   df<-as.data.frame(extract)
#   total<-sum(df,na.rm=TRUE)
#   print(paste0("farming_systems_",i,"sums to",total))
# }
#
# ##method 2 - by farming system, extract all 'values of the cells of a Raster* object that are covered by a polygon' and then sum
# m2_pop_farming_systems_1_extract<-extract(africa_pop_raster,farming_systems_1)#
# m2_pop_farming_systems_1_sum_by_feature<-lapply(m2_pop_farming_systems_1_extract, function(x){if (!is.null(x)) sum(x, na.rm = TRUE)})
# m2_pop_farming_systems_1_sum_by_feature_unlisted<-unlist(m2_pop_farming_systems_1_sum_by_feature)
# m2_pop_farming_systems_1_total_1<-sum(m2_pop_farming_systems_1_sum_by_feature_unlisted)
#
# m2_pop_farming_systems_2_extract<-extract(africa_pop_raster,farming_systems_2)#
# m2_pop_farming_systems_2_sum_by_feature<-lapply(m2_pop_farming_systems_2_extract, function(x){if (!is.null(x)) sum(x, na.rm = TRUE)})
# m2_pop_farming_systems_2_sum_by_feature_unlisted<-unlist(m2_pop_farming_systems_2_sum_by_feature)
# m2_pop_farming_systems_2_total_1<-sum(m2_pop_farming_systems_2_sum_by_feature_unlisted)
#
# for (i in 1:length(farming_systems_list)){
#   extract<-extract(africa_pop_raster,farming_systems_list[i])
#   lapply_to_result<-lapply(extract, function(x){if (!is.null(x)) sum(x, na.rm = TRUE)})
#   unlist<-unlist(lapply_to_result)
#   total<-sum(unlist)
#   print(paste0("farming_systems_",i,"sums to",total))
# }
#
# ## xx why different answers for methods 1 and 2???
#
#
#


## method 3 - extract everything from all farming systems at once
# zonal(africa_pop_raster, farming_systems_raster, fun=sum) - RasterLayers cannot be processed in memory.
library(raster)
farming_system_pop <- raster::extract(africa_pop_raster, farming_systems, fun = sum, na.rm = TRUE, df = TRUE)
farming_systems$africa_pop <- farming_system_pop
farming_systems_df <- data.frame(farming_systems)
farming_system_u0.64ha_pop <- raster::extract(pop_u0.64ha, farming_systems, fun = sum, na.rm = TRUE, df = TRUE)
farming_systems$pop_u0.64ha <- farming_system_u0.64ha_pop
farming_systems_df <- data.frame(farming_systems)
farming_systems_smallholder_sum <- sum(farming_systems_df$pop_u0.64ha, na.rm = TRUE)
farming_systems_population_sum <- sum(farming_systems_df$africa_pop, na.rm = TRUE)


# plot(africa_pop_raster)
# plot(farming_systems, bg="transparent", add=TRUE)

farming_systems_pop_summary <- summarise(group_by(data.frame(farming_systems), DESC), sum_africa_pop = sum(africa_pop), sum_pop_u0.64ha = sum(pop_u0.64ha))
farming_systems_pop_summary$proportion_u0.64ha <- farming_systems_pop_summary$sum_pop_u0.64ha / farming_systems_pop_summary$sum_africa_pop


# find number of smallholder farmers in each country

extent(countries_shp)
extent(africa_pop_raster)
countries_shp_df <- data.frame(countries_shp)
# ##work out whether geom ID corresponds to ID_0 in dataframe
# test=countries_shp[which(countries_shp_df$ID_0==1),]
# g<-geom(test)
# nrow(g)
# nrow(countries_shp_df[which(countries_shp_df$ID_0==1),])
# ###'Object' in geom corresponds to UID in dataframe

# plot(africa_pop_raster)
# plot(countries_shp, bg="transparent", add=TRUE)

## dissolve country shapes
#
# dissolve <- function(SPDF,attribute){
#   rownames(SPDF@data) <- sapply(slot(SPDF, "polygons"), function(x) slot(x, "ID"))
#   Temp <- gUnaryUnion(SPDF,attribute)
#   IDlist <- data.frame(ID=sapply(slot(Temp, "polygons"), function(x) slot(x, "ID")))
#   rownames(IDlist)  <- IDlist$ID
#   SpatialPolygonsDataFrame(Temp,IDlist)
# }
# countries_shp_dissolve<-dissolve(SPDF = countries_shp,attribute = countries_shp$NAME_0)

## extract raster from country shapes
countries_shp_pop <- raster::extract(africa_pop_raster, countries_shp, fun = sum, na.rm = TRUE, df = TRUE)
countries_shp_pop_df <- data.frame(countries_shp)
countries_shp$africa_pop <- countries_shp_pop
countries_shp_df <- data.frame(countries_shp)
countries_shp_u0.64ha_pop <- raster::extract(pop_u0.64ha, countries_shp, fun = sum, na.rm = TRUE, df = TRUE)
countries_shp$pop_u0.64ha <- countries_shp_u0.64ha_pop
countries_shp_df <- data.frame(countries_shp)
group_by(countries_shp_df, NAME_0)

countries_pop_summary <- summarise(group_by(data.frame(countries_shp), NAME_0),
  sum_pop_u0.64ha = sum(pop_u0.64ha),
  sum_africa_pop = sum(africa_pop)
)
countries_pop_summary$proportion_u0.64ha <- countries_pop_summary$sum_pop_u0.64ha / countries_pop_summary$sum_africa_pop

smallholder_population_country_sum <- sum(countries_pop_summary$sum_pop_u0.64ha, na.rm = TRUE)
population_country_sum <- sum(countries_pop_summary$sum_africa_pop, na.rm = TRUE)

# write outputs
writeRaster(pop_u0.64ha, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/pop_u0.64ha.tif", overwrite = TRUE)
write_csv(farming_systems_pop_summary, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/farming_systems_pop_summary.csv")
write_csv(countries_pop_summary, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/countries_pop_summary.csv")
jpeg("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/pop_field_size_map.png")
pop_field_size_map
dev.off()
save.image("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/global_env.RData")




#####################################################################

# Question 1 - Where are smallholder farmers vulnerable?



# source existing data and variables
# plot(africa_pop_raster)
# plot(pop_u0.64ha)
# plot(farming_systems)

# create new data and variables
# Lgpdelt: more than 5% decrease in length of growing period (LGP).
# Lpgflip: LGP flips from more than 120 days to less than 120 days (threshold 1).
# Rcgd: reliable crop growing days (RGCD) flip from more than 90 days to less than 90 days (threshold 2)
# Tmax: maximum temperature flips from <30 deg C to > 30 deg C (threshold 4)
# Tgrow: maximum temperature during the growing season flips from <30 deg C to > 30 deg C (threshold 5).
# Rdaydec: rainfall per rainy day decreases > 10% (threshold 6)
# Rdayinc: rainfall per rainy day increases > 10% (threshold 7)
# CV: coefficient of variability of rainfall currently greater than 21%
# Tmean: mean annual temperature flips from < 8 deg C to >8 deg C (threshold 3)

Lgpdelt <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_lgpdelt2/sta.adf")
Lpgflip <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_lgpflip2/sta.adf")
Rcgd <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_rcgd2/sta.adf")
Tmax <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_tmax2/sta.adf")
Tgrow <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_tgrow2/sta.adf")
Rdaydec <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_rdaydec2/sta.adf")
Rdayinc <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_rdayinc2/sta.adf")
CV <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_cv2/sta.adf")
Tmean <- raster("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Raw data real atlas/An, Polly, Phillip Vulnerability /Vulnerability report vulnerability data phil and an notenberg/domains_documented/dom_tmean2/sta.adf")

hist(Lgpdelt,
  main = "Distribution of vulnerability values",
  xlab = "value", ylab = "Frequency",
  col = "springgreen"
)
maxValue(Lgpdelt)
unique(Lgpdelt)


## resample
## use the dis(aggregate) or resample functions to make resolution same (this also makes origin the same)
Lgpdelt_resampled <- resample(Lgpdelt, africa_pop_raster, "ngb")




# ###############################
#
# # method for making a table and raster with sum u0.64ha pop by vulnerability class
# ## extract dummy data
# uga <- raster::getData('GADM', country = 'UGA', level = 0)
# r1 <- Lgpdelt_resampled %>% raster::crop(., uga) %>% raster::mask(., uga)
# r2 <-  pop_u0.64ha  %>% raster::crop(., uga) %>% raster::mask(., uga)
#
# ## Making a table (summarizing by r1)
# tb1 <- r1 %>% rasterToPoints %>% as.data.frame %>% as_tibble
# tb2 <- r2 %>% rasterToPoints %>% as.data.frame %>% as_tibble
#
# tbl <- merge(tb1, tb2) %>% as_tibble()
# smm <- tbl %>%
#   group_by(sta) %>%
#   dplyr::summarise(value = sum(AFR_PPP_2020_adj_v2, na.rm = TRUE)) %>%
#   ungroup()
#
# ## Creating a raster
# fnl <- r1
# fnl[which(fnl[] == 111)] <- as.numeric(smm[1,2])
# fnl[which(fnl[] == 121)] <- as.numeric(smm[2,2])
# fnl[which(fnl[] == 211)] <- as.numeric(smm[3,2])
# fnl[which(fnl[] == 212)] <- as.numeric(smm[4,2])
# fnl[which(fnl[] == 221)] <- as.numeric(smm[5,2])
# fnl[which(fnl[] == 222)] <- as.numeric(smm[6,2])
# #
# # options(scipen = 999)
# # plot(fnl)
# #
# #
#
# tm_shape(fnl) +
#   tm_layout(legend.frame = TRUE,
#             legend.outside = FALSE, legend.position = c("left","top"),
#             main.title = "Total no. of farmers per vulnerability domain",
#             main.title.size = 1,
#             main.title.position = "centre",
#             legend.bg.color = "#ffffff",
#             legend.text.size =1)+
#   tm_raster(title = "" ,
#             col="sta",palette = palette_4,  n=6,
#             labels=c("HLH, 15,196","HHH, 80,618","LLL, 2,384,112","HLL, 4,134,570","LHL, \n4,648,838","HHL, \n17,789,102"),
#             legend.hist = FALSE, legend.is.portrait = TRUE, style="cat") +
#   tmap_style("natural")
#


# #Vulnerability code (exposure, sensitivity and adaptive capacity) and total number of african smallholders
# ####################################


# making a table (summarising by Lgpdelt)
Lgpdelt_resampled_2 <- raster::rasterToPoints(Lgpdelt_resampled)
Lgpdelt_resampled_3 <- as.data.frame(Lgpdelt_resampled_2)
pop_u0.64ha_2 <- raster::rasterToPoints(pop_u0.64ha)
pop_u0.64ha_3 <- as.data.frame(pop_u0.64ha_2)
africa_pop_raster_points <- raster::rasterToPoints(africa_pop_raster)
africa_pop_raster_points2 <- as.data.frame(africa_pop_raster_points)

Lgpdelt_u0.64ha_pop_table <- merge(Lgpdelt_resampled_3, pop_u0.64ha_3)
Lgpdelt_u0.64ha_pop_table_2 <- group_by(Lgpdelt_u0.64ha_pop_table, sta)
Lgpdelt_u0.64ha_pop_table_3 <- dplyr::summarise(Lgpdelt_u0.64ha_pop_table_2, value = sum(AFR_PPP_2020_adj_v2, na.rm = TRUE)) # xx
Lgpdelt_u0.64ha_pop_table_4 <- ungroup(Lgpdelt_u0.64ha_pop_table_3)
smallholder_population_in_Lgpdelt_sum <- sum(Lgpdelt_u0.64ha_pop_table_4$value)

Lgpdelt_africa_pop_raster_points2_table <- merge(Lgpdelt_resampled_3, africa_pop_raster_points2)
Lgpdelt_africa_pop_raster_points2_table_2 <- group_by(Lgpdelt_africa_pop_raster_points2_table, sta)
Lgpdelt_africa_pop_raster_points2_table_3 <- dplyr::summarise(Lgpdelt_africa_pop_raster_points2_table_2, value = sum(AFR_PPP_2020_adj_v2, na.rm = TRUE)) # xx
Lgpdelt_africa_pop_raster_points2_table_4 <- ungroup(Lgpdelt_africa_pop_raster_points2_table_3)
africa_population_in_Lgpdelt_sum <- sum(Lgpdelt_africa_pop_raster_points2_table_4$value)


# Creating a raster
Lgpdelt_u0.64ha_raster <- Lgpdelt_resampled
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 111)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[1, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 112)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[2, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 121)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[3, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 122)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[4, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 211)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[5, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 212)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[6, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 221)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[7, 2])
Lgpdelt_u0.64ha_raster[which(Lgpdelt_u0.64ha_raster[] == 222)] <- as.numeric(Lgpdelt_u0.64ha_pop_table_4[8, 2])

# options(scipen = 999)
# plot(Lgpdelt_u0.64ha_raster)
# tmap_options(max.raster = c(plot = 140349575, view = 140349575))
# pal8 <- c("#33A02C", "#B2DF8A", "#FDBF6F", "#1F78B4", "#999999", "#E31A1C", "#E6E6E6", "#A6CEE3")
hist(Lgpdelt_u0.64ha_raster,
  main = "Distribution of vulnerability values",
  xlab = "value", ylab = "Frequency",
  col = "springgreen"
)
maxValue(Lgpdelt_u0.64ha_raster)
unique(Lgpdelt_u0.64ha_raster)

Lgpdelt_u0.64ha_map <- tm_shape(Lgpdelt_u0.64ha_raster) +
  tm_layout(
    legend.frame = TRUE,
    legend.outside = FALSE, legend.position = c("left", "top"),
    main.title = "Total number smallholder farmers per vulnerability domain",
    main.title.size = 1,
    main.title.position = "centre",
    legend.bg.color = "#ffffff",
    legend.text.size = 0.7
  ) +
  tm_raster(
    title = "",
    col = "layer", palette = palette_4, n = 8,
    labels = c(
      "LLL, 45843062", "LLH, 13709238", "LHL, 39989192", "LHH, 54137330",
      "HLL, 99266317", "HLH, 29958274", "HHL, 129520832", "HHH, 43557124"
    ),
    legend.hist = FALSE, legend.is.portrait = TRUE,
    style = "cat"
  ) +
  tmap_style("white") +
  tm_shape(countries_shp_cropped) +
  tm_borders("black", lwd = .5)

# Lgpdelt_u0.64ha_raster_test<-levelplot(Lgpdelt_u0.64ha_raster, margin=FALSE, par.settings = theme3,main="test")



# write files

write(smallholder_population_in_Lgpdelt_sum, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/smallholder_population_in_Lgpdelt_sum.txt")
writeRaster(Lgpdelt_u0.64ha_raster, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/Lgpdelt_u0.64ha_raster.tif", overwrite = TRUE)
writeRaster(Lgpdelt_resampled, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/Lgpdelt_resampled.tif", overwrite = TRUE)
write_csv(Lgpdelt_u0.64ha_pop_table_4, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/Lgpdelt_u0.64ha_pop_table_4.csv")
jpeg("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/Lgpdelt_u0.64ha_map.png")
Lgpdelt_u0.64ha_map
dev.off()
save.image("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/global_env.RData")
totals_df <- data.frame(
  variable = c(
    "WorldPop raster", "WorldPop raster + farming systems",
    "WorldPop raster + countries", "WorldPop raster + Lgpdelt"
  ),
  total_population = c(
    total_africa_pop, farming_systems_population_sum,
    population_country_sum, africa_population_in_Lgpdelt_sum
  ),
  smallholder_population = c(
    total_u0.64ha_pop, farming_systems_smallholder_sum,
    smallholder_population_country_sum, smallholder_population_in_Lgpdelt_sum
  )
)
write.csv(totals_df, "/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/populations_summary.csv")
save.image("/Users/kthickson/Library/Mobile Documents/com~apple~CloudDocs/Katie files 10-07-18.16/5. Career : pursuits/Practice/CIAT 2019/Adaptation atlas/Analysis/global_env_how_many_smallholders_and_where.RData")
