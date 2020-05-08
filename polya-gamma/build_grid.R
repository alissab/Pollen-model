library(reshape2)
require(rdist)
require(rgeos)
require(ggplot2)
require(sp)
require(rgdal)
library(raster)
require(fields)

#### READ MAP DATA ####
# getting data ready
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
+lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
+towgs84=0,0,0"

na_shp <- readOGR("../data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))
lake_shp <- readOGR("../data/map-data/Great_Lakes.shp", "Great_Lakes")
lake_shp <- sp::spTransform(lake_shp, proj_out)


#### READ IN MODEL DATA AND OUTPUT ####
out = readRDS('polya-gamma-posts_test.RDS')
dat = readRDS('polya-gamma-dat.RDS')

# note that locations were scaled to fit the model
# unscaling to think in meters, then will rescale again before prediction
rescale = dat$rescale
locs_pollen <- dat$locs*rescale 
names(locs_pollen) <- c("x", "y")

y = dat$y

#### CONSTRUCT GRID ####
# function to construct a raster grid
# take in bounding box for grid, resolution in m, and projection
build_grid <- function(veg_box, resolution = 24000, proj = '+init=epsg:3175') {
  raster::raster(xmn = veg_box[1],
                 xmx = veg_box[3],
                 ymn = veg_box[2],
                 ymx = veg_box[4],
                 resolution = resolution,
                 crs = proj)
}

bbox_tran <- function(x, coord_formula = '~ x + y', from, to) {
  
  sp::coordinates(x) <- formula(coord_formula)
  sp::proj4string(x) <- sp::CRS(from)
  bbox <- as.vector(sp::bbox(sp::spTransform(x, CRSobj = sp::CRS(to))))
  return(bbox)
}
# get bounding box from pollen record coordinates
# pol_box <- bbox_tran(locs_pollen, '~ x + y',
#                      '+init=epsg:3175',
#                      '+init=epsg:3175')
pol_box <- bbox_tran(locs_pollen, '~ x + y',
                     proj_out,
                     proj_out)

xlim = c(pol_box[1]-24000, pol_box[3]+24000)
ylim = c(pol_box[2]-24000, pol_box[4]+24000)

# build the raster grid
# 40 km grid cells
reconst_grid <- build_grid(pol_box,
                           resolution = 60000,
                           proj = proj_out)

# want to work with a data frame not a raster
reconst_grid = as.data.frame(reconst_grid, xy=TRUE)

locs_grid = reconst_grid[,1:2]

saveRDS(locs_grid, 'data/grid.RDS')
