library(neotoma)
library(analogue)
library(Bchron)

# USA Contiguous albers equal area
proj_out <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'


# get datasets withing bounded region
all.datasets <- get_dataset(loc = c(-100, 44, -70, 49),
                            datasettype = "pollen")

map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = get_site(all.datasets),
             aes(x = long, y = lat), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = 44,
            lat1 = 49, 
            xlim = c(-100, -70),
            ylim = c(44, 49))

# if (file.exists("../data/all.downloads.rds")) {
#   all.downloads <- readRDS("../data/all.downloads.rds")
# } else {
  all.downloads <- get_download(all.datasets, verbose = FALSE)
  saveRDS(all.downloads, "../data/all.downloads.rds")
# }

all.cores <- compile_taxa(all.downloads, "P25")
compiled.cores  <- compile_downloads(all.cores)

radio.years <- (compiled.cores$date.type %in% "Radiocarbon years BP") &
  (compiled.cores$age > 71 ) &
  (compiled.cores$age < 46401)
sryears <- sum(radio.years, na.rm = TRUE)
# BChronCalibrate is in the BChron package:
calibrated <- BchronCalibrate(compiled.cores$age[radio.years],
                              ageSds = rep(100, sryears),
                              calCurves = rep("intcal13", sryears))
#  we want the weighted means from "calibrated"
wmean.date <- function(x) sum(x$ageGrid*x$densities / sum(x$densities))

compiled.cores$age[radio.years] <- sapply(calibrated, wmean.date)
# compiled.cores <- na.omit(compiled.cores)

hist(compiled.cores$age)

compiled.cores.sub <- subset(compiled.cores, subset=(age>0)&(age<200))

map <- map_data("world")
ggplot(data = data.frame(map), aes(long, lat)) + 
  geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
  geom_point(data = compiled.cores.sub,
             aes(x = long, y = lat), color = 2) +
  xlab("Longitude West") + 
  ylab("Latitude North") +
  coord_map(projection = "albers",
            lat0 = 44,
            lat1 = 53, 
            xlim = c(-100, -60),
            ylim = c(44, 49))


sp::coordinates(compiled.cores) <- ~long+lat
sp::proj4string(compiled.cores) <- sp::CRS('+init=epsg:4326')
albers <- sp::spTransform(compiled.cores, proj_out)
xy = data.frame(coordinates(albers))
colnames(xy) = c('x', 'y')

compiled.cores = data.frame(xy, compiled.cores)
compiled.cores = compiled.cores[,which(colnames(compiled.cores)!= 'optional')]

compiled.counts = compiled.cores[,13:ncol(compiled.cores)]
compiled.props  = compiled.cores[,13:ncol(compiled.cores)]/rowSums(compiled.cores[,13:ncol(compiled.cores)])
compiled.locs   = compiled.cores[,c('x', 'y')]

# library(maptools)
# source('../source/pie_maps.R')
# 
# xlim = c(-100, -60)
# ylim = c(44, 49)
# 
# xlo = min(xlim)
# xhi = max(xlim)
# ylo = min(ylim)
# yhi = max(ylim)
# 
# shift = 24000
# 
# par(mfrow=c(1,1))
# pieMap(proportions = compiled.props,
#        centers  = compiled.locs,
#        restrict = FALSE,
#        inputRestricted = FALSE,
#        xlim   = c(xlo+shift, xhi-shift),
#        ylim   = c(ylo+shift, yhi-shift),
#        radius = 25000,
#        scale  = 1,
#        xlab   = 'x',
#        ylab   = 'y',
#        add_legend=FALSE,
#        main_title='',
#        cont_shp=cont_shp)


# taxa.keep = c('Acer', 'Alnus','Betula', 'Fagus', 'Ostrya.Carpinus', 'Picea', 'Pinus', 'Quercus', 'Tsuga', 'Ulmus')
taxa.keep = c('Acer', 'Alnus','Betula', 'Cyperaceae', 'Fagus', 'Ostrya.Carpinus', 
              'Picea', 'Pinus', 'Quercus', 'Tsuga', 'Ulmus')
taxa.nontree = c('Other', 'Prairie.Forbs', 'Poaceae')

compiled.other = compiled.counts[, which((!(colnames(compiled.counts) %in% taxa.keep)) & 
                                           (!(colnames(compiled.counts) %in% taxa.nontree)))]

compiled.counts.sub = data.frame(compiled.counts[, which(colnames(compiled.counts) %in% taxa.keep)], 
                 Other = rowSums(compiled.other, na.rm=TRUE))

compiled.props.sub = compiled.counts.sub/rowSums(compiled.counts.sub)


compiled.counts.locs = data.frame(compiled.locs, compiled.counts.sub)

library(data.table)
dat = as.data.table(compiled.counts.locs)[, lapply(.SD, sum, na.rm=TRUE), by = list(x, y)]
dat[,3:ncol(dat)] = round(dat[,3:ncol(dat)])


# mdat = melt(compiled.counts.locs, id.vars=c('x', 'y'))
# dat  = dcast(mdat, x + y ~ variable, function(x) sum(x, na.rm=TRUE))
# dat[,3:ncol(dat)] = round(dat[,3:ncol(dat)])

N_taxa = ncol(compiled.counts.sub)
saveRDS(dat, paste0('../data/', N_taxa, 'taxa_dat.RDS'))

# map <- map_data("world")
# ggplot(data = data.frame(map), aes(long, lat)) + 
#   geom_polygon(aes(group=group), color = "steelblue", alpha = 0.2) +
#   geom_point(data = get_site(all.datasets),
#              aes(x = long, y = lat), color = 2) +
#   xlab("Longitude West") + 
#   ylab("Latitude North") +
#   coord_map(projection = "albers",
#             lat0 = 44,
#             lat1 = 53, 
#             xlim = c(-100, -60),
#             ylim = c(44, 49))
