
require(raster)
require(geosphere)
require(tidyverse)
require(sp)
require(rgdal)
require(rdist)
require(rgeos)
require(ggplot2)
require(invgamma)

## code in part from https://romanabashin.com/how-to-generate-regularly-spaced-points-on-a-map-inside-a-polygon/

#### MCMC SAMPLER ####
source('mcmc.R')

#### DATA PREP ####
dat <- readRDS("../data/pollen_data.RData")

proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
  +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
  +towgs84=0,0,0"

na_shp <- readOGR("../data/map-data/NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))

dat_coords <- dat[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)
coords = coordinates(dat_coords_t)
coords <- data.frame(coords)


##### CREATING A GRID #####
# make regularly spaces points across spatial extent which will be used for
# making predictions
x_coords <- seq(min(coords[,1]), max(coords[,1]), by = 50000)
y_coords <- seq(min(coords[,2]), max(coords[,2] + 10000), by = 50000)
grid_coords <- data.frame(cbind(rep(x_coords, times = 25), rep(y_coords, each = 70)))
names(grid_coords) <- c("x", "y")
site_coords <- coords

# associate each site location with the closest of the grid points
# use rdist funtion to calculate distances between grid points and sites
# it'll give you a distance matrix
# find the minimum distance for each site
# use the associated grid coordinate for each site
site_coords <- as.matrix(site_coords)
grid_coords <- as.matrix(grid_coords)
dist <- cdist(X = site_coords, Y = grid_coords)

# for each row (site), which column contains the minimum value?
# the name of that column (minus "V") is the row in the grid_coords you should use

# make sure all distances are less than ~35km, otherwise something is wrong
min_dist <- apply(dist, MARGIN = 1, FUN = "min")

# find column number associated with minimum value for each row
min_dist <- apply(dist, MARGIN = 1, FUN = "which.min")

# associate sites with coordinates from grid_coords that are nearest
grid_coords <- data.frame(grid_coords)
grid_coords$ref <- as.integer(row.names(grid_coords))
min_dist <- data.frame(min_dist)
colnames(min_dist) <- "ref"
min_dist$site_ref <- 1:139
min_coords <- left_join(grid_coords, min_dist, by = "ref")

# create df with all grid points, with observations and NAs
dat$site_ref <- as.integer(row.names(dat))
dat_grid_coords <- left_join(min_coords, dat, by = "site_ref")
dat_grid_coords <- dat_grid_coords[,c('x','y','site_ref','Alnus','Betula','Ulmus')]
dat_grid_coords <- dat_grid_coords %>% group_by(x, y) %>% 
  summarise(site_ref = min(site_ref), Alnus = sum(Alnus), Betula = sum(Betula), 
            Ulmus = sum(Ulmus))

# make sure grid points are correct
site_coords <- data.frame(site_coords)
grid_coords <- data.frame(grid_coords)
new_sites <- left_join(min_dist, grid_coords, by = "ref")

ggplot(data = grid_coords) +
  geom_point(aes(x = x, y = y), alpha = 0.3) +
  geom_point(data = site_coords, aes(x = x, y = y),
             size = 2, alpha = 0.5, col = "red") +
  geom_point(data = new_sites, aes(x = x, y = y),
             size = 2, alpha = 0.7, col = "blue") +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(400000, 2200000)) +
  scale_x_continuous(limits = c(-800000, 3600000)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()


# then clip grid coordinates to land (don't estimate pollen in the ocean or lakes)
lakes_shp <- readOGR("../data/map-data/Great_Lakes.shp", "Great_Lakes")
lakes_shp <- sp::spTransform(lakes_shp, proj_out)
shp <- rgeos::gDifference(cont_shp, lakes_shp)
plot(shp)

dat_grid_coords_sp <- SpatialPointsDataFrame(coords = dat_grid_coords[,c('x','y')],
                                             data = dat_grid_coords[,c('site_ref','Alnus',
                                                                       'Betula','Ulmus')])
cropped <- raster::intersect(dat_grid_coords_sp, shp)
# warning messgae from above line that needs to be dealt with
cropped <- data.frame(cropped)

ggplot(data = cropped) +
  geom_point(aes(x = x, y = y), alpha = 0.3) +
  geom_point(data = site_coords, aes(x = x, y = y),
             size = 2, alpha = 0.5, col = "red") +
  geom_point(data = new_sites, aes(x = x, y = y),
             size = 2, alpha = 0.7, col = "blue") +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lakes_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(400000, 2200000)) +
  scale_x_continuous(limits = c(-800000, 3600000)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()




##### WRITING THE MODEL #####

require(tidyverse)
require(mvnfast)
require(splines)
require(pgdraw)
require(fields)
require(geoR)

K <- 1000
locs = coords[,c('x', 'y')]
y = as.data.frame(dat[,c('Alnus','Betula','Ulmus')])

## calculate the Matern correlation using parameters theta on the log scale
correlation_function <- function(D, theta) {
  geoR::matern(D, exp(theta[1]), exp(theta[2]))
}

#### RUNNING THE MODEL & SAVING OUTPUT####
rescale=1e3
locs_scaled = locs/rescale
out <- mcmc_mu(y,
            locs_scaled,
            K = 5000,
            message = 100,
            mean_nu = -1,
            sd_nu = 0.3,
            mean_range = 5,
            sd_range = 1)
# out <- mcmc(y,
#             locs_scaled,
#             K = 1000,
#             message = 100,
#             mean_nu = -1,
#             sd_nu = 0.3,
#             mean_range = 5,
#             sd_range = 1)
# out <- mcmc_cov(y,
#             locs_scaled,
#             K = 500,
#             message = 100,
#             mean_nu = -1,
#             sd_nu = 0.3,
#             mean_range = 5,
#             sd_range = 1)
saveRDS(out, 'polya-gamma-posts.RDS')

dat = list(y=y,
     locs=locs_scaled,
     rescale=rescale)
saveRDS(dat, 'polya-gamma-dat.RDS')




# out <- readRDS("tipton_mod_with_data.RData")


#### SUMMARIZE PROCESS PARAMETER ####
# taxon 1
t1 <- data.frame(out[["eta"]][,,1])
t1 <- t1 %>% gather()
t1 <- t1 %>% group_by(key) %>% summarise_at("value", 
                                            funs(mean(value), 
                                                 sd(value), 
                                                 quantile(value, 0.025), 
                                                 quantile(value, 0.975)))
t1$key <- gsub("X", "", t1$key)
t1$key <- as.integer(t1$key)
t1 <- arrange(t1, key)

# taxon 2 
t2 <- data.frame(out[["eta"]][,,2])
t2 <- t2 %>% gather()
t2 <- t2 %>% group_by(key) %>% summarise_at("value", 
                                            funs(mean(value), 
                                                 sd(value), 
                                                 quantile(value, 0.025), 
                                                 quantile(value, 0.975)))
t2$key <- gsub("X", "", t2$key)
t2$key <- as.integer(t2$key)
t2 <- arrange(t2, key)


# plot model estimates
require(raster)
require(rgdal)
require(ggplot2)

# # getting data ready
# proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
#   +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# # WGS84
# proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
#   +towgs84=0,0,0"
# na_shp <- readOGR("NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
# na_shp <- sp::spTransform(na_shp, proj_out)
# cont_shp <- subset(na_shp,
#                    (NAME_0 %in% c("United States of America", "Mexico", "Canada")))

# project lat/long from dataframe
# dat <- readRDS("pollen_data.RData")

# dat_coords <- dat[, c("long","lat")]
# names(dat_coords) <- c("x", "y")
# coordinates(dat_coords) <- dat_coords
# sp::proj4string(dat_coords) <- proj_WGS84
# dat_coords_t <- sp::spTransform(dat_coords, proj_out)
# coords = coordinates(dat_coords_t)

dat1 <- data.frame(coords, t1[,c("mean","sd")])
dat_plot1 <- dat1
dat_plot1$mean <- with(dat_plot1, ifelse(mean <= -4, -4, 
                                         ifelse(mean >= 1, 1, mean)))

dat2 <- data.frame(coords, t2[,c("mean","sd")])
dat_plot2 <- dat2
dat_plot2$mean <- with(dat_plot2, ifelse(mean >= 20, 20, mean))

ggplot(data = dat_plot2) +
  geom_point(aes(x = x, y = y, fill = mean), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) +
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Taxon 2\nEta estimate") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()



#### CONVERTING ETAS TO PROPORTIONS AND PLOTTING ####
# use eta summaries (t1 and t2 from above) and combine taxa into one df
eta1 <- t1[,"mean"]
names(eta1) <- "Alnus"
eta2 <- t2[,"mean"]
names(eta2) <- "Betula"
etas <- cbind(eta1, eta2)

# function to convert eta to pi (proportions)
expit <- function(x) {
  1 / (1 + exp(-x))
}

eta_to_pi <- function(eta) {
  ## convert eta to a probability vector pi
  ## can make this more general by first checking if vector vs. matrix and then
  ## calculating the response
  N <- nrow(eta)
  J <- ncol(eta) + 1
  pi <- matrix(0, N, J)
  stick <- rep(1, N)
  for (j in 1:(J - 1)) {
    pi[, j] <- expit(eta[, j]) * stick
    stick <- stick - pi[, j]
  }
  pi[, J] <- stick
  return(pi)
}

pis <- eta_to_pi(etas)
# pis <- pis %>% mutate(sum = rowSums(.))  # check to make sure it worked


# PLOT PROPORTIONS
pis <- data.frame(coords, pis)

alnus <- ggplot(data = pis) +
  geom_point(aes(x = x, y = y, fill = X1), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) +
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Alnus\nproportions") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

betula <- ggplot(data = pis) +
  geom_point(aes(x = x, y = y, fill = X2), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) +
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Betula\nproportions") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

ulmus <- ggplot(data = pis) +
  geom_point(aes(x = x, y = y, fill = X3), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) +
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Ulmus\nproportions") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

require(gridExtra)
grid.arrange(alnus, betula, ulmus, nrow = 3, ncol = 1)


# CALCULATE PROPORTIONS OF OBSERVED COUNTS
props <- dat$y[,c("Alnus","Betula","Ulmus")]
props <- props %>% mutate(sum = rowSums(.))
props <- props %>% mutate(aprop = Alnus/sum, bprop = Betula/sum, uprop = Ulmus/sum)

# PLOT OBSERVED OVER MODEL ESTIMATED PROPORTIONS
par(mfrow = c(3, 1))

plot(pis$X1, props$aprop)
abline(0, 1, col = "red")
plot(pis$X2, props$bprop)
abline(0, 1, col = "red")
plot(pis$X3, props$uprop)
abline(0, 1, col = "red")

