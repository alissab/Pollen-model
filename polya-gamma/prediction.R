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
build_grid <- function(veg_box, resolution = 8000, proj = '+init=epsg:3175') {
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
pol_box <- bbox_tran(locs_pollen, '~ x + y',
                     '+init=epsg:3175',
                     '+init=epsg:3175')

# build the raster grid
# 40 km grid cells
reconst_grid <- build_grid(pol_box,
                           resolution = 40000,
                           proj = '+init=epsg:3175')

# want to work with a data frame not a raster
reconst_grid = as.data.frame(reconst_grid, xy=TRUE)

locs_grid = reconst_grid[,1:2]

# # make regularly spaces points across extent
# x_coords <- seq(min(locs_pollen[,1]), max(locs_pollen[,1]), by = 50000)
# y_coords <- seq(min(locs_pollen[,2]), max(locs_pollen[,2] + 10000), by = 50000)
# locs_grid <- data.frame(cbind(rep(x_coords, times = 25), rep(y_coords, each = 70)))
# names(locs_grid) <- c("x", "y")

N_grid = nrow(locs_grid)
N_cores = nrow(locs_pollen)

#### DISTANCE MATRICES ####
D_pollen <- fields::rdist(locs_pollen/rescale)# N_cores x N_cores
# D_pollen <- rdist(as.matrix(locs_pollen/rescale))# N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0

# D_grid   <- fields::rdist(locs_grid/rescale) # N_locs x N_locs
D_grid   <- fields::rdist(locs_grid/rescale) # N_locs x N_locs
any(D_grid == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_grid <- ifelse(D_grid == 0, 0.007, D_grid)  # remove them
diag(D_grid) <- 0

D_inter   <- fields::rdist(as.matrix(locs_grid/rescale), as.matrix(locs_pollen/rescale)) # N_locs x N_cores
D_inter   <- cdist(as.matrix(locs_grid/rescale), as.matrix(locs_pollen/rescale)) # N_locs x N_cores
any(D_inter == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
# D_inter <- ifelse(D_inter == 0, 0.007, D_inter)  # remove them

foo=as.matrix(locs_grid/rescale)
plot(foo[,1],foo[,2])
bar=as.matrix(locs_pollen/rescale)
points(bar[,1], bar[,2], col="blue", pch=19)

check = apply(D_inter, 2, function(x) which.min(x))

#### PREDICTIONS ####
N_iter = length(out$tau)
J = dim(out$eta)[3] + 1

burn = 20
N_keep = N_iter-burn+1

tau   = out$tau[burn:N_iter]
theta = out$theta[,burn:N_iter,]
omega = out$omega[burn:N_iter,,]
eta   = out$eta[burn:N_iter,,]
mu    = out$mu[burn:N_iter,]

# theta = theta[1,]

# check process estimates
eta_mean = apply(eta, c(2,3), median)
eta_mean = data.frame(locs_pollen, eta_mean)
plot(eta_mean[,'y'], eta_mean[,'X1'])
eta_melt = melt(eta_mean, id.vars=c('x', 'y'))

# this looks okay
ggplot(data=eta_melt) + 
  geom_point(aes(x=x, y=y, fill=value, color=value)) + 
  facet_wrap(~variable) +
  scale_colour_gradientn(colours = terrain.colors(10)) + 
  scale_fill_gradientn(colours = terrain.colors(10)) + 
  coord_equal()

## calculate the Matern correlation using parameters theta on the log scale
correlation_function <- function(D, theta) {
  geoR::matern(D, exp(theta[1]), exp(theta[2]))
}

x = seq(0, max(D_grid), length=1000)
taxa = c('Acer', 'Alnus','Betula', 'Fagus', 'Ostrya.Carpinus', 'Ulmus')
cov_df = data.frame(taxon=character(0),
                    distance=numeric(0),
                    covar=numeric(0))
for (j in 1:(J-1)){
 cov_df = rbind(cov_df, 
                data.frame(taxon = rep(taxa[j], length(x)),
                           distance = x, 
                           covar = mean(tau)^2*geoR::matern(x, exp(colMeans(out$theta[j,,]))[1], exp(colMeans(out$theta[j,,]))[2])))
}


ggplot(data=cov_df) + 
  geom_line(aes(x=distance, y=covar, color=taxon)) +
  theme_bw() +
  xlim(c(0,500))

eta_preds <- array(0, dim=c(N_grid, J-1, N_keep))
for (i in 1:N_keep){
  # for (i in 1:1){
  print(i)
  for (j in 1:(J-1)){
    
    # but tau^2 cancels in the calculation below
    Sigma_pol   = correlation_function(D_pollen, theta[j,i, ])
    Sigma_inter = correlation_function(D_inter, theta[j,i, ])
    Sigma_pol_inv2 <- solve(Sigma_pol) 
  
    y_SB = eta[i,,j] - mu[i,j] 
    #eta_preds[,j,i] = mu2[j] + Sigma_inter %*% Sigma_pol_inv2 %*% y_SB  
    eta_preds[,j,i] = mu[i,j] + Sigma_inter %*% Sigma_pol_inv2 %*% y_SB 
  }
}

# # check spatial correlation
# # want this to not decline as quickly as it does now
# # expect correlation at 500 km
# x = seq(0, 3000, length=500)
# plot(x, correlation_function(x, colMeans(theta)))
# plot(x[1:50], correlation_function(x[1:50], colMeans(theta)))

# checking first iter
eta_preds_real = eta_preds[,,1]
eta_preds_real = data.frame(locs_grid, eta_preds_real)
eta_preds_real_melt = melt(eta_preds_real, id.vars=c('x', 'y'))
ggplot(data=eta_preds_real_melt) + 
  geom_point(aes(x=x, y=y, fill=value, color=value)) + 
  facet_wrap(~variable) +
  scale_colour_gradientn(colours = terrain.colors(10)) + 
  scale_fill_gradientn(colours = terrain.colors(10)) + 
  coord_equal()

# check process estimates on grid
eta_preds_mean = apply(eta_preds, c(1,2), median)
eta_preds_mean = data.frame(locs_grid, eta_preds_mean)

# 
eta_preds_melt = melt(eta_preds_mean, id.vars=c('x', 'y'))

# something is fucked up here
ggplot(data=eta_preds_melt) + 
  geom_point(aes(x=x, y=y, fill=value, color=value)) + 
  facet_wrap(~variable) +
  scale_colour_gradientn(colours = terrain.colors(10)) + 
  scale_fill_gradientn(colours = terrain.colors(10)) + 
  coord_equal()

#### PREDICTED PROPROTIONS ####
# # function to convert eta to pi (proportions)
expit <- function(x) {
  1 / (1 + exp(-x))
}
# 
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

pi_preds = array(NA, dim=c(N_grid, J, N_keep))
for (i in 1:N_keep){
  pi_preds[,,i] <- eta_to_pi(eta_preds[,,i])
  # pis <- pis %>% mutate(sum = rowSums(.))  # check to make sure it worked
}

pi_mean = apply(pi_preds, c(1,2), mean, na.rm=TRUE)
colnames(pi_mean) = c('Acer', 'Alnus','Betula', 'Fagus', 'Ostrya.Carpinus', 'Ulmus')


preds = data.frame(locs_grid, pi_mean)
preds_melt = melt(preds, id.vars=c('x', 'y'))

ggplot(data=preds_melt) + 
  geom_point(aes(x=x, y=y, fill=value, color=value)) + 
  scale_colour_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  scale_fill_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  facet_wrap(~variable) + 
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  geom_path(data = lake_shp, aes(x = long, y = lat, group = group)) +
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
  coord_equal()




library(maptools)
source('../source/pie_maps.R')

xlo = min(locs_grid[,'x'])
xhi = max(locs_grid[,'x'])
ylo = min(locs_grid[,'y'])
yhi = max(locs_grid[,'y'])

shift = 24000

par(mfrow=c(1,1))
pieMap(proportions = pi_mean, 
       centers  = locs_grid,
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 25000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y', 
       add_legend=FALSE,
       main_title='',
       cont_shp=cont_shp)

# # plot simulated r values on a map
# r_dat = data.frame(locs_grid, pi_mean)
# r_melt = melt(r_dat, id.var=c('x', 'y'))
# 
# sim_plot <- ggplot(data = r_melt) +
#   geom_point(aes(x = x, y = y, fill = value), alpha = 0.7, pch = 21, size = 1) +
#   scale_fill_gradient(low = "white", high = "forestgreen") +
#   geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 1) +
#   geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
#   scale_y_continuous(limits = c(300000, 1900000)) + 
#   scale_x_continuous(limits = c(-800000, 2760000)) +
#   labs(fill = "sim data") +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         line = element_blank(),
#         legend.title = element_text(size = 16),
#         legend.text = element_text(size = 14),
#         plot.title = element_blank()) +
#   coord_fixed() +
#   facet_wrap(.~variable)
# print(sim_plot)

# ########### Calculate Mi ###################
# Mi <- data.frame(matrix(0, N_cores, J-1))  # DEFINING COUNTS - REWRITING MULTINOMIAL AS PRODUCT OF BINOMIALS (STICK BREAKING ALGORITHM)
# sumY <- rep(0, times = N_cores)
# cumsumY <- data.frame(matrix(0, N_cores, J-1))
# 
# for(i in 1: N_cores){
#   sumY[i] <- sum(y[i, ])
# }
# 
# for(i in 1: N_cores){
#   cumsumY[i,] <- c(0, cumsum(y[i,1:(J-2)]))
# }
# 
# for(i in 1: N_cores){
#   Mi[i,] <- sumY[i] - cumsumY[i,]     # WHAT DOES MI REPRESENT?
# }
# 
# ####################initialize kappa###################
# kappa <- data.frame(matrix(0, N_cores, J-1))
# for (i in 1: N_cores) {
#   kappa[i,] <- y[i, 1:(J-1)] - Mi[i, ] / 2  # WHAT IS KAPPA FOR?
# }
