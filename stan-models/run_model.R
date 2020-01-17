library(rstan)
require(MASS)
library(reshape2)
require(raster)
require(rgdal)

dat <- readRDS("../pollen_data.RData")
dist <- readRDS("../dist_matrix.RData")


# getting data ready
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
  +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
  +towgs84=0,0,0"
na_shp <- readOGR("../NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))

# project lat/long from dataframe
dat_coords <- dat[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)
coords = coordinates(dat_coords_t)

# prepare data
mu <- rep(10, 139)
eta <- 0.5
lambda <- 20
mcov <- eta*exp(-dist*lambda)
g <- mvrnorm(mu = mu, Sigma = mcov)
data <- data.frame(g, x = dat$long, y = dat$lat)

gp_data <- list(
  N = 139, 
  dist = dist,
  g = g
) 

fit1 <- stan(
  file = "GP.stan",  # Stan program
  data = gp_data,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 2000,            # total number of iterations per chain
  cores = 1,              # number of cores (could use one per chain)
  refresh = 0             # no progress shown
)

print(fit1, pars=c("mu", "eta", "lambda", "lp__"), probs=c(.1,.5,.9))


###########################################################################################################################
# simulating data for model testing
Sum_counts <- dat$sum_counts
eta <- c(2,3,5)
lambda <- c(20, 25, 15)

mu <- cbind(rep(-1, 139), rep(0, 139), rep(1, 139))

mcov=array(NA, c(NJ, NI, NI))
for (j in 1:NJ){
  mcov[j,,] = eta[j] * exp(-dist*lambda[j])
}

# mcov <- eta * exp(-dist*lambda)
g <- matrix(nrow = 139, ncol = 3)
NI <- nrow(mu)
NJ <- ncol(mu)

for(j in 1:NJ){
  g[1:NI,j] <- mvrnorm(mu = mu[1:NI,j], Sigma = mcov[j,,])
}

r1 <- exp(g[,1]) / (exp(g[,1]) + exp(g[,2]) + exp(g[,3]))
r2 <- exp(g[,2]) / (exp(g[,1]) + exp(g[,2]) + exp(g[,3]))
r3 <- exp(g[,3]) / (exp(g[,1]) + exp(g[,2]) + exp(g[,3]))
r <- cbind(r1, r2, r3)

Y <- matrix(nrow = 139, ncol = 3)
for(i in 1:NI){
  Y[i,1:NJ] <- rmultinom(n = 1, prob = r[i,1:NJ], size = Sum_counts[i])
}


gp_data <- list(
  N = 139,
  N_taxa = 3,
  dist = dist,
  y = Y
) 

fit1 <- stan(
  file = "GP_count.stan",  # Stan program
  data = gp_data,    # named list of data
  chains = 1,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 5000,            # total number of iterations per chain
  cores = 1#,              # number of cores (could use one per chain)
  #refresh = 0             # no progress shown
)

print(fit1, pars=c("mu", "eta", "lambda", "lp__"), probs=c(.1,.5,.9))

saveRDS(gp_data, 'gp_data.RDS')
saveRDS(fit1, 'output/GP_count-output.RDS')


r_dat = data.frame(coords, r)
r_melt = melt(r_dat, id.var=c('x', 'y'))


# plot simulated r values on a map
sim_plot <- ggplot(data = r_melt) +
  geom_point(aes(x = x, y = y, fill = value), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "sim data") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed() +
  facet_wrap(.~variable)
print(sim_plot)

# props = r
# colnames(props) = c('alnus', 'betula', 'ulmus')
# dat = data.frame(coordinates(dat_coords_t), region=seq(1,nrow(props)), props)

library(maptools)
source('../pie_maps.R')

xlo = min(coords[,'x'])
xhi = max(coords[,'x'])
ylo = min(coords[,'y'])
yhi = max(coords[,'y'])

shift = 24000

pieMap(proportions = r, 
       centers  = coordinates(dat_coords_t),
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 30000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y', 
       add_legend=TRUE,
       main_title='',
       cont_shp=cont_shp)

pars = rstan::extract(fit1)

g_fit = pars$g
g_fit = t(apply(g_fit, c(2,3), mean))
r_fit = exp(g_fit)/rowSums(exp(g_fit))

pieMap(proportions = r_fit, 
       centers  = coordinates(dat_coords_t),
       restrict = FALSE,
       inputRestricted = FALSE,
       xlim   = c(xlo+shift, xhi-shift),
       ylim   = c(ylo+shift, yhi-shift),
       radius = 30000,
       scale  = 1,
       xlab   = 'x',
       ylab   = 'y', 
       add_legend=TRUE,
       main_title='',
       cont_shp=cont_shp)
