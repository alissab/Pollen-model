out = readRDS('tipton_mods_output.RDS')
rescale = 1e6

# getting data ready
proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
  +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
  +towgs84=0,0,0"

na_shp <- readOGR("NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))



# data saved as Y_with_NAs.Rdata
# called "cropped" in above code
Y_with_NAs <- readRDS("Y_with_NAs.RData")
# y <- Y_with_NAs[,c('Alnus','Betula','Ulmus')]
locs_grid <- Y_with_NAs[,c('x','y')]
# K <- 1000

# read in pollen data
dat <- readRDS("pollen_data.RData")
dat_coords <- dat[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)
coords = coordinates(dat_coords_t)

# make regularly spaces points across extent
x_coords <- seq(min(coords[,1]), max(coords[,1]), by = 50000)
y_coords <- seq(min(coords[,2]), max(coords[,2] + 10000), by = 50000)
grid_coords <- data.frame(cbind(rep(x_coords, times = 25), rep(y_coords, each = 70)))
names(grid_coords) <- c("x", "y")
site_coords <- coords
site_coords <- as.matrix(site_coords)
locs_pollen = site_coords[,c('x', 'y')]
y = as.data.frame(dat[,c('Alnus','Betula','Ulmus')])

N_grid = nrow(locs_grid)
N_cores = nrow(locs_pollen)

D_pollen <- fields::rdist(locs_pollen)/rescale # N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0

D_grid   <- fields::rdist(locs_grid)/rescale # N_locs x N_locs
any(D_grid == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_grid <- ifelse(D_grid == 0, 0.007, D_grid)  # remove them
diag(D_grid) <- 0

D_inter   <- fields::rdist(locs_grid, locs_pollen)/rescale # N_locs x N_cores
any(D_inter == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_inter <- ifelse(D_inter == 0, 0.007, D_inter)  # remove them
diag(D_inter) <- 0

N_iter = length(out$tau)
J = dim(out$eta)[3] + 1


########### Calculate Mi ###################
Mi <- data.frame(matrix(0, N_cores, J-1))  # DEFINING COUNTS - REWRITING MULTINOMIAL AS PRODUCT OF BINOMIALS (STICK BREAKING ALGORITHM)
sumY <- rep(0, times = N_cores)
cumsumY <- data.frame(matrix(0, N_cores, J-1))

for(i in 1: N_cores){
  sumY[i] <- sum(y[i, ])
}

for(i in 1: N_cores){
  cumsumY[i,] <- c(0, cumsum(y[i,1:(J-2)]))
}

for(i in 1: N_cores){
  Mi[i,] <- sumY[i] - cumsumY[i,]     # WHAT DOES MI REPRESENT?
}

####################initialize kappa###################
kappa <- data.frame(matrix(0, N_cores, J-1))
for (i in 1: N_cores) {
  kappa[i,] <- y[i, 1:(J-1)] - Mi[i, ] / 2  # WHAT IS KAPPA FOR?
}

tau = out$tau
theta = out$theta
omega = out$omega
eta = out$eta

## calculate the Matern correlation using parameters theta on the log scale
correlation_function <- function(D, theta) {
  geoR::matern(D, exp(theta[1]), exp(theta[2]))
}


eta_preds <- array(0, dim=c(N_grid, J-1, N_iter))
for (i in 1:N_iter){
  print(i)
  Sigma_pol   = tau[i]^2*correlation_function(D_pollen, theta[i, ])
  Sigma_inter = tau[i]^2*correlation_function(D_inter, theta[i, ])
  
  Sigma_pol_chol <- chol(Sigma_pol)  # WHY DO THEY DO THIS?
  Sigma_pol_inv <- chol2inv(Sigma_pol_chol) 
  
  mu = c(0,0)
  
  for (j in 1:(J-1)){
    y_SB = eta[i,,j] - mu[j]#mu_tilde
    
    eta_preds[,j,i] = mu[j] + Sigma_inter %*% solve(Sigma_pol) %*% y_SB  
  }
}


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

pi_preds = array(NA, dim=c(N_grid, J, N_iter))
for (i in 1:N_iter){
  pi_preds[,,i] <- eta_to_pi(eta_preds[,,i])
  # pis <- pis %>% mutate(sum = rowSums(.))  # check to make sure it worked
}

pi_mean = apply(pi_preds, c(1,2), mean, na.rm=TRUE)
colnames(pi_mean) = c('T1', 'T2', 'T3')

library(maptools)
source('pie_maps.R')

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
       add_legend=TRUE,
       main_title='',
       cont_shp=cont_shp)

# plot simulated r values on a map
r_dat = data.frame(locs_grid, pi_mean)
r_melt = melt(r_dat, id.var=c('x', 'y'))

sim_plot <- ggplot(data = r_melt) +
  geom_point(aes(x = x, y = y, fill = value), alpha = 0.7, pch = 21, size = 1) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 1) +
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
