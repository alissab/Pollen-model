
require(raster)
require(geosphere)
require(tidyverse)
require(sp)
require(rgdal)
require(rdist)
require(rgeos)
require(ggplot2)

## code in part from https://romanabashin.com/how-to-generate-regularly-spaced-points-on-a-map-inside-a-polygon/

#### DATA PREP ####
dat <- readRDS("../pollen_data.RData")

proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
  +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
  +towgs84=0,0,0"

na_shp <- readOGR("../NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
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
lakes_shp <- readOGR("../Great_Lakes.shp", "Great_Lakes")
lakes_shp <- sp::spTransform(lakes_shp, proj_out)
shp <- rgeos::gDifference(cont_shp, lakes_shp)
plot(shp)

dat_grid_coords_sp <- SpatialPointsDataFrame(coords = dat_grid_coords[,c('x','y')],
                                             data = dat_grid_coords[,c('site_ref','Alnus',
                                                                       'Betula','Ulmus')])
cropped <- raster::intersect(dat_grid_coords_sp, shp)
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

# REGRESSION MCMC - SAMPLER FUNCTION
mcmc <- function (y, locs, K, message = 100,
                  adapt = K / 2, ## when to stop tuning the MCMC
                  ## default priors
                  alpha_tau = 0.1,
                  beta_tau = 0.1,
                  mean_nu = -1,
                  sd_nu = 1,
                  mean_range = 0,
                  sd_range = 10) {
  require(pgdraw)
  require(mvnfast)
  ## define priors
  theta_mean <- c(mean_range, mean_nu)
  theta_sd <- c(sd_range, sd_nu)
  N <- nrow(y)
  J <- ncol(y)
  D <- fields::rdist(locs)
  ########### Calculate Mi ###################
  Mi <- data.frame(matrix(0, N, J-1))
  sumY <- rep(0, times = N)
  cumsumY <- data.frame(matrix(0, N, J-1))
  
  for(i in 1: N){
    sumY[i] <- sum(y[i, ])
  }
  
  for(i in 1: N){
    cumsumY[i,] <- c(0, cumsum(y[i,1:(J-2)]))
  }
  
  for(i in 1: N){
    Mi[i,] <- sumY[i] - cumsumY[i,]
  }
  
  ####################initialize kappa###################
  kappa <- data.frame(matrix(0, N, J-1))
  for (i in 1: N) {
    kappa[i,] <- y[i, 1:(J-1)] - Mi[i, ] / 2 
  }
  
  ## initialize tau2 using an inverse-gamma(1, 1) prior
  tau2 <- rep(0, K)
  tau <- rep(0, K)
  tau2[1] <- min(1 / rgamma(1, alpha_tau, beta_tau), 10)
  tau[1] <- sqrt(tau2[1])
  
  ## initialize theta
  theta <- matrix(0, K, 2)
  theta[1, ] <- c(rmvn(1, theta_mean, diag(theta_sd)))
  
  ## construct Sigma covariance matrix
  # MAKE SURE IT'S INVERTIBLE
  diag(D) <- NA
  any(D == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
  D <- ifelse(D == 0, 0.007, D)  # remove them
  diag(D) <- 0
  
  Sigma <- tau2[1] * correlation_function(D, theta[1, ])
  Sigma_chol <- chol(Sigma)
  Sigma_inv <- chol2inv(Sigma_chol)
  
  ## don't save eta due to memory constraints for larger datasets
  # eta <- matrix(0, N, J-1)
  eta <- array(0, dim = c(K, N, J-1))
  eta[1, , ] <- t(rmvn(J-1, rep(0, N), Sigma_chol, isChol = TRUE))
  # eta[1, , ] <- t(rmvn(J-1, rep(0, N), Sigma, isChol = FALSE))
  
  
  ##
  ## initialize omega
  ##

  # omega <- matrix(0, N, J-1)
  omega <- array(0, dim=c(K, N, J-1))
  for (i in 1:N) {
    for (j in 1:(J-1)) {
      if (Mi[i, j] != 0) {
        omega[1, i, j] <- pgdraw(Mi[i, j], eta[1, i, j])
      }
    }
  }

  Omega <- vector(mode = "list", length = J-1)
  for (j in 1:(J - 1)) {
    Omega[[j]] <- diag(omega[1, , j])
  }
  
  ##
  ## tuning variables for adaptive MCMC
  ##
  theta_batch <- matrix(0, 50, 2)
  theta_accept <- 0
  theta_accept_batch <- 0
  lambda_theta <- 0.05
  Sigma_theta_tune <- 1.8 * diag(2) - .8
  Sigma_theta_tune_chol <- chol(Sigma_theta_tune)
  
  ## function to update the parameter
  k = K
  accept = theta_accept_batch
  lambda = lambda_theta
  batch_samples = theta_batch
  Sigma_tune = Sigma_theta_tune
  Sigma_tune_chol = Sigma_theta_tune_chol
  
  updateTuningMV <- function(k, accept, lambda, batch_samples,
                             Sigma_tune, Sigma_tune_chol) {
    arr <- c(0.44, 0.35, 0.32, 0.25, 0.234)
    # std::vector<double> acceptance_rates (arr, arr + sizeof(arr) / sizeof(arr[0]))
    dimension <- nrow(batch_samples)
    if (dimension >= 5) {
      dimension <- 5
    }
    d <- ncol(batch_samples)
    batch_size <- nrow(batch_samples)
    optimal_accept <- arr[dimension]
    times_adapted <- floor(k / 50)
    gamma1 <- 1.0 / ((times_adapted + 3.0)^0.8)
    gamma2 <- 10.0 * gamma1
    adapt_factor <- exp(gamma2 * (accept - optimal_accept))
    ## update the MV scaling parameter
    lambda_out <- lambda * adapt_factor
    ## center the batch of MCMC samples
    batch_samples_tmp <- batch_samples
    for (j in 1:d) {
      mean_batch = mean(batch_samples[, j])
      for (i in 1:batch_size) {
        batch_samples_tmp[i, j] <- batch_samples[i, j] - mean_batch
      }
    }
    
    ## 50 is an MCMC batch size, can make this function more general later...
    Sigma_tune_out <- Sigma_tune + gamma1 *
      (t(batch_samples) %*% batch_samples / (50.0-1.0) - Sigma_tune)
    Sigma_tune_chol_out <- chol(Sigma_tune)
    accept_out <- 0.0
    batch_samples_out <- matrix(0, batch_size, d)
    return(list(
      batch_samples = batch_samples_out,
      Sigma_tune = Sigma_tune_out,
      Sigma_tune_chol = Sigma_tune_chol_out,
      lambda = lambda_out,
      accept = accept_out
    ))
  }
  
  message("Starting MCMC, fitting for ", K, " iterations")
  
  for (k in 2:K) {
    print(paste("Iteration", k, sep=" "))
    # if (k %% message == 0) {
    #   message("Fitting iteration ", k, " out of ", K)
    # }
    ##
    ## sample Omega
    ##
    for (i in 1:N) {
      for (j in 1:(J-1)) {
        # if(!is.na(Mi[i, j])) {
        if(Mi[i, j] != 0){
          omega[k, i, j] <- pgdraw(Mi[i, j], eta[k, i, j])
        }
        else {
          omega[k, i, j] <- 0
        }
      }
    }
    
    for (j in 1:(J-1)) {
      Omega[[j]] <- diag(omega[k, , j])
    }
    
    ##
    ## sample eta
    ##
    for (j in 1:(J-1)) {
      ## can make this much more efficient
      Sigma_tilde <- chol2inv(chol(Sigma_inv + Omega[[j]]))
      mu_tilde <- c(Sigma_tilde %*% (Sigma_inv %*% rep(0, N) + kappa[, j]))
      eta[k, , j] <- rmvn(1, mu_tilde, Sigma_tilde)  # k
    }
    
    ##
    ## sample spatial correlation parameters theta
    ##
    theta_star <- rmvn(
      n = 1,
      mu = theta[k-1, ],
      sigma = lambda_theta * Sigma_theta_tune_chol,
      isChol = TRUE
    )
    Sigma_star <- tau[k-1]^2 * correlation_function(D, theta_star)
    
    Sigma_chol_star <- chol(Sigma_star)
    Sigma_inv_star <- chol2inv(Sigma_chol_star)
    
    # mh1 = NA
    mh1 <- sum(
      sapply(
        1:(J-1),
        function(j) {
          dmvn(eta[k, , j], rep(0, N), Sigma_chol_star, isChol = TRUE, log = TRUE)
        }
      )#, na.rm = TRUE  # ADDED THIS ARGUMENT
    ) +
      
      ## prior
      sum(dnorm(theta_star, theta_mean, theta_sd, log = TRUE))
    mh2 <- sum(
      sapply(
        1:(J-1),
        function(j) {
          dmvn(eta[k, , j], rep(0, N), Sigma_chol, isChol = TRUE, log = TRUE)
        }
      )#, na.rm = TRUE  # ADDED THIS
    ) +
      ## prior
      sum(dnorm(theta, theta_mean, theta_sd, log = TRUE))
    mh <- exp(mh1 - mh2)
    
    if(mh > runif(1, 0, 1)) {
      theta[k, ] <- theta_star
      Sigma <- Sigma_star
      Sigma_chol <- Sigma_chol_star
      Sigma_inv <- Sigma_inv_star
      if(k <= adapt) {
        theta_accept_batch <- theta_accept_batch + 1 / 50
      } else {
        theta_accept <- theta_accept + 1 / (K - adapt)
        theta[k, ] <- theta[k-1, ]
      }
    } else {
      theta[k, ] <- theta[k- 1, ]
    }
    
    ## adapt the tuning
    if (k <= adapt) {
      theta_batch[k %% 50, ] <- theta[k, ]#}
      if (k %% 50 == 0) {
        out_tuning <- updateTuningMV(
          k = k,
          accept = theta_accept_batch,
          lambda = lambda_theta,
          batch_samples = theta_batch,
          Sigma_tune = Sigma_theta_tune,
          Sigma_tune_chol = Sigma_theta_tune_chol
        )
        theta_accept_batch <- out_tuning$accept
        lambda_theta_tune <- out_tuning$lambda
        theta_batch <- out_tuning$batch_samples
        Sigma_theta_tune <- out_tuning$Sigma_tune
        Sigma_theta_tune_chol <- out_tuning$Sigma_tune_chol
      }
    }
    
    ##
    ## sample spatial process variance tau^2
    ##
    devs <- kappa - eta[k, , ]
    SS <- sum(devs * (tau[k-1]^2 * Sigma_inv %*% as.matrix(devs)))
    tau2[k] <- 1 / rgamma(1, N * (J - 1) / 2 + alpha_tau,
                          SS / 2 + beta_tau)
    tau[k] <- sqrt(tau2[k])
    Sigma <- tau[2]^2 * correlation_function(D, theta[2, ])
    Sigma_chol <- chol(Sigma)
    Sigma_inv <- chol2inv(Sigma_chol)
  }
  
  message("Acceptance rate for theta is ", theta_accept)
  
  return(
    list(
      tau = tau,
      theta = theta,
      eta = eta,
      omega = omega
    )
  )
}


#### RUNNING THE MODEL ####
locs_scaled = locs/1e6
out <- mcmc(y, locs_scaled, K = 50, message = 100)
# may need to run "mcmc" function multiple times before it works

saveRDS(out, 'tipton_mods_output.RDS')
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
props <- dat[,c("Alnus","Betula","Ulmus")]
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

