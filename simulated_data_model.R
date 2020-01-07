
require(MASS)
require(runjags)

# modeling simulated data
# can you recover correct params?
dat <- readRDS("pollen_data.RData")
dist <- readRDS("dist_matrix.RData")

Sum_counts <- dat$sum_counts
eta <- 1
lambda <- 0.3
mu <- cbind(rep(-1, 139), rep(0, 139), rep(1, 139))
mcov <- eta * exp(-dist/lambda)
g <- matrix(nrow = 139, ncol = 3)
NI <- nrow(mu)
NJ <- ncol(mu)

for(j in 1:NJ){
  g[1:NI,j] <- mvrnorm(mu = mu[1:NI,j], Sigma = mcov)
}

r1 <- exp(g[,1]) / (exp(g[,1]) + exp(g[,2]) + exp(g[,3]))
r2 <- exp(g[,2]) / (exp(g[,1]) + exp(g[,2]) + exp(g[,3]))
r3 <- exp(g[,3]) / (exp(g[,1]) + exp(g[,2]) + exp(g[,3]))
r <- cbind(r1, r2, r3)

Y <- matrix(nrow = 139, ncol = 3)
for(i in 1:NI){
  Y[i,1:NJ] <- rmultinom(n = 1, prob = r[i,1:NJ], size = Sum_counts[i])
}


cat("
    model{
    # PRIORS
    for(j in 1:NJ){
    spp.mu[j] ~ dnorm(0, 1e-6)
    lambda[j] ~ dunif(0.01, 10)
    eta[j] ~ dunif(1e-6, 10)
    }
    
    for(i in 1:NI){
    for(j in 1:NJ){
    mu[i,j] <- spp.mu[j]
    }}
    
    # COVARIANCE MATRIX
    for(i in 1:NI){
    for(c in i:NI){
    for(j in 1:NJ){
    M.cov[i,c,j] <- eta[j] * exp(-dist[i,c] / lambda[j])
    }}}
    
    for(i in 1:NI){
    for(c in (i+1):NI){
    for(j in 1:NJ){
    M.cov[c,i,j] <- M.cov[i,c,j]
    }}}
    
    M.tau[1:NI,1:NI,1] <- inverse(M.cov[,,1])
    M.tau[1:NI,1:NI,2] <- inverse(M.cov[,,2])
    M.tau[1:NI,1:NI,3] <- inverse(M.cov[,,3])
    
    # LIKELIHOOD
    for(j in 1:NJ){
    g[1:NI,j] ~ dmnorm(mu[1:NI,j], M.tau[,,j])
    }
    
    for(i in 1:NI){
    for(j in 1:NJ){
    r[i,j] <- exp(g[i,j]) / (exp(g[i,1]) + exp(g[i,2]) + exp(g[i,3]))
    }}
    
    for(i in 1:NI){
    Y[i,1:NJ] ~ dmulti(r[i,1:NJ], Sum_counts[i])
    }
    
    #data# Y, Sum_counts, NJ, NI, dist
    #monitor# lambda, eta, g, spp.mu
    }", 
    file = 'simulate.txt')

results <- run.jags(model = 'simulate.txt', burnin = 2000, sample = 5000, 
                      adapt = 2000, n.chains = 3, method = 'parallel')

plot(results, 'trace', vars = c("lambda", "eta", "spp.mu"))
plot(results, 'crosscorr', vars = c("lambda", "eta", "spp.mu"))
summ <- summary(results)

# extract g estimates from summary
g_summ <- summ[grep("g\\[", row.names(summ)), ]
g1 <- g_summ[1:139,]
g2 <- g_summ[140:278,]
g3 <- g_summ[279:417,]

# plot simulated gs against model estimated gs
plot(g[,1], g1[,4], xlab = "simulated g", ylab = "model estimated g")
plot(g[,2], g2[,4])
plot(g[,3], g3[,4])



# PLOTTING PREDICTIONS
require(raster)
require(rgdal)
require(ggplot2)
require(gridExtra)

proj_out <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5
  +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"

# WGS84
proj_WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84
  +towgs84=0,0,0"

# project lat/long from dataframe
dat_coords <- dat[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)

# connect locations with estimates of 'g' for each species
alnus_g <- data.frame(cbind(g1[,4], dat_coords_t@coords))
bet_g <- data.frame(cbind(g2[,4], dat_coords_t@coords))
ulm_g <- data.frame(cbind(g3[,4], dat_coords_t@coords))
colnames(alnus_g)[1] <- "g"
colnames(bet_g)[1] <- "g"
colnames(ulm_g)[1] <- "g"

# read in North America shape files
na_shp <- readOGR("NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))

# plot gs on a map for each taxon
alnus <- ggplot(data = alnus_g) +
  geom_point(aes(x = x, y = y, fill = g), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Alnus \ng param") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

bet <- ggplot(data = bet_g) +
  geom_point(aes(x = x, y = y, fill = g), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Betula \ng param") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

ulm <- ggplot(data = ulm_g) +
  geom_point(aes(x = x, y = y, fill = g), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Ulmus \ng param") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

grid.arrange(alnus, bet, ulm)



# plotting proportions for each taxon
dat_coords <- dat[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)

aln_r <- data.frame(cbind(r[,1], dat_coords_t@coords))
bet_r <- data.frame(cbind(r[,2], dat_coords_t@coords))
ulm_r <- data.frame(cbind(r[,3], dat_coords_t@coords))
colnames(aln_r)[1] <- "r"
colnames(bet_r)[1] <- "r"
colnames(ulm_r)[1] <- "r"

alnus <- ggplot(data = aln_r) +
  geom_point(aes(x = x, y = y, fill = r), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Alnus \nprop (r)") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

bet <- ggplot(data = bet_r) +
  geom_point(aes(x = x, y = y, fill = r), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Betula \nprop (r)") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

ulm <- ggplot(data = ulm_r) +
  geom_point(aes(x = x, y = y, fill = r), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "Ulmus \nprop (r)") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

grid.arrange(alnus, bet, ulm)
