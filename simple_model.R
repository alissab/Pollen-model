
require(MASS)
require(runjags)

# modeling simulated data
# can you recover correct params?
dat <- readRDS("pollen_data.RData")
dist <- readRDS("dist_matrix.RData")

# read in North America shape files
na_shp <- readOGR("NA_States_Provinces_Albers.shp", "NA_States_Provinces_Albers")
na_shp <- sp::spTransform(na_shp, proj_out)
cont_shp <- subset(na_shp,
                   (NAME_0 %in% c("United States of America", "Mexico", "Canada")))


Sum_counts <- dat$sum_counts
mu <- rep(10, 139)
NI <- length(mu)
eta <- 2
lambda <- 20

x = seq(0, max(dist), length=100)
y = eta*exp(-x*lambda)
plot(x, y, type="l", xlim=c(0,1))

mcov <- eta * exp(-dist*lambda)
g <- matrix(nrow = 139, ncol = 1)
g <- mvrnorm(mu = mu, Sigma = mcov)
hist(g)

# plot simulated gs to make sure there's spatial correlation
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

g_plot <- data.frame(cbind(g, dat_coords_t@coords))


# plot simulated g values on a map
sim_plot <- ggplot(data = g_plot) +
  geom_point(aes(x = x, y = y, fill = g), alpha = 0.7, pch = 21, size = 6) +
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
  coord_fixed()
print(sim_plot)

cat("
    model{
    # PRIORS
    lambda ~ dunif(5, 50)
    eta ~ dunif(1e-6, 10)
    alpha ~ dnorm(0, 0.001)


    for(i in 1:NI){
      mu[i] <- alpha
    }
 
   
    # COVARIANCE MATRIX
    for(i in 1:NI){
    for(c in i:NI){
    M.cov[i,c] <- eta * exp(-dist[i,c]*lambda)
    }}
    
    for(i in 1:NI){
    for(c in (i+1):NI){
    M.cov[c,i] <- M.cov[i,c]
    }}
    
    M.tau[1:NI,1:NI] <- inverse(M.cov[,])

    # LIKELIHOOD
    g[1:NI] ~ dmnorm(mu[1:NI], M.tau[,])

    #data# NI, dist, g
    #monitor# lambda, eta, alpha
    }", 
    file = 'simple.txt')

simple <- run.jags(model = 'simple.txt', burnin = 2000, sample = 5000, 
                    adapt = 2000, n.chains = 1)#, method = 'parallel')
plot(simple, 'trace', vars = c("lambda", "eta", "alpha"))
simp_summ <- summary(simple)
g_summ <- simp_summ[3:141,]
hist(g_summ[,4])

# model not able to recover eta, lambda

# plot model estimated g values on a map
plot_g <- data.frame(cbind(g_summ[,4], dat_coords_t@coords))
names(plot_g)[1] <- "g"

est_plot <- ggplot(data = plot_g) +
  geom_point(aes(x = x, y = y, fill = g), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) + 
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "model est") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

grid.arrange(sim_plot, est_plot)
