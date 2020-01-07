
require(MASS)
require(runjags)

# modeling simulated data
# can you recover correct params?
dat <- readRDS("pollen_data.RData")
dist <- readRDS("dist_matrix.RData")

Sum_counts <- dat$sum_counts
mu <- rep(10, 139)
NI <- length(mu)
eta <- 0.02
lambda <- 0.3
mcov <- eta * exp(-dist/lambda)
g <- matrix(nrow = 139, ncol = 1)
g <- mvrnorm(mu = mu, Sigma = mcov)
hist(g)

# plot simulated gs to make sure there's spatial correlation
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

cat("
    model{
    # PRIORS
    lambda ~ dunif(0.01, 10)
    eta ~ dunif(1e-6, 10)

    # COVARIANCE MATRIX
    for(i in 1:NI){
    for(c in i:NI){
    M.cov[i,c] <- eta * exp(-dist[i,c] / lambda)
    }}
    
    for(i in 1:NI){
    for(c in (i+1):NI){
    M.cov[c,i] <- M.cov[i,c]
    }}
    
    M.tau[1:NI,1:NI] <- inverse(M.cov[,])

    # LIKELIHOOD
    g[1:NI] ~ dmnorm(mu[1:NI], M.tau[,])

    #data# NI, dist, mu
    #monitor# lambda, eta, g
    }", 
    file = 'simple.txt')

simple <- run.jags(model = 'simple.txt', burnin = 2000, sample = 5000, 
                    adapt = 2000, n.chains = 3, method = 'parallel')
plot(simple, 'trace', vars = c("lambda", "eta"))
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
