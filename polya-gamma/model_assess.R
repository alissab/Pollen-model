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
N_cores = nrow(locs_pollen)

#### DISTANCE MATRICES ####
D_pollen <- fields::rdist(locs_pollen/rescale)# N_cores x N_cores
# D_pollen <- rdist(as.matrix(locs_pollen/rescale))# N_cores x N_cores
any(D_pollen == 0, na.rm = TRUE)   # check if there are off-diagonal zeros
D_pollen <- ifelse(D_pollen == 0, 0.007, D_pollen)  # remove them
diag(D_pollen) <- 0


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
  geom_line(aes(x=distance, y=covar, color=taxon), size=2) +
  theme_bw() +
  xlim(c(0,500)) + 
  xlab("Distance (km)") + 
  ylab("Covariance")
ggsave("../figs/polya-gamma/covariance_vs_distance.pdf")#, device="pdf", type="cairo")

###############################################################################################################################
## trace
###############################################################################################################################

# tau
tau_melt = data.frame(iter=seq(1, N_keep), value=tau)
ggplot() + geom_line(data=tau_melt, aes(x=iter, y=value))

# theta

# mu
mu_melt = melt(mu)
colnames(mu_melt) = c('iter', 'taxon', 'value')
mu_melt$taxon = taxa[mu_melt$taxon]

ggplot(data=mu_melt) + 
  geom_line(aes(x=iter, y=value, color=taxon)) +
  theme_bw()
ggsave('../figs/trace_mu.png', device="png", type="cairo")

###############################################################################################################################
## maps
###############################################################################################################################

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

# pis <- eta_to_pi(etas)

pis = array(NA, dim=c(N_cores, J, N_keep))
for (i in 1:N_keep){
  pis[,,i] <- eta_to_pi(eta[i,,])
  # pis <- pis %>% mutate(sum = rowSums(.))  # check to make sure it worked
}

pi_mean = apply(pis, c(1,2), mean, na.rm=TRUE)
colnames(pi_mean) = c('Acer', 'Alnus','Betula', 'Fagus', 'Ostrya.Carpinus', 'Ulmus')

preds = data.frame(locs_pollen, pi_mean)
preds_melt = melt(preds, id.vars=c('x', 'y'))

props = y/rowSums(y)
dat_melt = melt(data.frame(locs_pollen, props), id.vars=c('x', 'y'))

all_melt = rbind(data.frame(dat_melt, type="data"), 
                 data.frame(preds_melt, type="preds"))

breaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
all_melt$value_binned = cut(all_melt$value, breaks, include.lowest=TRUE, labels=FALSE)

breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })


ggplot() + 
  geom_point(data=all_melt, aes(x=x, y=y, color=factor(value_binned), fill=factor(value_binned))) + 
  scale_fill_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
  scale_colour_manual(values = tim.colors(length(breaks)), labels=breaklabels, name='Proportion', drop=FALSE) + 
  # scale_colour_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  # scale_fill_gradientn(colours = tim.colors(10), limits=c(0,1)) + 
  facet_grid(variable~type) + 
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
ggsave('../figs/all_binned.png', device="png", type="cairo")

###############################################################################################################################
## observed versus predicted
###############################################################################################################################

# foo = dcast(all_melt, x+y+variable~type, fun.aggregate=mean)

foo = merge(preds_melt, dat_melt, by=c("x", "y", "variable"))
ggplot(data=foo) +
  geom_point(aes(x=value.y, y=value.x)) +
  theme_bw() + 
  coord_equal() +
  xlim(c(0,1)) + 
  ylim(c(0,1)) +
  xlab("Observed proportions") +
  ylab("Predicted proportions") +
  geom_abline(intercept=0, slope=1) + 
  facet_wrap(~variable)
ggsave('../figs/props_obs_vs_preds.png', device="png", type="cairo")
