
require(brms)
require(MASS)
require(ggplot2)
require(gridExtra)
require(raster)
require(rgdal)

dat <- readRDS("pollen_data.RData")
dist <- readRDS("dist_matrix.RData")

# simulating data for model testing
mu1 <- rep(10, 139)
mu2 <- rep(5, 139)
mu3 <- rep(15, 139)
eta <- 0.5
lambda <- 20
mcov <- eta*exp(-dist*lambda)
g1 <- mvrnorm(mu = mu1, Sigma = mcov)
g2 <- mvrnorm(mu = mu2, Sigma = mcov)
g3 <- mvrnorm(mu = mu3, Sigma = mcov)

r1 <- exp(g1) / (exp(g1) + exp(g2) + exp(g3))
r2 <- exp(g2) / (exp(g1) + exp(g2) + exp(g3))
r3 <- exp(g3) / (exp(g1) + exp(g2) + exp(g3))
r <- cbind(r1, r2, r3)

NI <- nrow(r)
NJ <- ncol(r)
sum_count <- as.integer(rlnorm(139, meanlog=5, sdlog=1))
Y <- matrix(nrow = 139, ncol = 3)
for(i in 1:NI){
  Y[i,1:NJ] <- rmultinom(n = 1, prob = r[i,1:NJ], size = sum_count[i])
}


# MULTIVARIATE MODEL FOR PROCESS G RESPONSE
mv_mod <- brm(mvbind(g1,g2,g3) ~ gp(x, y), data = all, chains = 1)
summary(mv_mod)

# model predictions
# object pred is a 3-d array with 
# nrow = number of sites
# dimension 2 = stats (mean, sd, CIs)
# dimension 3 = number of taxa
pred_mv <- predict(mv_mod)
p1 <- data.frame(cbind(pred_mv[ ,1,1], dat_coords_t@coords))
names(p1)[1] <- "pred"
p2 <- data.frame(cbind(pred_mv[ ,1,2], dat_coords_t@coords))
names(p2)[1] <- "pred"
p3 <- data.frame(cbind(pred_mv[ ,1,3], dat_coords_t@coords))
names(p3)[1] <- "pred"

plot(p1$pred, g1, xlab = "Predictions", ylab = "Simulated data")
plot(p2$pred, g2, xlab = "Predictions", ylab = "Simulated data")
plot(p3$pred, g3, xlab = "Predictions", ylab = "Simulated data")




# MULTINOMIAL MODEL FOR COUNT RESPONSE
# prepare data
count <- data.frame(c1 = Y[,1] , c2 = Y[,2], c3 = Y[,3], x = dat$long, 
                    y = dat$lat, sum_count = sum_count)
count$c <- with(count, cbind(c1, c2, c3))

# run model
mn_mod <- brm( bf( c | trials(sum_count) ~ gp(x, y)),
            family = multinomial(), data = count, chains = 1)

# summarize results
summary(mn_mod)
pred_mn <- predict(mn_mod)
c1 <- data.frame(cbind(pred_mn[ ,1,1], dat_coords_t@coords))
names(c1)[1] <- "pred"
c2 <- data.frame(cbind(pred_mn[ ,1,2], dat_coords_t@coords))
names(c2)[2] <- "pred"
c3 <- data.frame(cbind(pred_mn[ ,1,3], dat_coords_t@coords))
names(c3)[3] <- "pred"

plot(c1$pred, Y[,1], xlab = "Predictions", ylab = "Simulated data")
plot(c2$pred, Y[,2], xlab = "Predictions", ylab = "Simulated data")
plot(c3$pred, Y[,3], xlab = "Predictions", ylab = "Simulated data")



# UNIVARIATE GAUSSIAN PROCESS USING G PROCESS RESPONSE 
# prepare data
mu <- rep(10, 139)
eta <- 0.5
lambda <- 20
mcov <- eta*exp(-dist*lambda)
g <- mvrnorm(mu = mu, Sigma = mcov)
data <- data.frame(g, x = dat$long, y = dat$lat)

# run model
mod <- brm(g ~ gp(x, y), data = data)

# summarize
pred <- predict(mod)
plot <- data.frame(cbind(pred[,1], dat_coords_t@coords))
names(plot)[1] <- "pred"

# prediction vs. simulated data scatterplot
pred_sim <- data.frame(pred[,1], g)
names(pred_sim) <- c("predictions", "simulated_data")
plot(pred_sim$predictions, pred_sim$simulated_data, 
     xlab = "Predictions", ylab = "Simulated data")



# STAN CODE FOR UNIVARIATE MODEL OF G:
"
functions {
  
  /* compute a latent Gaussian process
  * Args:
    *   x: array of continuous predictor values
  *   sdgp: marginal SD parameter
  *   lscale: length-scale parameter
  *   zgp: vector of independent standard normal variables 
  * Returns:  
    *   a vector to be added to the linear predictor
  */ 
    vector gp(vector[] x, real sdgp, vector lscale, vector zgp) { 
      int Dls = rows(lscale);
      int N = size(x);
      matrix[N, N] cov;
      if (Dls == 1) {
        // one dimensional or isotropic GP
        cov = cov_exp_quad(x, sdgp, lscale[1]);
      } else {
        // multi-dimensional non-isotropic GP
        cov = cov_exp_quad(x[, 1], sdgp, lscale[1]);
        for (d in 2:Dls) {
          cov = cov .* cov_exp_quad(x[, d], 1, lscale[d]);
        }
      }
      for (n in 1:N) {
        // deal with numerical non-positive-definiteness
        cov[n, n] += 1e-12;
      }
      return cholesky_decompose(cov) * zgp;
    }
  
  /* Spectral density function of a Gaussian process
  * Args:
    *   x: array of numeric values of dimension NB x D
  *   sdgp: marginal SD parameter
  *   lscale: vector of length-scale parameters
  * Returns: 
    *   numeric values of the function evaluated at 'x'
  */
    vector spd_cov_exp_quad(vector[] x, real sdgp, vector lscale) {
      int NB = dims(x)[1];
      int D = dims(x)[2];
      int Dls = rows(lscale);
      vector[NB] out;
      if (Dls == 1) {
        // one dimensional or isotropic GP
        real constant = square(sdgp) * (sqrt(2 * pi()) * lscale[1])^D;
        real neg_half_lscale2 = -0.5 * square(lscale[1]);
        for (m in 1:NB) {
          out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
        }
      } else {
        // multi-dimensional non-isotropic GP
        real constant = square(sdgp) * sqrt(2 * pi())^D * prod(lscale);
        vector[Dls] neg_half_lscale2 = -0.5 * square(lscale);
        for (m in 1:NB) {
          out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m])));
        }
      }
      return out;
    }
  /* compute an approximate latent Gaussian process
  * Args:
    *   X: Matrix of Laplacian eigen functions at the covariate values
  *   sdgp: marginal SD parameter
  *   lscale: vector of length-scale parameters
  *   zgp: vector of independent standard normal variables 
  *   slambda: square root of the Laplacian eigen values
  * Returns:  
    *   a vector to be added to the linear predictor
  */ 
    vector gpa(matrix X, real sdgp, vector lscale, vector zgp, vector[] slambda) { 
      vector[cols(X)] diag_spd = sqrt(spd_cov_exp_quad(slambda, sdgp, lscale));
      return X * (diag_spd .* zgp);
    }
}
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  // data related to GPs
  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Kgp_1;
  int<lower=1> Dgp_1;  // GP dimension
  // covariates of the GP
  vector[Dgp_1] Xgp_1[N];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept;
  // GP standard deviation parameters
  vector<lower=0>[Kgp_1] sdgp_1;
  // GP length-scale parameters
  vector<lower=0>[1] lscale_1[Kgp_1];
  // latent variables of the GP
  vector[N] zgp_1;
  real<lower=0> sigma;  // residual SD
}
transformed parameters {
}
model {
  // initialize linear predictor term
  vector[N] mu = Intercept + rep_vector(0, N) + gp(Xgp_1, sdgp_1[1], lscale_1[1], zgp_1);
  // priors including all constants
  target += student_t_lpdf(Intercept | 3, 10, 10);
  target += student_t_lpdf(sdgp_1 | 3, 0, 10)
  - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zgp_1 | 0, 1);
  target += inv_gamma_lpdf(lscale_1[1] | 0.620388, 0.000501);
  target += student_t_lpdf(sigma | 3, 0, 10)
  - 1 * student_t_lccdf(0 | 3, 0, 10);
  // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(Y | mu, sigma);
  }
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
"


# STAN CODE FOR MULTINOMIAL MODEL
"
// generated with brms 2.10.0
functions {
  
  /* multinomial-logit log-PMF
  * Args: 
    *   y: array of integer response values
  *   mu: vector of category logit probabilities
  * Returns:  
    *   a scalar to be added to the log posterior 
  */ 
    real multinomial_logit_lpmf(int[] y, vector mu) {
      return multinomial_lpmf(y | softmax(mu));
    }
  
  /* compute a latent Gaussian process
  * Args:
    *   x: array of continuous predictor values
  *   sdgp: marginal SD parameter
  *   lscale: length-scale parameter
  *   zgp: vector of independent standard normal variables 
  * Returns:  
    *   a vector to be added to the linear predictor
  */ 
    vector gp(vector[] x, real sdgp, vector lscale, vector zgp) { 
      int Dls = rows(lscale);
      int N = size(x);
      matrix[N, N] cov;
      if (Dls == 1) {
        // one dimensional or isotropic GP
        cov = cov_exp_quad(x, sdgp, lscale[1]);
      } else {
        // multi-dimensional non-isotropic GP
        cov = cov_exp_quad(x[, 1], sdgp, lscale[1]);
        for (d in 2:Dls) {
          cov = cov .* cov_exp_quad(x[, d], 1, lscale[d]);
        }
      }
      for (n in 1:N) {
        // deal with numerical non-positive-definiteness
        cov[n, n] += 1e-12;
      }
      return cholesky_decompose(cov) * zgp;
    }
  
  /* Spectral density function of a Gaussian process
  * Args:
    *   x: array of numeric values of dimension NB x D
  *   sdgp: marginal SD parameter
  *   lscale: vector of length-scale parameters
  * Returns: 
    *   numeric values of the function evaluated at 'x'
  */
    vector spd_cov_exp_quad(vector[] x, real sdgp, vector lscale) {
      int NB = dims(x)[1];
      int D = dims(x)[2];
      int Dls = rows(lscale);
      vector[NB] out;
      if (Dls == 1) {
        // one dimensional or isotropic GP
        real constant = square(sdgp) * (sqrt(2 * pi()) * lscale[1])^D;
        real neg_half_lscale2 = -0.5 * square(lscale[1]);
        for (m in 1:NB) {
          out[m] = constant * exp(neg_half_lscale2 * dot_self(x[m]));
        }
      } else {
        // multi-dimensional non-isotropic GP
        real constant = square(sdgp) * sqrt(2 * pi())^D * prod(lscale);
        vector[Dls] neg_half_lscale2 = -0.5 * square(lscale);
        for (m in 1:NB) {
          out[m] = constant * exp(dot_product(neg_half_lscale2, square(x[m])));
        }
      }
      return out;
    }
  /* compute an approximate latent Gaussian process
  * Args:
    *   X: Matrix of Laplacian eigen functions at the covariate values
  *   sdgp: marginal SD parameter
  *   lscale: vector of length-scale parameters
  *   zgp: vector of independent standard normal variables 
  *   slambda: square root of the Laplacian eigen values
  * Returns:  
    *   a vector to be added to the linear predictor
  */ 
    vector gpa(matrix X, real sdgp, vector lscale, vector zgp, vector[] slambda) { 
      vector[cols(X)] diag_spd = sqrt(spd_cov_exp_quad(slambda, sdgp, lscale));
      return X * (diag_spd .* zgp);
    }
}
data {
  int<lower=1> N;  // number of observations
  int<lower=2> ncat;  // number of categories
  int Y[N, ncat];  // response array
  int trials[N];  // number of trials
  // data related to GPs
  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Kgp_muc2_1;
  int<lower=1> Dgp_muc2_1;  // GP dimension
  // covariates of the GP
  vector[Dgp_muc2_1] Xgp_muc2_1[N];
  // data related to GPs
  // number of sub-GPs (equal to 1 unless 'by' was used)
  int<lower=1> Kgp_muc3_1;
  int<lower=1> Dgp_muc3_1;  // GP dimension
  // covariates of the GP
  vector[Dgp_muc3_1] Xgp_muc3_1[N];
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  // temporary intercept for centered predictors
  real Intercept_muc2;
  // GP standard deviation parameters
  vector<lower=0>[Kgp_muc2_1] sdgp_muc2_1;
  // GP length-scale parameters
  vector<lower=0>[1] lscale_muc2_1[Kgp_muc2_1];
  // latent variables of the GP
  vector[N] zgp_muc2_1;
  // temporary intercept for centered predictors
  real Intercept_muc3;
  // GP standard deviation parameters
  vector<lower=0>[Kgp_muc3_1] sdgp_muc3_1;
  // GP length-scale parameters
  vector<lower=0>[1] lscale_muc3_1[Kgp_muc3_1];
  // latent variables of the GP
  vector[N] zgp_muc3_1;
}
transformed parameters {
}
model {
  // initialize linear predictor term
  vector[N] muc2 = Intercept_muc2 + rep_vector(0, N) + gp(Xgp_muc2_1, sdgp_muc2_1[1], lscale_muc2_1[1], zgp_muc2_1);
  // initialize linear predictor term
  vector[N] muc3 = Intercept_muc3 + rep_vector(0, N) + gp(Xgp_muc3_1, sdgp_muc3_1[1], lscale_muc3_1[1], zgp_muc3_1);
  // linear predictor matrix
  vector[ncat] mu[N];
  for (n in 1:N) {
    mu[n] = [0, muc2[n], muc3[n]]';
  }
  // priors including all constants
  target += student_t_lpdf(sdgp_muc2_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zgp_muc2_1 | 0, 1);
  target += inv_gamma_lpdf(lscale_muc2_1[1] | 0.620388, 0.000501);
  target += student_t_lpdf(sdgp_muc3_1 | 3, 0, 10)
    - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(zgp_muc3_1 | 0, 1);
  target += inv_gamma_lpdf(lscale_muc3_1[1] | 0.620388, 0.000501);
  // likelihood including all constants
  if (!prior_only) {
    for (n in 1:N) {
      target += multinomial_logit_lpmf(Y[n] | mu[n]);
    }
  }
}
generated quantities {
  // actual population-level intercept
  real b_muc2_Intercept = Intercept_muc2;
  // actual population-level intercept
  real b_muc3_Intercept = Intercept_muc3;
}
"

# MAPPING MODEL PREDICTIONS FOR UNIVARIATE MODEL
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
# project lat/long from dataframe
dat_coords <- dat[, c("long","lat")]
names(dat_coords) <- c("x", "y")
coordinates(dat_coords) <- dat_coords
sp::proj4string(dat_coords) <- proj_WGS84
dat_coords_t <- sp::spTransform(dat_coords, proj_out)

# plot model predictions vs. simulated data on maps
pred_plot <- ggplot(data = plot) +
  geom_point(aes(x = x, y = y, fill = pred), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) +
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "g pred") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()

sim_plot <- data.frame(cbind(g, dat_coords_t@coords))
sim_plot <- ggplot(data = sim_plot) +
  geom_point(aes(x = x, y = y, fill = g), alpha = 0.7, pch = 21, size = 6) +
  scale_fill_gradient(low = "white", high = "forestgreen") +
  geom_point(aes(x = x, y = y), color = 'black', pch = 21, size = 6) +
  geom_path(data = cont_shp, aes(x = long, y = lat, group = group)) +
  scale_y_continuous(limits = c(300000, 1900000)) +
  scale_x_continuous(limits = c(-800000, 2760000)) +
  labs(fill = "g sim") +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        line = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_blank()) +
  coord_fixed()
grid.arrange(sim_plot, pred_plot)



# when i was trying to fit the model using inverse-distances as a spatial weighting
dist2 <- ifelse(dist < 0.05, 0, dist) # to limit range of inv_dist values
fun <- function(x) 1/x
inv_dist <- fun(dist2)
inv_dist[inv_dist == Inf] <- NA
fun2 <- function(x) x/max(x, na.rm = TRUE)
inv_dist <- fun2(inv_dist)
inv_dist <- ifelse(is.na(inv_dist), 1, inv_dist)
