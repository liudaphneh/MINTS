######################################
# Full MINTS Run for Enrollment Data #
######################################
library(tidyverse)
library(truncnorm)
library(imputeTS)
library(nimble)
library(aspline)
library(splines2)

library(mvtnorm)
library(tmvtnorm)

library(parallel)


##### import data ####
dat.enroll <- read.csv("ner_ger_no_regions_20210805.csv")
dat.enroll$year <- as.numeric(dat.enroll$year)
dat.enroll <- as.data.frame(dat.enroll)

# remove pre-1970 years
dat.enroll <- dat.enroll[which(dat.enroll$year >= 1970), ]

# remove all countries that have no observations
no.data <- dat.enroll %>% group_by(code) %>%
  mutate(num.NA = sum(is.na(NER)), num.NA.GER = sum(is.na(GER))) %>%
  filter(num.NA == length(unique(dat.enroll$year))) %>%
  filter(num.NA.GER == length(unique(dat.enroll$year))) %>%
  select(code) %>%
  unique() %>%
  as.matrix() %>% as.numeric()

dat.enroll <- filter(dat.enroll, !(code %in% no.data))

# correct the two observations for the Seychelles so that NER <= GER
# (likely a data issue due to small population size in denominator;
#   truncate the observed NER to be the observed GER)
dat.enroll$NER[8804] <- dat.enroll$GER[8804]
dat.enroll$NER[8811] <- dat.enroll$GER[8811]

val.1 <- dat.enroll


##### MCMC inputs #####
# if r == 0.8, use a smaller number of knots for asplines; otherwise value of r does not matter
r <- 0.1
f.aspline.degree <- 1

# proposal function for Metropolis-Hastings step
# uses estimated variances and covariances from small univariate run
Sigma.prop <- matrix(c(0.001475944, -0.001511349,
                       -0.001511349,  0.001577437), nrow = 2, byrow = TRUE)
MH.scaling <- 0.02
Sigma.prop <-  MH.scaling *Sigma.prop

n.iter <- 40000
n.chain <- 5
n.burnin <- 35000

##### fit f and h using complete cases #####
training.cc <- filter(val.1, !is.na(NER), !is.na(GER))
x.spline <- training.cc$GER
y.spline <- training.cc$NER

# option for smaller number of knots
small.knots <- round(nrow(training.cc)/4)

# fit A-splines for f
if(r < 0.8){
  aspline.1.list <- tryCatch({
    aspline.1.ridge <- aspline(x.spline, y.spline, degree = f.aspline.degree)
    aspline.1.fit <- lm(y.spline ~ bSpline(x.spline,
                                           knots = aspline.1.ridge$knots_sel[[which.min(aspline.1.ridge$ebic)]],
                                           degree = f.aspline.degree))
    aspline.1.GER <- predict(aspline.1.fit, newdata = data.frame(x.spline = training.cc$GER))
    list(aspline.1.ridge, aspline.1.fit, aspline.1.GER)
  },
  warning = function(cond){
    n.knots <- 11
    aspline.1.ridge <- aspline(x.spline, y.spline, degree = f.aspline.degree,
                               knots = seq(min(x.spline), max(x.spline), length = n.knots)[-c(1, n.knots)])
    aspline.1.fit <- lm(y.spline ~ bSpline(x.spline,
                                           knots = aspline.1.ridge$knots_sel[[which.min(aspline.1.ridge$ebic)]],
                                           degree = f.aspline.degree))
    aspline.1.GER <- predict(aspline.1.fit, newdata = data.frame(x.spline = training.cc$GER))
    return(list(aspline.1.ridge, aspline.1.fit, aspline.1.GER))
  }
  )
  aspline.1.ridge <- aspline.1.list[[1]]
  aspline.1.fit <- aspline.1.list[[2]]
  aspline.1.GER <- aspline.1.list[[3]]
} else{
  # for r = 80%, use smaller number of knots
  n.knots <- small.knots
  aspline.1.ridge <- aspline(x.spline, y.spline, degree = f.aspline.degree,
                             knots = seq(min(x.spline), max(x.spline), length = n.knots)[-c(1, n.knots)])
  aspline.1.fit <- lm(y.spline ~ bSpline(x.spline,
                                         knots = aspline.1.ridge$knots_sel[[which.min(aspline.1.ridge$ebic)]],
                                         degree = f.aspline.degree))
  aspline.1.GER <- predict(aspline.1.fit, newdata = data.frame(x.spline = training.cc$GER))
}

# check fit of estimated aspline curve
# plot(training.cc$GER, training.cc$NER, ylab = "y", xlab = "x", main = "f")
# lines(training.cc$GER[order(training.cc$GER)], aspline.1.GER[order(training.cc$GER)], col = "red")


# fit A-splines for h
abs.diff.NER <- abs(training.cc$NER - aspline.1.GER)
x.resid <- training.cc$GER
y.resid <- abs.diff.NER

n.knots <- ifelse(r == 0.8, small.knots, 42)
resid.ridge <- aspline(x.resid, y.resid, degree = 1,
                       knots = seq(min(x.resid), max(x.resid), length = n.knots)[-c(1, n.knots)])
resid.fit <- lm(y.resid ~ bSpline(x.resid,
                                  knots = resid.ridge$knots_sel[[which.min(resid.ridge$ebic)]],
                                  degree = 1))
resid.GER <- predict(resid.fit, newdata = data.frame(x.resid = training.cc$GER))

# check fit
# plot(training.cc$GER, abs.diff.NER,
#      xlab = "GER", ylab = "Abs(NER - Spline-predicted NER)",
#      main = "Abs NER Residuals from Spline Curve")
# lines(training.cc$GER[order(training.cc$GER)],
#       resid.GER[order(training.cc$GER)], col = "red")


##### set up indicators ####
# indicators of missingness
idx.NER.R <- is.na(val.1$NER)
idx.GER.R <- is.na(val.1$GER)

# indicator for country, not including t = 1
country.idx <- rep(1:length(setdiff(val.1$code, 0)),
                   each = length(unique(val.1$year)) - 1)

# indicator for what country it is, including t = 1
country.idx.t1 <- rep(1:length(setdiff(val.1$code, 0)),
                      each = length(unique(val.1$year)))


##### functions to store one chain's results #####
# function to initialize list that stores one chain's results
# - called at beginning of sample.block.fxn
# - includes observed values for NER and GER and initial values from chain.inits
store.iter.parallel <- function(n.iter, n.burnin, chain.inits, val.1, idx.NER.R, idx.GER.R){
  n <- n.iter - n.burnin + 3

  accept.count <- rep(0, n)
  accept.rate <- rep(0, n)

  rho <- rep(NA, n)
  beta <- rep(NA, n)
  sigma2.N0 <- rep(NA, n)
  sigma2.N <- matrix(NA, nrow = nrow(val.1), ncol = n)
  sigma2.G <-  rep(NA, n)
  tau.G <- sigma2.G

  alpha.c <- matrix(NA, nrow = length(unique(val.1$code)), ncol = n)
  gamma.c <- matrix(NA, nrow = length(unique(val.1$code)), ncol = n)

  mu.drift <- rep(NA, n)
  sigma2.drift <- rep(NA, n)
  mu.0 <- rep(NA, n)
  sigma2.0 <- rep(NA, n)

  NER.0 <- matrix(NA, nrow = length(unique(val.1$code)), ncol = n)
  GER.0 <- matrix(NA, nrow = length(unique(val.1$code)), ncol = n)

  # NER and GER matrices include the observed values
  NER <- matrix(NA, nrow = nrow(val.1), ncol = n)
  GER <- matrix(NA, nrow = nrow(val.1), ncol = n)
  NER[which(!idx.NER.R), ] <- val.1$NER[which(!idx.NER.R)]
  GER[which(!idx.GER.R), ] <- val.1$GER[which(!idx.GER.R)]

  # fill in initial values
  rho[1] <- chain.inits$rho
  beta[1] <- chain.inits$beta
  sigma2.N0[1] <- chain.inits$sigma2.N0
  sigma2.N[,1] <- chain.inits$sigma2.N
  sigma2.G[1] <- chain.inits$sigma2.G
  tau.G[1] <- chain.inits$tau.G
  alpha.c[,1] <- chain.inits$alpha.c
  gamma.c[,1] <- chain.inits$gamma.c
  mu.drift[1] <- chain.inits$mu.drift
  sigma2.drift[1] <- chain.inits$sigma2.drift
  mu.0[1] <- chain.inits$mu.0
  sigma2.0[1] <- chain.inits$sigma2.0
  NER.0[,1] <- chain.inits$NER.0
  GER.0[,1] <- chain.inits$GER.0
  NER[,1] <- chain.inits$NER
  GER[,1] <- chain.inits$GER

  c.results <- list(accept.count = accept.count, accept.rate = accept.rate,
                    rho = rho, sigma2.N0 = sigma2.N0, sigma2.N = sigma2.N,
                    alpha.c = alpha.c, beta = beta, NER = NER,
                    gamma.c = gamma.c, sigma2.G = sigma2.G, tau.G = tau.G, GER = GER,
                    mu.drift = mu.drift, sigma2.drift = sigma2.drift,
                    mu.0 = mu.0, sigma2.0 = sigma2.0,
                    NER.0 = NER.0, GER.0 = GER.0,
                    NER.R = idx.NER.R, GER.R = idx.GER.R)

  return(c.results)
}

# function to set up list that stores initial values for all chains
store.init.parallel <- function(n.chain, val.1){
  n <- 1

  accept.count <- rep(0, n)
  accept.rate <- rep(0, n)

  rho <- rep(NA, n)
  sigma2.N0 <- rep(NA, n)
  sigma2.N <- matrix(NA, nrow = nrow(val.1),
                     ncol = n)
  alpha.c <- matrix(NA, nrow = length(unique(val.1$code)),
                    ncol = n)
  mu.0 <- rep(NA, n)
  sigma2.0 <- rep(NA, n)
  beta <- rep(NA, n)
  NER <- matrix(NA, nrow = nrow(val.1),
                ncol = n)
  NER.0 <- matrix(NA, nrow = length(unique(val.1$code)),
                  ncol = n)

  gamma.c <- matrix(NA, nrow = length(unique(val.1$code)),
                    ncol = n)
  mu.drift <- rep(NA, n)
  sigma2.drift <- rep(NA, n)
  sigma2.G <-  rep(NA, n)
  tau.G <- sigma2.G
  GER <- matrix(NA, nrow = nrow(val.1),
                ncol = n)
  GER.0 <- matrix(NA, nrow = length(unique(val.1$code)),
                  ncol = n)

  c.results <- list(accept.count = accept.count, accept.rate = accept.rate,
                    rho = rho, sigma2.N0 = sigma2.N0, sigma2.N = sigma2.N,
                    alpha.c = alpha.c, beta = beta, NER = NER,
                    gamma.c = gamma.c, sigma2.G = sigma2.G, tau.G = tau.G, GER = GER,
                    mu.drift = mu.drift, sigma2.drift = sigma2.drift,
                    mu.0 = mu.0, sigma2.0 = sigma2.0,
                    NER.0 = NER.0, GER.0 = GER.0)

  return(rep(list(c.results), n.chain))
}



##### set up initial values #####
mcmc.inits <- store.init.parallel(n.chain, val.1)

# inverse for f: given a Y value, find X based on aspline.1.fit
f.inv <- function(y.search){
  abs.resid.search <- abs(aspline.1.GER - y.search)

  # use the first that shows up if there are multiple matches
  min.idx <- which(abs.resid.search == min(abs.resid.search))[1]

  # x value corresponding to closest y value using training set
  return(training.cc$GER[min.idx])
}


ger.init <- rep(NA, length = nrow(val.1))
ner.init <- rep(NA, length = nrow(val.1))
for(nm in setdiff(val.1$code, 0)){
  nm.subset <- filter(val.1, code == nm)
  cd <- nm.subset$code[1]
  init.idx <- which(val.1$code == nm)[which(is.na(nm.subset$GER))]
  ner.init.idx <- which(val.1$code == nm)[which(is.na(nm.subset$NER))]

  # store initial values for this country
  ger.interp <- nm.subset$GER
  ner.interp <- nm.subset$NER

  # observed years
  obs.GER.idx <- which(!is.na(nm.subset$GER))
  obs.NER.idx <- which(!is.na(nm.subset$NER))

  # if 0 data points, use mean of other countries
  if((sum(!is.na(nm.subset$GER)) == 0) & (sum(!is.na(nm.subset$NER)) == 0)){
    ger.interp <- rep(mean(val.1$GER, na.rm = TRUE), nrow(nm.subset))
    ner.interp <- rep(mean(val.1$NER, na.rm = TRUE), nrow(nm.subset))
  }

  if((sum(!is.na(nm.subset$GER)) == 1) & (sum(!is.na(nm.subset$NER)) == 0)){
    # if only 1 data point total (GER), use constant extrapolation
    ger.interp <- rep(nm.subset$GER[which(!is.na(nm.subset$GER))], nrow(nm.subset))
    ner.interp[obs.GER.idx] <-  predict(aspline.1.fit, newdata = data.frame(x.spline = ger.interp[obs.GER.idx]))
    ner.interp <- rep(ner.interp[obs.GER.idx], nrow(nm.subset))

  } else if((sum(!is.na(nm.subset$NER)) == 1) & (sum(!is.na(nm.subset$GER)) == 0)){
    # if only 1 data point total (NER), use constant extrapolation
    ner.interp <- rep(nm.subset$NER[which(!is.na(nm.subset$NER))], nrow(nm.subset))
    ger.interp[obs.NER.idx] <- sapply(ner.interp[obs.NER.idx], f.inv)
    ger.interp <- rep(ger.interp[obs.NER.idx], nrow(nm.subset))

  } else if((sum(!is.na(nm.subset$GER)) != 0) & (sum(!is.na(nm.subset$NER)) == 0)){
    # we have >1 GER, but no NER
    # first set NER points to f(observed GER)
    # then linear interpolate for both
    ner.interp[obs.GER.idx] <- predict(aspline.1.fit, newdata = data.frame(x.spline = ger.interp[obs.GER.idx]))
    ner.interp <- na_interpolation(ner.interp)
    ger.interp <- na_interpolation(ger.interp)

  } else if((sum(!is.na(nm.subset$NER)) != 0) & (sum(!is.na(nm.subset$GER)) == 0)){
    # we have >1 NER, but no GER
    # first set GER points to f.inv(observed NER)
    # then linear interplate for both
    ger.interp[obs.NER.idx] <- sapply(ner.interp[obs.NER.idx], f.inv)
    ger.interp <- na_interpolation(ger.interp)
    ner.interp <- na_interpolation(ner.interp)

  } else if(!((sum(!is.na(nm.subset$GER)) == 0) & (sum(!is.na(nm.subset$NER)) == 0))){
    # otherwise have at least 1 NER and 1 GER for the country
    # identify first year with observed GER or NER
    first.GER <- obs.GER.idx[1]
    first.NER <- obs.NER.idx[1]
    first.obs <- min(first.GER, first.NER)

    # identify last year with observed GER or NER
    last.GER <- obs.GER.idx[length(obs.GER.idx)]
    last.NER <- obs.NER.idx[length(obs.NER.idx)]
    last.obs <- max(last.GER, last.NER)

    # if we only have one year observed for both, constant extrapolation
    if(first.obs == last.obs){
      ger.interp <- rep(nm.subset$GER[which(!is.na(nm.subset$GER))], nrow(nm.subset))
      ner.interp <- rep(nm.subset$NER[which(!is.na(nm.subset$NER))], nrow(nm.subset))

    } else{
      # we have more than one year observed

      #### GER: interpolate and extrapolate
      # if first observation is GER, use value of GER to extrapolate GER
      # if first observation is NER, find f^{-1}(NER) and use that to extrapolate GER
      const <- ifelse(first.obs == first.GER, nm.subset$GER[first.GER], f.inv(nm.subset$NER[first.NER]))
      ger.interp[1:first.obs] <- const

      # if last observation is GER, use value of GER to extrapolate GER
      # if last observation is NER, find f^{-1}(NER) and use that to extrapolate GER
      const <- ifelse(last.obs == last.GER, nm.subset$GER[last.GER], f.inv(nm.subset$NER[last.NER]))
      ger.interp[last.obs:nrow(nm.subset)] <- const

      # linear interpolation for years [first.obs + 1, last.obs - 1]
      if(first.obs < last.obs){
        ger.interp[first.obs:last.obs] <- na_interpolation(ger.interp[first.obs:last.obs])
      }


      #### NER: interpolate and extrapolate
      # if first observation is NER, use value of NER to extrapolate NER
      # if first observation is GER, find f(GER) and use that to extrapolate NER
      const <- ifelse(first.obs == first.NER, nm.subset$NER[first.NER],
                      predict(aspline.1.fit, newdata = data.frame(x.spline = nm.subset$GER[first.GER])))
      ner.interp[1:first.obs] <- const

      # if last observation is NER, use value of NER to extrapolate NER
      # if last observation is GER, find f(GER) and use that to extrapolate NER
      const <- ifelse(last.obs == last.NER, nm.subset$NER[last.NER],
                      predict(aspline.1.fit, newdata = data.frame(x.spline = nm.subset$GER[last.GER])))
      ner.interp[last.obs:nrow(nm.subset)] <- const

      # linear interpolation for years [first.obs + 1, last.obs - 1]
      if(first.obs < last.obs){
        ner.interp[first.obs:last.obs] <- na_interpolation(ner.interp[first.obs:last.obs])
      }
    }
  }

  # require initial values to satisfy NER <= GER + eps
  idx.to.change <- which(ner.interp > (ger.interp + 5))
  ner.interp[idx.to.change] <- ger.interp[idx.to.change] - 0.1

  if(length(which(ner.interp < 0)) > 0){
    idx.to.change <- which(ner.interp < 0)
    ner.interp[idx.to.change] <- 0
  }
  if(length(which(ger.interp < 0)) > 0){
    idx.to.change <- which(ger.interp < 0)
    ger.interp[idx.to.change] <- 0
  }

  ger.init[init.idx] <- ger.interp[which(is.na(nm.subset$GER))]
  ner.init[ner.init.idx] <- ner.interp[which(is.na(nm.subset$NER))]
}


# fill in initial values
set.seed(123)
# SD for initial values of NER and GER
sd.init <- 2
for(c in 1:n.chain){
  # observed values
  mcmc.inits[[c]][["NER"]][which(!idx.NER.R), ] <- val.1$NER[which(!idx.NER.R)]
  mcmc.inits[[c]][["GER"]][which(!idx.GER.R), ] <- val.1$GER[which(!idx.GER.R)]

  # initial values for missing NER
  ner.init.c <- ner.init[which(idx.NER.R)] + rnorm(length(ner.init[which(idx.NER.R)]), mean = 0, sd = sd.init)
  ner.init.c <- ifelse(ner.init.c < 0, 0, ner.init.c)
  mcmc.inits[[c]][["NER"]][which(idx.NER.R), 1] <- ner.init.c

  # initial values for missing GER
  ger.init.c <- ger.init[which(idx.GER.R)] + rnorm(length(ger.init[which(idx.GER.R)]), mean = 0, sd = sd.init)
  ger.init.c <- ifelse(ger.init.c < 0, 0, ger.init.c)
  mcmc.inits[[c]][["GER"]][which(idx.GER.R), 1] <- ger.init.c

  # initial values for time 0 (set at initial values for time 1 plus noise)
  ner.idx <- which(val.1$year == unique(val.1$year)[1])
  mcmc.inits[[c]][["NER.0"]][, 1]  <- mcmc.inits[[c]][["NER"]][ner.idx, 1] + rnorm(length(ner.idx), mean = 0, sd = 2)
  mcmc.inits[[c]][["GER.0"]][, 1]  <- mcmc.inits[[c]][["GER"]][ner.idx, 1] + rnorm(length(ner.idx), mean = 0, sd = 2)

  mcmc.inits[[c]][["rho"]][1] <- rtruncnorm(1, a = 0, b = 1,
                                            mean = 0.5, sd = 0.3)
  mcmc.inits[[c]][["sigma2.N0"]][1] <- var(mcmc.inits[[c]][["NER"]][,1])/var(mcmc.inits[[c]][["GER"]][,1]) + rnorm(1, 0, 0.1)
  mcmc.inits[[c]][["sigma2.N"]][ ,1] <- rnorm(1, 80, 10)
  mcmc.inits[[c]][["alpha.c"]][,1] <- rnorm(length(mcmc.inits[[c]][["alpha.c"]][,1]), 0, 1)
  mcmc.inits[[c]][["beta"]][1] <- runif(1, 0.5, 1)

  mcmc.inits[[c]][["gamma.c"]][,1] <- rnorm(length(mcmc.inits[[c]][["gamma.c"]][,1]), 3, 1)
  mcmc.inits[[c]][["sigma2.G"]][1] <- rnorm(1, var(diff(mcmc.inits[[c]][["GER"]][,1]))/2,
                                            sd = 10)
  mcmc.inits[[c]][["tau.G"]][1] <- 1/mcmc.inits[[c]][["sigma2.G"]][1]
  mcmc.inits[[c]][["mu.drift"]][1] <- rnorm(1, 2, 0.5)
  mcmc.inits[[c]][["sigma2.drift"]][1] <- rnorm(1, 1, 0.1)
  mcmc.inits[[c]][["mu.0"]][1] <- rnorm(1, 2, 0.5)
  mcmc.inits[[c]][["sigma2.0"]][1] <- rnorm(1, 1, 0.1)
}



##### boundaries for NER.0, GER.0  #####
NER.0.up <- rep(100, length(unique(val.1$code)))
GER.0.up <- rep(170, length(unique(val.1$code)))
for(nm in unique(val.1$code)){
  nm.idx <- which(unique(val.1$code) == nm)
  nm.subset <- filter(val.1, code == nm)

  if(sum(!is.na(nm.subset$GER)) > 0){
    GER.0.up[nm.idx] <- min(nm.subset$GER, na.rm = TRUE)
  }
  if(sum(!is.na(nm.subset$NER)) > 0){
    NER.0.up[nm.idx] <- min(nm.subset$NER, na.rm = TRUE)
  }
}


#### set up priors for GER model ####
train.GER.diff <- val.1 %>%
  group_by(code) %>%
  mutate(diff_G = c(0, diff(GER))) %>%
  filter(year > 1970) %>%
  ungroup() %>%
  select(code, diff_G)

# sample variance, sample mean, and SE(sample mean) for differences in GER
train.GER.var <- var(train.GER.diff$diff_G, na.rm = TRUE)
train.GER.mean <- mean(train.GER.diff$diff_G, na.rm = TRUE)
train.GER.mean.SE <- sd(train.GER.diff$diff_G, na.rm = TRUE)/sqrt(sum(!is.na(train.GER.diff$diff_G)))

# country-specific means
train.GER.mean.c <- train.GER.diff %>% group_by(code) %>%
  mutate(n_c = sum(!is.na(diff_G))) %>%
  filter(n_c > 0) %>%
  mutate(mean_c = mean(diff_G, na.rm = TRUE)) %>%
  select(code, mean_c) %>% unique()

# prior for sigma2.G, parameterized as precision using Gamma instead of InvGamma
alpha.G.prior <- 2
beta.G.prior <- train.GER.var

mu.drift.mu.prior <- train.GER.mean
mu.drift.var.prior <- train.GER.mean.SE^2
sigma2.drift.alpha.prior <- 2
sigma2.drift.beta.prior <- var(train.GER.mean.c$mean_c)


#### set up priors for NER model ####
train.NER.diff <- val.1 %>%
  group_by(code) %>%
  mutate(diff_N = c(0, diff(NER))) %>%
  filter(year > 1970) %>%
  ungroup() %>%
  select(code, diff_N)

# sample variance, sample mean, and SE(sample mean) for differences in NER
train.NER.var <- var(train.NER.diff$diff_N, na.rm = TRUE)
train.NER.mean <- mean(train.NER.diff$diff_N, na.rm = TRUE)
train.NER.mean.SE <- sd(train.NER.diff$diff_N, na.rm = TRUE)/sqrt(sum(!is.na(train.NER.diff$diff_N)))

# country-specific means
train.NER.mean.c <- train.NER.diff %>% group_by(code) %>%
  mutate(n_c = sum(!is.na(diff_N))) %>%
  filter(n_c > 0) %>%
  mutate(mean_c = mean(diff_N, na.rm = TRUE)) %>%
  select(code, mean_c) %>% unique()

# prior for sigma2.N0
alpha.prior <- 2
beta.prior <- train.NER.var

# priors for mu.0
mu.0.mu.prior <- 0
mu.0.var.prior <- train.NER.mean.SE^2
sigma2.0.alpha.prior <- 2
sigma2.0.beta.prior <- var(train.GER.mean.c$mean_c)


#### set up priors for (N0, G0) ####
early.cc <- filter(val.1, year <= 1980, !is.na(NER), !is.na(GER))
# if not enough CC, use less restrictive year range
if(nrow(early.cc) < 3){
  early.cc <- filter(val.1, year <= 2000, !is.na(NER), !is.na(GER))
}
if(nrow(early.cc) < 3){
  early.cc <- filter(val.1, year <= 2010, !is.na(NER), !is.na(GER))
}

Mu.0 <- c(mean(early.cc$NER), mean(early.cc$GER))
Sigma.0 <- matrix(c(var(early.cc$NER),
                    cov(early.cc$NER, early.cc$GER),
                    cov(early.cc$NER, early.cc$GER),
                    var(early.cc$GER)), nrow = 2, byrow = TRUE)
cor.0 <- Sigma.0[1,2]/(sqrt(Sigma.0[1,1])*sqrt(Sigma.0[2,2]))
sigma1.0 <- sqrt(Sigma.0[1,1])
sigma2.0 <- sqrt(Sigma.0[2,2])


#### time indices ####
# indices for last time point T, out of val.1 rows
T.idx <- which(val.1$year == unique(val.1$year)[length(unique(val.1$year))])

# indices for first time point t = 1, out of val.1 rows
first.idx <- which(val.1$year == unique(val.1$year)[1])

#### inputs for sample.block.fxn ####
# are we modeling heteroscedasticity?
het <- TRUE

# for h(GER) to make sure we don't get negative values or values very close to 0 since we divide by h
resid.fit.obs <- predict(resid.fit, x.resid = training.cc$GER)
min.obs <- min(resid.fit.obs[which(resid.fit.obs > 0)])
min.obs <- ifelse(min.obs < 0.01, 0.01, min.obs)

cd.unique <- unique(val.1$code)
yr.unique <- unique(val.1$year)
ger.miss.idx <- which(idx.GER.R)
ner.miss.idx <- which(idx.NER.R)
cd.idx.0 <- seq(1, nrow(val.1), by = length(yr.unique))
cd.alpha.idx.0 <- 1:length(cd.unique)


#### sampling function ####
# returns one chain's values
sample.block.fxn <- function(c){
  chain.inits <- mcmc.inits[[c]]
  chain.results <- store.iter.parallel(n.iter, n.burnin, chain.inits,
                                       val.1, idx.NER.R, idx.GER.R)

  for(iter in 2:n.iter){
    # set up i:
    if(iter == 2){
      i <- 2
      i.last <- 1

    } else if(iter == 3){
      i <- 3
      i.last <- 2
    } else if(iter <= n.burnin){
      # if we're past first iteration but still in burn-in period, store this as i = 3 repeatedly
      i <- 3
      i.last <- 3

    } else{
      # e.g. iter = n.burnin + 1, stored as i = 4
      i <- iter - n.burnin + 3
      i.last <- i-1
    }

    x <- chain.results[["NER"]][ ,i.last]
    x.GER <- chain.results[["GER"]][ ,i.last]

    f.x.GER <- suppressWarnings(predict(aspline.1.fit, data.frame(x.spline = x.GER)))
    # for extrapolation outside the range of NER, truncate
    f.x.GER[which(f.x.GER < 0)] <- 0
    f.x.GER[which(f.x.GER > 100)] <- 100

    alpha.c.vector <- sapply(country.idx, function(x){chain.results[["alpha.c"]][x, i.last]})
    gamma.c.vector <- sapply(country.idx, function(x){chain.results[["gamma.c"]][x, i.last]})
    alpha.c.vector.t1 <- sapply(country.idx.t1, function(x){chain.results[["alpha.c"]][x, i.last]})
    gamma.c.vector.t1 <- sapply(country.idx.t1, function(x){chain.results[["gamma.c"]][x, i.last]})

    rho <- chain.results[["rho"]][i.last]
    sigma2.N0 <- chain.results[["sigma2.N0"]][i.last]
    sigma2.N <- chain.results[["sigma2.N"]][,i.last]
    beta <- chain.results[["beta"]][i.last]
    sigma2.G <- chain.results[["sigma2.G"]][i.last]
    tau.G <- chain.results[["tau.G"]][i.last]
    mu.drift <- chain.results[["mu.drift"]][i.last]
    sigma2.drift <- chain.results[["sigma2.drift"]][i.last]
    mu.0 <- chain.results[["mu.0"]][i.last]
    sigma2.0 <- chain.results[["sigma2.0"]][i.last]
    NER.0 <- chain.results[["NER.0"]][,i.last]
    GER.0 <- chain.results[["GER.0"]][,i.last]


    ########### STEP 1: GLOBAL PARAMETERS ###########
    ############ sample sigma2.N
    # find h(G): to avoid negative and 0 values from A-spline, truncate to smallest positive h(GER) from observed data
    if(het){
      h.GER <- suppressWarnings(predict(resid.fit, data.frame(x.resid = x.GER)))
      h.GER[which(h.GER <= min.obs)] <- min.obs
    } else{
      h.GER <- rep(1, length(x.GER))
    }

    sigma2.N0 <- rinvgamma(n = 1, shape = alpha.prior + (length(x)/2),
                           scale = beta.prior + 0.5*sum(((x[-first.idx] - rho*x[-T.idx] -
                                                            alpha.c.vector - beta*f.x.GER[-first.idx])^2) / h.GER[-first.idx]) +
                             0.5*sum((x[first.idx] - rho*NER.0 - alpha.c.vector.t1[first.idx] -
                                        beta*f.x.GER[first.idx])^2  / h.GER[first.idx]  ) )
    sigma2.N <- sigma2.N0*h.GER
    chain.results[["sigma2.N0"]][i] <- sigma2.N0
    chain.results[["sigma2.N"]][,i] <- sigma2.N


    ############ block update (beta, rho)
    br.new <- rtmvnorm(1, mean = c(beta, rho),
                       sigma = Sigma.prop,
                       lower = c(-Inf, 0),
                       upper = c(Inf, 1))
    beta.new <- br.new[,1]
    rho.new <- br.new[,2]

    exp.new <- sum((1/sigma2.N[-first.idx])*(x[-first.idx] - alpha.c.vector - rho.new*x[-T.idx] - beta.new*f.x.GER[-first.idx])^2) + sum((1/sigma2.N[first.idx])*(x[first.idx] - alpha.c.vector.t1[first.idx] - rho.new*NER.0 - beta.new*f.x.GER[first.idx])^2)
    exp.old <- sum((1/sigma2.N[-first.idx])*(x[-first.idx] - alpha.c.vector - rho*x[-T.idx] - beta*f.x.GER[-first.idx])^2) + sum((1/sigma2.N[first.idx])*(x[first.idx] - alpha.c.vector.t1[first.idx] - rho*NER.0 - beta*f.x.GER[first.idx])^2)

    alpha.br  <- exp(-0.5 * (exp.new - exp.old))
    if(alpha.br > 1){ alpha.br <- 1 }
    if(alpha.br >= runif(1)){
      beta <- beta.new
      rho <- rho.new
      chain.results[["accept.count"]][i] <- chain.results[["accept.count"]][i.last] + 1
    } else{
      chain.results[["accept.count"]][i] <- chain.results[["accept.count"]][i.last]
    }
    chain.results[["accept.rate"]][i] <- chain.results[["accept.count"]][i] / iter

    chain.results[["beta"]][i] <- beta
    chain.results[["rho"]][i] <- rho

    ############ sample sigma2.G
    tau.G <- rgamma(n = 1, shape = alpha.G.prior + (length(x.GER)/2),
                    rate = beta.G.prior + 0.5*sum((x.GER[-first.idx] - x.GER[-T.idx] - gamma.c.vector)^2) +
                      0.5*sum((x.GER[first.idx] - GER.0 - gamma.c.vector.t1[first.idx])^2))
    sigma2.G <- 1/tau.G
    chain.results[["tau.G"]][i] <- tau.G
    chain.results[["sigma2.G"]][i] <- sigma2.G


    ############ sample mu_drift, sigma^2_drift
    mu.drift.var.post <- ( (1/mu.drift.var.prior) + (length(cd.unique)/sigma2.drift) )^(-1)
    mu.drift.mu.post <- mu.drift.var.post*( (mu.drift.mu.prior/mu.drift.var.prior) +
                                              sum(chain.results[["gamma.c"]][, i-1])/sigma2.drift )
    mu.drift <- rnorm(1, mean = mu.drift.mu.prior, sd = sqrt(mu.drift.var.post))
    chain.results[["mu.drift"]][i] <- mu.drift

    sigma2.drift <- rinvgamma(n = 1, shape = sigma2.drift.alpha.prior + 0.5*length(cd.unique),
                              scale = sigma2.drift.beta.prior + 0.5*sum((chain.results[["gamma.c"]][, i-1] - mu.drift)^2) )
    chain.results[["sigma2.drift"]][i] <- sigma2.drift


    ############ sample mu_0, sigma^2_0
    mu.0.var.post <- ( (1/mu.0.var.prior) + (length(cd.unique)/sigma2.0) )^(-1)
    mu.0.mu.post <- mu.0.var.post*( (mu.0.mu.prior/mu.0.var.prior) + sum(chain.results[["alpha.c"]][, i-1])/sigma2.0 )
    mu.0 <- rnorm(1, mean = mu.0.mu.post, sd = sqrt(mu.0.var.post))
    chain.results[["mu.0"]][i] <- mu.0

    sigma2.0 <- rinvgamma(n = 1, shape = sigma2.0.alpha.prior + 0.5*length(cd.unique),
                          scale = sigma2.0.beta.prior + 0.5*sum((chain.results[["alpha.c"]][, i-1] - mu.0)^2) )
    chain.results[["sigma2.0"]][i] <- sigma2.0


    ########### STEP 2: COUNTRY-SPECIFIC PARAMETERS ###########
    ############ sample alpha_c and gamma_c
    for(cd in cd.unique){
      cd.idx <- which(val.1$code == cd)
      cd.alpha.idx <- which(cd.unique == cd)
      cd.x <- x[cd.idx]
      cd.x.GER <- x.GER[cd.idx]
      cd.f.x.GER <- f.x.GER[cd.idx]
      cd.first.idx <- first.idx[cd.alpha.idx]
      cd.T.idx <- T.idx[cd.alpha.idx]

      ##### sample alpha_c
      alpha.var.post <- ( (1/sigma2.0) + sum((1/sigma2.N[cd.idx])) )^(-1)
      alpha.mu.post <- alpha.var.post * ( (mu.0/sigma2.0)  + sum( (1/sigma2.N[cd.idx[-1]])*(cd.x[-1] - rho*cd.x[-length(cd.x)] - beta*cd.f.x.GER[-1]) ) + sum( (1/sigma2.N[cd.idx[1]])*(cd.x[1] - rho*NER.0[cd.alpha.idx] - beta*cd.f.x.GER[1]) ) )

      chain.results[["alpha.c"]][cd.alpha.idx, i] <- rnorm(1, mean = alpha.mu.post,
                                                           sd = sqrt(alpha.var.post))

      ##### sample gamma_c
      gamma.var.post <- ( (1/sigma2.drift) + length(yr.unique)/sigma2.G )^(-1)
      gamma.mu.post <- gamma.var.post * ( (mu.drift/sigma2.drift) + sum(cd.x.GER[-1] - cd.x.GER[-length(cd.x.GER)])/sigma2.G + sum(cd.x.GER[1] - GER.0[cd.alpha.idx])/sigma2.G )
      chain.results[["gamma.c"]][cd.alpha.idx, i] <- rnorm(1, mean = gamma.mu.post,
                                                           sd = sqrt(gamma.var.post))
    }


    ########### STEP 3: IMPUTATION ###########
    ############ impute GER
    # draw missing values of GER from imputation model in descending order of time
    for(j in ger.miss.idx[length(ger.miss.idx):1]){
      cd <- val.1$code[j]
      cd.alpha.idx <- which(cd.unique == cd)

      if(j %in% T.idx){
        chain.results[["GER"]][j, i] <- rtruncnorm(1, a = 0, b = Inf,
                                                   mean = x.GER[j-1] + chain.results[["gamma.c"]][cd.alpha.idx, i],
                                                   sd = sqrt(sigma2.G))

      } else if(!(j %in% first.idx)){
        chain.results[["GER"]][j, i] <- rtruncnorm(1, a = 0, b = Inf,
                                                   mean = 0.5*(x.GER[j-1] + chain.results[["GER"]][j+1, i]),
                                                   sd = sqrt(sigma2.G/2))
      } else {
        chain.results[["GER"]][j, i] <- rtruncnorm(1, a = 0, b = Inf,
                                                   mean = 0.5*(GER.0[cd.alpha.idx] + chain.results[["GER"]][j+1, i]),
                                                   sd = sqrt(sigma2.G/2))
      }
    }

    #### impute GER_0
    GER.0.post.var <- ( (1/sigma2.G) + (1/Sigma.0[2,2]) )^(-1)
    GER.0.post.mean <- GER.0.post.var * ( (chain.results[["GER"]][cd.idx.0, i] - chain.results[["gamma.c"]][cd.alpha.idx.0, i])/sigma2.G + (Mu.0[2]/Sigma.0[2,2]) )

    chain.results[["GER.0"]][ , i] <- rtruncnorm(1, a = 0, b = GER.0.up,
                                                 mean = GER.0.post.mean,
                                                 sd = sqrt(GER.0.post.var))
    # rtruncnorm gives NA when a = b so treat separately
    check.bounds <- GER.0.up == 0
    if(sum(check.bounds) > 0){
      chain.results[["GER.0"]][which(check.bounds), i] <- 0
    }

    GER.0 <- chain.results[["GER.0"]][,i]


    ############ impute NER
    # draw missing values of NER from imputation model, sampling in descending order of time
    f.GER.fit <- suppressWarnings(predict(aspline.1.fit, data.frame(x.spline = chain.results[["GER"]][, i])))
    # for extrapolation outside the range of NER, truncate
    f.GER.fit[which(f.GER.fit < 0)] <- 0
    f.GER.fit[which(f.GER.fit > 100)] <- 100

    for(j in ner.miss.idx[length(ner.miss.idx):1]){
      cd <- val.1$code[j]
      cd.alpha.idx <- which(cd.unique == cd)

      if(j %in% T.idx){
        cd.f.GER <- as.numeric(f.GER.fit[j])

        chain.results[["NER"]][j, i] <- rtruncnorm(1, a = 0,
                                                   b = min(chain.results[["GER"]][j,i], 100),
                                                   mean = rho*x[j-1] + chain.results[["alpha.c"]][cd.alpha.idx , i] + beta*cd.f.GER,
                                                   sd = sqrt(sigma2.N[j]))

      } else if(!(j %in% first.idx)){
        cd.f.GER <- as.numeric(f.GER.fit[j])
        cd.f.GER.tp1 <- as.numeric(f.GER.fit[j+1])

        NER.var.post <- ((rho^2/sigma2.N[j+1]) + (1/sigma2.N[j]))^(-1)
        NER.mean.post <- NER.var.post*( ((rho*(chain.results[["NER"]][j+1, i] - chain.results[["alpha.c"]][cd.alpha.idx, i] - beta*cd.f.GER.tp1))/sigma2.N[j+1]) + ((rho*x[j-1] + chain.results[["alpha.c"]][cd.alpha.idx, i] + beta*cd.f.GER)/sigma2.N[j]) )

        chain.results[["NER"]][j, i] <- rtruncnorm(1, a = 0,
                                                   b = min(chain.results[["GER"]][j,i], 100),
                                                   mean = NER.mean.post,
                                                   sd = sqrt(NER.var.post))

        # for countries like Oman where GER is 0, rtruncnorm gives NA when a = b so treat separately
        if(min(chain.results[["GER"]][j,i], 100) == 0){
          chain.results[["NER"]][j, i] <- 0
        }

      } else {
        cd.f.GER <- as.numeric(f.GER.fit[j])
        cd.f.GER.tp1 <- as.numeric(f.GER.fit[j+1])

        NER.var.post <- ((rho^2/sigma2.N[j+1]) + (1/sigma2.N[j]))^(-1)
        NER.mean.post <- NER.var.post*( ((rho*(chain.results[["NER"]][j+1, i] -
                                                 chain.results[["alpha.c"]][cd.alpha.idx, i] -
                                                 beta*cd.f.GER.tp1))/sigma2.N[j+1]) +
                                          ((rho*NER.0[cd.alpha.idx] + chain.results[["alpha.c"]][cd.alpha.idx, i] +
                                              beta*cd.f.GER)/sigma2.N[j]) )

        chain.results[["NER"]][j, i] <- rtruncnorm(1, a = 0,
                                                   b = min(chain.results[["GER"]][j,i], 100),
                                                   mean = NER.mean.post,
                                                   sd = sqrt(NER.var.post))

        # for countries like Oman where GER is 0, rtruncnorm gives NA when a = b so treat separately
        if(min(chain.results[["GER"]][j,i], 100) == 0){
          chain.results[["NER"]][j, i] <- 0
        }
      } # end else statement

    } # end NER for loop

    #### impute NER_0
    cd.f.GER.tp1 <- suppressWarnings(predict(aspline.1.fit,
                                             data.frame(x.spline = chain.results[["GER"]][cd.idx.0, i])))
    # for extrapolation outside the range of NER, truncate
    cd.f.GER.tp1[which(cd.f.GER.tp1 < 0)] <- 0
    cd.f.GER.tp1[which(cd.f.GER.tp1 > 100)] <- 100

    NER.0.post.var <- ( (rho^2 / sigma2.N[cd.idx.0]) + (1 / ((1-cor.0^2)*Sigma.0[1,1]) ) )^(-1)
    NER.0.post.mean <- NER.0.post.var * ( rho*(chain.results[["NER"]][cd.idx.0, i] - chain.results[["alpha.c"]][cd.alpha.idx.0 , i] - beta*cd.f.GER.tp1)/sigma2.N[cd.idx.0]  + (Mu.0[1] + (sigma1.0/sigma2.0) * cor.0 * (GER.0[cd.alpha.idx.0] - Mu.0[2]) )/((1 - cor.0^2)*Sigma.0[1,1] )   )

    chain.results[["NER.0"]][ , i] <- rtruncnorm(1, a = 0,
                                                 b = apply(cbind(GER.0, NER.0.up), 1, min),
                                                 mean = NER.0.post.mean,
                                                 sd = sqrt(NER.0.post.var))
    # rtruncnorm gives NA when a = b so treat separately
    check.bounds <- apply(cbind(GER.0, NER.0.up), 1, min) == 0
    if(sum(check.bounds) > 0){
      chain.results[["NER.0"]][which(check.bounds), i] <- 0
    }

  }

  return(chain.results)
}


##### run model in parallel #####
cl <- makeCluster(n.chain)

clusterEvalQ(cl, {
  library(tidyverse)
  library(truncnorm)
  library(imputeTS)
  library(nimble)
  library(aspline)
  library(splines2)

  library(mvtnorm)
  library(tmvtnorm)
})

clusterExport(cl,
              varlist = c("val.1", "mcmc.inits",
                          "idx.GER.R", "idx.NER.R", "ger.miss.idx", "ner.miss.idx",
                          "n.iter", "n.burnin",
                          "store.iter.parallel", "sample.block.fxn",
                          "country.idx", "country.idx.t1",
                          "aspline.1.fit", "resid.fit",
                          "Sigma.prop", "MH.scaling",
                          "resid.fit.obs", "min.obs", "het",
                          "alpha.prior", "beta.prior",
                          "alpha.G.prior", "beta.G.prior",
                          "mu.0.mu.prior", "mu.0.var.prior",
                          "sigma2.0.alpha.prior", "sigma2.0.beta.prior",
                          "mu.drift.mu.prior", "mu.drift.var.prior",
                          "sigma2.drift.alpha.prior", "sigma2.drift.beta.prior",
                          "Mu.0", "Sigma.0", "cor.0", "sigma1.0", "sigma2.0",
                          "T.idx", "first.idx",
                          "NER.0.up", "GER.0.up",
                          "cd.unique", "yr.unique",
                          "cd.idx.0", "cd.alpha.idx.0",
                          "r"))


MINTS.results <- parLapply(cl, 1:n.chain, sample.block.fxn)

stopCluster(cl)
