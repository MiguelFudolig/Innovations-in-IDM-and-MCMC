library(deSolve)
library(dplyr)
library(rstan)
library(outbreaks)
library(bayesplot)
library(tidyverse)

#setwd("C:/Users/migsf/Google Drive/Migs Files/UNL life/Research/HCC")
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

knitr::opts_chunk$set(cache = TRUE, echo = TRUE, message = FALSE, warning = FALSE)



#flu virus data

read_csv("2008data.csv") -> flu2008
head(flu2008)
flu2008$others <- flu2008$a_h1 + flu2008$a_h3 + flu2008$a_subtyping_not_performed + flu2008$a_unable_to_subtype + flu2008$b + flu2008$h3n2v
#where nebraska is
# flu2008 <- flu2008[flu2008$region=="Region 7",(names(flu2008) %in% c("week", "total_specimens","a_2009_h1n1","others"))]
flu2008 <- flu2008 %>% filter(region=="Region 7") %>% select(week, total_specimens,a_2009_h1n1,others) %>% filter(a_2009_h1n1> 0)
flu2008$time <- 1:nrow(flu2008) 

mod1_stat <- '
functions {
  real[] SIR(real t,  // time
  real[] y,           // system state {susceptible,infected,recovered}
  real[] theta,       // parameters
  real[] x_r,
  int[] x_i) {
  
  real dy_dt[4];
  
  dy_dt[1] = - theta[1] * y[1] * y[2];
  dy_dt[2] = theta[1] * y[1] * y[2] - theta[2] * y[2];
  dy_dt[3] = theta[2] * y[2];
  dy_dt[4] = 0;
  return dy_dt;
  }
  
  }
  data {
  int<lower = 1> n_obs;       // number of days observed
  int<lower = 1> n_theta;     // number of model parameters
  int<lower = 1> n_difeq;     // number of differential equations
  int<lower = 1> n_pop;       // population 
  int y[n_obs];           // data, total number of infected individuals each day
  real t0;                // initial time point (zero)
  real ts[n_obs];         // time points observed 
  real p;                // vaccination rate
  }
  
  transformed data {
  real x_r[0];
  int x_i[0];
  }
  
  parameters {
  real<lower = 0> theta[n_theta]; // model parameters {beta,gamma}
  real<lower = 0> rho;   // underreporting rate
  real<lower = 0, upper = 1-p> S0;
  }
  
  transformed parameters{
  real y_hat[n_obs, n_difeq]; // solution from the ODE solver
  real y_init[n_difeq];     // initial conditions for both fractions of S and I
  
  y_init[1] = S0;
  y_init[2] = 1 - S0 - p;
  y_init[3] = 0;
  y_init[4] = p;
  y_hat = integrate_ode_rk45(SIR, y_init, t0, ts, theta, x_r, x_i);
  
  }
  
  model {
  real lambda[n_obs];      //poisson parameter
  
  //priors
  theta[1] ~ lognormal(0,0.5);
  theta[2] ~ gamma(0.04,0.02);  //Assume mean infectious period = 3 days/0.42857 weeks 
  S0 ~ beta(0.5, 0.5);
  rho ~ normal(79,20);
  

  //likelihood
  for (i in 1:n_obs){
  lambda[i] = y_hat[i,2]*n_pop/rho;
  }
  y ~ poisson(lambda);
  }
  
  generated quantities {
  real R_0;      // Basic reproduction number
  R_0 = theta[1]/theta[2]*(1-p);
  }
'

m1 <- stan_model(model_code = mod1_stat)


## implementing HMC
pop <- floor((1.8E6 + 3E6 + 6E6 + 2.8E6))
flu_data = list(n_obs = max(flu2008$time),
                n_theta = 2,
                n_difeq = 4,
                n_pop = pop,
                y = flu2008$others,
                t0 = 0,
                ts = flu2008$time,
	    p=0.5
                )


n_chains=4
n_warmups=500
n_iter=100500
n_thin=10
set.seed(1234)
# Set initial values:
ini_1 = function(){
  list(theta=c(runif(1,0.3,0.5),runif(1,0.3,0.5)), 
       S0=runif(1,0.49,0.499),rho=runif(1,77,85))  
}
#ini_1 = function(){
#	list(theta=c(2.6,2),S0=0.999)
#}
parameters = c("y_hat", "y_init", "theta",  "R_0", "S0", "rho")

time.start_nuts1 <- Sys.time()
nuts_fit_1 = sampling(m1, data = flu_data, init = ini_1, chains = n_chains, warmup = n_warmups, iter = n_iter, thin=n_thin, seed=13219,
					sample_file="origstrain_samples.csv",diagnostic_file="origstrain_diagnostics.csv")
time.end_nuts1 <- Sys.time()
duration_nuts1<- time.end_nuts1 - time.start_nuts1
saveRDS(nuts_fit_1, "orig_singlestrain_RRfit_gp4chains.RDS")

nuts_fit_1_summary <- summary(nuts_fit_1, pars = c("lp__", "theta[1]", "theta[2]", "y_init[1]", "R_0","rho"))$summary
print(nuts_fit_1_summary,scientific=FALSE,digits=6)
posts_1 <-  rstan::extract(nuts_fit_1)

mod1_diagnostics <-rstan::get_sampler_params(nuts_fit_1)

# Check for divergent transitions
rstan::check_divergences(nuts_fit_1)

pdf("trace and pairs plots orig 1strain randomrho 4chains.pdf")
posterior_1 <- as.array(nuts_fit_1)
color_scheme_set("viridis")
# Markov chain traceplots
mcmc_trace(posterior_1, pars="lp__")
mcmc_trace(posterior_1, pars=c("theta[1]", "theta[2]", "y_init[1]"))
mcmc_trace(posterior_1, pars=c("R_0","rho"))

# Univariate and bivariate marginal posterior distributions
pairs(nuts_fit_1, pars = c("theta[1]", "theta[2]", "y_init[1]","R_0","rho"),
	 labels = c("beta", "gamma", "s(0)","Rep. Number", "URF","Vax rate"), 
      cex.labels=1.5, font.labels=9, condition = "accept_stat__")  

# Kernel density estimates of each Markov chain separately, overlaid
mcmc_dens_overlay(posterior_1, pars=c("theta[1]", "theta[2]", "y_init[1]","R_0","rho"))
dev.off()
#Central posterior uncertainty intervals
mcmc_intervals(posterior_1,pars = c("theta[1]", "theta[2]", "y_init[1]","R_0","rho"))

summary(do.call(rbind, mod1_diagnostics), digits = 2)
