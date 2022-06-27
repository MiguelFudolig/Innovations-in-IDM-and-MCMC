library(MASS)
library(tidyverse)
library(rstan)
library(rmutil)
library(pracma)
library(LaplacesDemon)
library(foreach)
library(doParallel)
# registerDoParallel(cores=16)


llik<- function(x,lambda=20){
  log(0.3*dgamma(x,shape=3,scale=1) + 0.7*dgamma(x,shape=25,scale=1))
}

#true distribution
true <- function(x){
  exp(llik(x))
}

## TEMPERED MCMC based on R code from Oganisian (2021)
##============================
### ------------------------- Helper Functions  -----------------------------###

vanilla_mh <- function(log_target, theta_curr, prop_sd){
  theta_star <- rnorm(1, theta_curr, prop_sd )
  
  log_accept_prob <- log_target(theta_star) - log_target(theta_curr)
  accept <- log(runif(1)) < log_accept_prob
  if(is.na(accept)) browser()
  if(accept){ 
    return(c(theta_star,accept)) 
  }
  else{ 
    return(c(theta_curr,accept))
    }
}

### tempered distribution for various temperature, temp.
ftemp <- function(x, temp) (true(x))^(1/temp)
log_target <- function(x) log(ftemp(x, temp=1))

### ------------------------- Parallel Tempering  ---------------------------###

temperedmcmc <- function(temperatures=100,cl=chainlength,sd=10,start=0){
  set.seed(10)
  
  iter <- cl
  tempv <- c(1,temperatures)
  
  n_temps <- length(tempv)
  temp_indx <- 1:n_temps
  
  
  theta_shell <- matrix(start, nrow=iter, ncol=n_temps)
  theta_shell_no_switch = matrix(start, nrow=iter, ncol=n_temps)
  accept_shell <- matrix(1, nrow=iter, ncol=n_temps)
  swap_shell <- numeric(length = iter)
  
  for( i in 2:iter){
    
    ## update chains (potentially in parallel )
    for(t in temp_indx){
      log_target <- function(x) log( ftemp(x, temp=tempv[t]) )
      
      prop_sd = ifelse(t==1, 1, sd)
      mh <-  vanilla_mh(log_target = log_target, 
                        theta_curr = theta_shell[i-1, t], 
                        prop_sd = prop_sd)
      theta_shell[i, t] <- mh[1]
      accept_shell[i,t] <- mh[2]
      
    }
    theta_shell_no_switch[i, ] = theta_shell[i, ]
    
    ## propose swap, from swap_idx[1] (chain j) to swap_idx[2] (chain k)
    swap_idx <- sample(temp_indx, 2, replace = T)
    cj <- swap_idx[1]
    ck <- swap_idx[2]
    theta_j <- theta_shell[ i , cj]
    theta_k <- theta_shell[ i , ck]
    
    log_target <- function(x) log(ftemp(x, temp=1))
    f1 <- tempv[cj]*( log_target( theta_j ) - log_target( theta_k ) )
    f2 <- tempv[ck]*( log_target( theta_k ) - log_target( theta_j )  )
    
    accept_prob <- min( c(1, exp(f1 + f2) ) )
    
    if( rbinom(1,1, accept_prob)==1 ){
      
      ## make the swap
      theta_shell[i, cj] <- theta_k
      theta_shell[i, ck] <- theta_j
      
      ## record the swap
      swap_shell[i] = ifelse(cj!=ck, 1, 0)
      
    }
  }
  
  theta_shelldf <- data.frame(tempered = theta_shell[,1],swap=swap_shell)
  acc_rates <- data.frame(cold=accept_shell[,1],hot=accept_shell[,2])
  list(cold=theta_shelldf,acceptance=acc_rates)
}

#=========================

#acceleration is based on having lagrangian as Kinetic + logP
acceleration <-function(x,v,lpos,parms){
  (1 + (v/parms[["scale"]])^2)^2/(2*(1-(v/parms[["scale"]])^2))* fderiv(f=lpos,x=x)*parms[["scale"]]^2
}

leapfrog<- function(eom,y0,v0,lpos,parms,timestep){
  dv <- eom(y0,v0,lpos,parms)
  vnew_2 <- v0+timestep*dv/2
  ynew <- y0+timestep*vnew_2
  vnew <- vnew_2+timestep*eom(ynew,vnew_2,lpos,parms)/2
  c(ynew,vnew)
}

#parms= parameters, eom=equationofmotion(acceleration)
#lpos= log of target distribution, #initial=initial y and v values
#traj=data frame of NA values with 4 columns: y, v,jump, #"energy"

run_cauchy_mm_js <- function(parms,eom,lpos,initial=init,traj=traj){
  #diff timestep for jumps and slow
  scale <- parms[["scale"]]
  N <- parms[["chainlength"]]
  v <- initial[,2]
  y<-initial[,1]
  etjmax <- parms[["jfixedsteps"]]
  etsmax <- parms[["sfixedsteps"]]
  timestepjmax <- parms[["jumptime"]]
  timestepsmax <- parms[["slowtime"]]
  
  
  traj[1,] <- c(y,v,1*(v>scale),log(1+t(v)%*%v/scale^2)- lpos(y))
  accept<-1
  iter <- seq(2,N)
  ct <- 0
  for (a in iter){
    energy_i <- -log(1+t(v)%*%v/scale^2) + lpos(y)
    y0 <- y
    v0<- v
    if (abs(v) > scale){ #jump
      mm=1 #jump motion mode
      etj <- sample(1:etjmax,1)
      timestepj <- runif(1,0,timestepjmax)
      for(t in seq(1,etj)){
        #leapfrog until etj steps
        ynew <- leapfrog(eom,y0,v0,lpos,parms,timestepj)[1]
        vnew <- leapfrog(eom,y0,v0,lpos,parms,timestepj)[2]
        if(t == etj){
          break
        }
        else{
          y0 <- ynew
          v0<-vnew
        }
      }
    }
    else {
      mm=0 #slow motion mode
      ets <- sample(1:etsmax,1)
      timesteps <- runif(1,0,timestepsmax)
      for(t in seq(1,ets)){
        #leapfrog until ets steps
        ynew <- leapfrog(eom,y0,v0,lpos,parms,timesteps)[1]
        vnew <- leapfrog(eom,y0,v0,lpos,parms,timesteps)[2]
        if(t == ets){
          break
        }
        else{
          y0 <- ynew
          v0<-vnew
        }
      }
    }
    
    if((is.nan(lpos(ynew))| sum(is.nan(vnew)) | sum(is.nan(ynew)))){
      prob=0
      energy_f = -Inf
    }
    else{
      # montecarlo/ detailed balance condition
      #t(-v)%*%(-v) = t(v) %*%v 
      energy_f <-   - log(1+t(v)%*%v/scale^2)  + lpos(ynew)
      prob <- min(1,exp(energy_f-energy_i))
    }
    if (runif(1) < prob){
      traj[a,] <- c(ynew,vnew,mm,-energy_f)
      y<-ynew
      accept = accept+1
    }
    else{
      
      traj[a,] <- c(y,v,mm,-energy_f)
    }
    
    
    v<-rcauchy(length(initial[,2]))*scale
  }
  
  
  
  list(samples=traj,totaliter = ct,acceptance=accept,accprop=accept/N) 
}

#getting true values
true_values <- tibble(x=seq(0,50,by=0.001)) %>% mutate(density=sapply(x,true))
true_values %>% mutate(cdf=cumsum(density)*0.001)->true_values


chainlength=150000
nchains=4
set.seed(57)
init_y <- runif(nchains,0,50)
print(init_y)


#tempered MCMC results
tmcmc1 <- temperedmcmc(temperatures=100,cl=chainlength,sd=1,start=init_y[1])
tmcmc2 <- temperedmcmc(temperatures=100,cl=chainlength,sd=1,start=init_y[2])
tmcmc3 <- temperedmcmc(temperatures=100,cl=chainlength,sd=1,start=init_y[3])
tmcmc4 <- temperedmcmc(temperatures=100,cl=chainlength,sd=1,start=init_y[4])


tmcmc_accrates <- bind_rows(tmcmc1[["acceptance"]],
                            tmcmc2[["acceptance"]],
                            tmcmc3[["acceptance"]],
                            tmcmc4[["acceptance"]])
temperedmcmc_values_mc <- bind_rows(tmcmc1[["cold"]],
                                    tmcmc2[["cold"]],
                                    tmcmc3[["cold"]],
                                    tmcmc4[["cold"]])




scaleset <- c(0.05,0.1,0.25,0.5)

foreach(i =1:length(scaleset),.packages = c("MASS","rmutil","tidyverse","pracma","LaplacesDemon"))%dopar%
  {
    totsamples <- data.frame(y1=NA,
                             v1=NA,
                             mode=NA,
                             energy=NA,
                             chain=NA,
                             scale=NA)
    scale = scaleset[i]
    print(paste("scale = ",scale)) #check if parallel works
    for (chain in 1:nchains){
      traj <- data.frame(y1=rep(NA,parms[["chainlength"]]),
                         v1=rep(NA,parms[["chainlength"]]),
                         mode=rep(NA,parms[["chainlength"]]),
                         energy=rep(NA,parms[["chainlength"]]))
      #setting parameters of cauchy flight
      parms <- list(chainlength=chainlength,
                    scale=scale,
                    jfixedsteps=20,
                    sfixedsteps=20,
                    jumptime=0.1,
                    slowtime=0.1)

      #initializing values
      init_v <-rcauchy(1)*parms[["scale"]]
      init <- matrix(c(init_y[chain],init_v),nrow=1)
      
      run_cauchy_mm_js(parms=parms,
                       eom=acceleration,
                       lpos=llik,
                       initial=init,
                       traj=traj)->cauchy_js
      
      
      samples <- cauchy_js$samples
      samples$chain=chain
      samples$scale=scale
      totsamples <- bind_rows(totsamples,samples)

    }
    
  totsamples <- totsamples[-1,]    

} ->parsamples

parsamples %>% bind_rows -> allsamples
print(dim(allsamples))
saveRDS(allsamples,file="gamma/gammasamples.rds")
saveRDS(temperedmcmc_values_mc,file="gammatemperedmcmc_multichain.rds")
saveRDS(tmcmc_accrates,file="gammatemperedmcmc_acceptance.rds")
