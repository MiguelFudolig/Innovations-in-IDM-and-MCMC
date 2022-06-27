library(MASS)
library(tidyverse)
library(rstan)
library(rmutil)
library(pracma)
library(LaplacesDemon)
library(foreach)
library(doParallel)
registerDoParallel(cores=32)

llik<- function(x){
  log(0.2*dnorm(x,sd=2)+0.5*dnorm(x,mean=50)+0.3*dnorm(x,mean=-50))
}

#true distribution
true <- function(x){
  exp(llik(x))
}

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

run_cauchy_mm_js <- function(parms,eom,lpos,initial=init,traj=traj){
  #diff timestep for jumps and slow
  scale <- parms[["scale"]]
  N <- parms[["chainlength"]]
  v <- initial[,2]
  y<-initial[,1]
  #establish max fixed steps for jump steps and slow/crawl steps
  etjmax <- parms[["jfixedsteps"]]
  etsmax <- parms[["sfixedsteps"]]
  #establish time step 
  timestepjmax <- parms[["jumptime"]]
  timestepsmax <- parms[["slowtime"]]
  
  #setting trajectory
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
        #leapfrog integrator
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
      # montecarlo/ metropolis-hastings for detailed balance
      energy_f <-   - log(1+t(v)%*%v/scale^2)  + lpos(ynew)
      prob <- min(1,exp(energy_f-energy_i))
    }
    if (runif(1) < prob){
      traj[a,] <- c(ynew,vnew,mm,-energy_f)
      y<-ynew
      accept = accept+1
    }
    else{
      
      # y<-y, keep y
      
      traj[a,] <- c(y,v,mm,-energy_f)
    }
    
    #resample velocitu
    v<-rcauchy(length(initial[,2]))*scale
  }
  list(samples=traj,totaliter = ct,acceptance=accept,accprop=accept/N) 
}




nexpt=5
chainlength=200000
nchains=4

set.seed(1029) #first half experiment seed
#set.seed(1021) #2nd half experiment seed
seeds <-sample.int(10000,size=nexpt)
expt=1
set.seed(seeds[expt])
init_y <- runif(nchains,-60,60)

print(init_y)
scaleset <- c(0.05,0.1,0.25,0.5,0.7,1.0,1.25,sqrt(2),1.5,2,5)
#run each scale in parallel
foreach(i =1:length(scaleset),.packages = c("MASS","rmutil","tidyverse","pracma","LaplacesDemon"))%dopar%
  {
    totsamples <- data.frame(y1=NA,
                             v1=NA,
                             mode=NA,
                             energy=NA,
                             chain=NA,
                             scale=NA)
    scale = scaleset[i]
    print(paste("scale = ",scale))
    for (chain in 1:nchains){
      traj <- data.frame(y1=rep(NA,parms[["chainlength"]]),
                         v1=rep(NA,parms[["chainlength"]]),
                         mode=rep(NA,parms[["chainlength"]]),
                         energy=rep(NA,parms[["chainlength"]]))
      
      parms <- list(chainlength=chainlength,
                    scale=scale,
                    jfixedsteps=40,
                    sfixedsteps=20,
                    jumptime=0.1,
                    slowtime=0.1)


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
allsamples$expt <- expt
print(dim(allsamples))
saveRDS(allsamples,file=paste("experiment/samplestrimodal_experiment",expt,".rds"))
