library(MASS)
library(tidyverse)
library(rstan)
library(rmutil)
library(bayesplot)
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
  # (1 + (v/parms[["scale"]])^2)^2/(2*(1-(v/parms[["scale"]])^2))* diff(x,delta,lpos)*parms[["scale"]]^2
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
  etjmax <- parms[["jfixedsteps"]]
  etsmax <- parms[["sfixedsteps"]]
  #db <- log(beta*gamma(beta)*sin(pi*beta/2)/pi)
  timestepjmax <- parms[["jumptime"]]
  timestepsmax <- parms[["slowtime"]]
  
  
  # traj[1,] <- c(y,v,log(1+t(v)%*%v)- lpos(y))
  traj[1,] <- c(y,v,1*(v>scale),log(1+t(v)%*%v/scale^2)- lpos(y))
  accept<-1
  iter <- seq(2,N)
  ct <- 0
  # print(iter)
  for (a in iter){
    energy_i <- -log(1+t(v)%*%v/scale^2) + lpos(y)
    y0 <- y
    v0<- v
    if (abs(v) > scale){ #jump
      mm=1 #jump motion mode
      etj <- sample(1:etjmax,1)
      timestepj <- runif(1,0,timestepjmax)
      for(t in seq(1,etj)){
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
      # print("nan")
      prob=0
      energy_f = -Inf
    }
    else{
      # montecarlo
      
      energy_f <-   - log(1+t(v)%*%v/scale^2)  + lpos(ynew)
      prob <- min(1,exp(energy_f-energy_i))
      # print(cbind(y,v,ynew,vnew,prob))
    }
    if (runif(1) < prob){
      traj[a,] <- c(ynew,vnew,mm,-energy_f)
      y<-ynew
      accept = accept+1
    }
    else{
      
      # y<-y
      
      traj[a,] <- c(y,v,mm,-energy_f)
    }
    
    
    v<-rcauchy(length(initial[,2]))*scale
  }
  
  
  
  list(samples=traj,totaliter = ct,acceptance=accept,accprop=accept/N) 
}




nexpt=5
chainlength=200000
nchains=4

set.seed(1029)#first half experiment seed
#set.seed(1021) #2nd half experiment seed
seeds <-sample.int(10000,size=nexpt)
expt=5
set.seed(seeds[expt])
init_y <- runif(nchains,-60,60)

print(init_y)
scaleset <- c(0.05,0.1,0.25,0.5,0.7,1.0,1.25,sqrt(2),1.5,2,5)

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
      # 
      # #KLD
      # h <- hist(samples$y1,plot=F,breaks=500)
      # h$density->cflmc_dens
      # h$mids -> cflmc_mids
      # mmtrue <- sapply(cflmc_mids,true)
      # 
      # KLD(py=mmtrue,px=h$density)$sum.KLD.py.px ->cflmc_kld
      # 
      # #KS
      # ks.test(samples$y1,target_cdf)->cflmcks 
      # 
      # 
      # c(scale=scale,
      #   bias=mean(samples$y1)-5,
      #   KLD=cflmc_kld,
      #   KS=as.numeric(cflmcks$statistic),
      #   acceptance=cauchy_js$accprop,
      #   chain=chain)->info
      # metrics <- bind_rows(metrics,info)
    }
    
  totsamples <- totsamples[-1,]    
    
    # totsamples %>% filter(scale==scale) %>% mutate(time=row_number()) %>% ggplot(aes(x=time,y=y1,color=chain)) + geom_line() +
    #   ggtitle(paste("Trace plot for 0.5 N(0,1) + 0.5 N(10,1) with scale=",scale))->tracep
    # ggsave(tracep,filename=paste("Traceplot_multichain_bimodal_scale=",scale,".png"),height=6,width=8,units="in")
    # 
    # totsamples %>% filter(scale==scale) %>% ggplot(aes(x=y1)) + geom_density() + 
    #   geom_line(data=true_values,aes(x=x,y=density),color="blue")+
    #   geom_density(data=temperedmcmc_values_mc,aes(x=tempered),color="red")+
    #   ggtitle(paste("Density plot for 0.5 N(0,1) + 0.5 N(10,1) with scale=",scale)) -> densityp
    # ggsave(densityp,filename=paste("Densityplot_multichain_bimodal_scale=",scale,".png"),height=6,width=8,units="in")
    
} ->parsamples

parsamples %>% bind_rows -> allsamples
allsamples$expt <- expt
print(dim(allsamples))
saveRDS(allsamples,file=paste("experiment/samplestrimodal_experiment",expt,".rds"))
# saveRDS(metrics,file="metricsbimodal_4chains.rds")
# saveRDS(temperedmcmc_values_mc,file="trimodal_temperedmcmc_lj.rds")
