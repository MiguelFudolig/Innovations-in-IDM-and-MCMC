library(tidyverse)
library(mcmcse)
library(LaplacesDemon)
# readRDS("metricsbimodal_4chains.rds")[-1,] -> metrics
readRDS("gammasamples.rds") %>% drop_na() -> sampleshcc
# readRDS("gammasamples_vnew2.rds") %>% drop_na() -> sampleshcc
# readRDS("gammasamples_vnew2.rds") %>% drop_na() -> sampleshcc

# readRDS("gammatemperedmcmc.rds") ->tmcmc_samples
readRDS("gammatemperedmcmc_multichain.rds") ->tmcmc_samples
readRDS("gammatemperedmcmc_acceptance.rds") ->tmcmc_acceptance


llik<- function(x,lambda=20){
  log(0.3*dgamma(x,shape=3,scale=1) + 0.7*dgamma(x,shape=25,scale=1))
}

#true distribution
true <- function(x){
  exp(llik(x))
}

true_values <- tibble(x=seq(0,50,by=0.01)) %>% mutate(density=sapply(x,true))
true_values %>% mutate(cdf=cumsum(density)*0.001)->true_values

truemean <- 0.3*3+0.7*25



sampleshcc %>%group_by(scale) %>%mutate(accept=c(1,1*(abs(diff(y1)) > 1e-7))) %>%ungroup->sampleshcc
tmcmc_samples %>%mutate(chain=rep(1:4,each=n()/4))->tmcmc_samples
# 
# set.seed(11)
# ensemble1 <- sample(unique(sampleshcc$scale),size=4)
# ensemble2 <- sample(unique(sampleshcc$scale),size=4)
# ensemble3 <- sample(unique(sampleshcc$scale),size=4)
# 

sampleshcc %>% group_by(scale) %>% mutate(time=row_number()) %>%
  ggplot(aes(x=time,y=y1)) + geom_line() +
  facet_wrap(~as.factor(scale),ncol=2)+
  ylim(c(0,50))+
  #scale_x_continuous(breaks=seq(0,8e5,by=5e4),name="Time") +
  ylab(expression(theta))+xlab("Time")+ geom_hline(yintercept=25,color="blue")+
  geom_hline(yintercept=3,color="red")

ggsave(filename="allscalestrace.png",
       width=8,
       height=8,
       units="in") 

sampleshcc   %>%
  #filter(scale%in%c(0.50)) %>%
    ggplot(aes(x=y1)) + 
  geom_line(data = true_values,
            aes(x=x,y=density),
            color="blue",
            size=1,
            linetype=1)+
    geom_density(linetype=5,size=1)+
  geom_density(data=tmcmc_samples,
               aes(x=tempered),
               color="red",
               linetype=3,size=1)+
  xlab(expression(theta)) + ylab("Density") +
  scale_linetype_discrete(name="Scale",
                          limits=as.factor(c(1,3,5)),
                          labels=c("True","TMCMC","CFLMC"))+
  facet_wrap(~as.factor(scale),ncol=2)+
  xlim(0,60)
  

ggsave(filename="allscales2.png",
       width=8,
       height=8,
       units="in") 
  
sampleshcc %>% group_by(scale,chain) %>% summarise(mean=mean(y1)) %>% round(3)->chainmeans
View(chainmeans)

#jumps

sampleshcc %>% group_by(scale) %>% filter(accept==1)%>% summarise(totjumps=sum(mode)/n())
sampleshcc %>% filter(scale==5) %>% head(20) %>% View

sink(file="gamma peaks.txt")

#overall means
print("cflmc overall means")

sampleshcc %>% 
  group_by(scale) %>% 
  summarise(mean=mean(y1),bias=mean(y1)-truemean,acceptance=sum(accept)/n()) %>% round(5)%>% as.data.frame

# sampleshcc %>% 
# #  group_by(scale) %>% 
#   summarise(mean=mean(y1),bias=mean(y1)-truemean,acceptance=sum(accept)/n()) %>% round(5)
# 

#peak at 2

print("peak at 2")
sampleshcc %>% filter(y1 <10)%>% 
  group_by(scale) %>%
  summarise(mean=mean(y1),meanbias=mean(y1)-3,sd = sd(y1),n=n()/(150000*4)) %>% round(5) %>% as.data.frame()

#peak at 24

print("peak at 24")
sampleshcc %>% filter(y1 >10)%>% 
  group_by(scale) %>%
  summarise(mean=mean(y1),meanbias=mean(y1)-25,sd = sd(y1),n=n()/(150000*4)) %>% round(5) %>% as.data.frame()


#temperedmcmc
print("tempered")
tmcmc_samples %>% summarise(mean=mean(tempered),sd=sd(tempered),bias=mean(tempered)-truemean)

print("tempered at 2")
tmcmc_samples %>% filter(tempered<10) %>% summarise(mean=mean(tempered),bias=mean(tempered)-3,sd=sd(tempered),n=n()/(150000*4))
print("tempered at 24")
tmcmc_samples %>% filter(tempered>10) %>% summarise(mean=mean(tempered),bias=mean(tempered)-25,sd=sd(tempered),n=n()/(150000*4))

tmcmcacc <- tmcmc_acceptance %>% summarise(coldrate = sum(cold)/n(),hotrate=sum(hot)/n())
tmcmcacc
sink(file=NULL)


#===========looking at 0.05 (best candidate) convergence ==========
sink("gamma diagnostics scale.txt")
#at zero


print("overall diagnostics")
nchains=4
grstat <- function(s)
{sampleshcc %>% filter(scale==s) ->peak
  
  peak %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(peak$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

essov <- function(s){
  sampleshcc %>%  filter(scale==s) ->peak
  ess(peak$y1)
}

print("Gelman-Rubin Statistic")
sampleshcc %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat))
print("ess ")
sampleshcc %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=essov))

print("at 2")
nchains=4
grstat2 <- function(s)
{sampleshcc %>% filter(y1 <10)%>% filter(scale==s) ->peak
  
  peak %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(peak$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

ess2 <- function(s){
  sampleshcc %>% filter(y1 <10)%>% filter(scale==s) ->peak
  ess(peak$y1)
}
print("Gelman-Rubin Statistic")
sampleshcc %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat2))
print("ess ")
sampleshcc %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=ess2))

print("at 24")

grstat24 <- function(s)
{sampleshcc %>% filter(y1 >10)%>% filter(scale==s) ->peak
  
  peak %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(peak$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

ess24 <- function(s){
  sampleshcc %>% filter(y1 >10)%>% filter(scale==s) ->peak
  ess(peak$y1)
}
print("Gelman-Rubin Statistic")
sampleshcc %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat24))
print("ess ")
sampleshcc %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=ess24))

print("tempered mcmc convergence")
rhat <- function(peak)
{ peak %>% group_by(chain) %>% summarize(mean=mean(tempered),n=n(),var=sd(tempered)^2) ->means
  mean(peak$tempered)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}
print("at 24")
tmcmc_samples %>% filter(tempered>10) %>%summarise(ess=ess(tempered))
print("rhat")
tmcmc_samples %>% filter(tempered>10) %>%rhat
print("at 2")
tmcmc_samples %>% filter(tempered<10) %>%summarise(ess=ess(tempered))
print("rhat")
tmcmc_samples %>% filter(tempered<10) %>%rhat
#=====================Kullback Leibler Divergence==============
print("KLD")
comp_kld <- function(s){
best <- sampleshcc %>% filter(scale==s)
cflmc_hist <- hist(best$y1,plot=F,breaks=500)
cflmc_dens <- cflmc_hist$density
cflmc_mids <- cflmc_hist$mids


mmtrue <- sapply(cflmc_mids,true)

KLD(py=mmtrue,px=cflmc_dens)$sum.KLD.py.px 
}

sampleshcc %>% summarise(cflmc_scales=unique(scale),kld=sapply(cflmc_scales,FUN=comp_kld))
#tempered mcmc

tmcmc_hist <- hist(tmcmc_samples$tempered,plot=F,breaks=500)
tmcmc_dens <- tmcmc_hist$density
tmcmc_mids <- tmcmc_hist$mids


mmtrue <- sapply(tmcmc_mids,true)
KLD(py=mmtrue,px=tmcmc_dens)$sum.KLD.py.px ->tmcmc_kld
c(tmcmc=tmcmc_kld)


#==========kolgomorov smirnov statistic============
print("KS")

target_cdf <- function(x){
  0.3*pgamma(x,shape=3,scale=1) + 0.7*pgamma(x,shape=25,scale=1)
}

comp_ks <- function(s){
  best <- sampleshcc %>% filter(scale==s)
  ks.test(best$y1,target_cdf)$statistic
}
  
sampleshcc %>% summarise(cflmc_scales=unique(scale),ks=sapply(cflmc_scales,FUN=comp_ks))

ks.test(tmcmc_samples$tempered,target_cdf)$statistic->tmcmcks
tmcmcks

sink(file=NULL)
