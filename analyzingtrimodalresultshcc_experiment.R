library(tidyverse)
library(mcmcse)
readRDS("samplestrimodal_experiment 1 .rds") %>% drop_na() ->exp1 
readRDS("samplestrimodal_experiment 2 .rds")%>% drop_na() ->exp2 
readRDS("samplestrimodal_experiment 3 .rds")%>% drop_na() ->exp3 
readRDS("samplestrimodal_experiment 4 .rds")%>% drop_na() ->exp4 
readRDS("samplestrimodal_experiment 5 .rds") %>% drop_na()->exp5 
readRDS("samplestrimodal_experiment 6 .rds") %>% drop_na() ->exp6 
readRDS("samplestrimodal_experiment 7 .rds") %>% drop_na() ->exp7 
readRDS("samplestrimodal_experiment 8 .rds") %>% drop_na() ->exp8 
readRDS("samplestrimodal_experiment 9 .rds") %>% drop_na() ->exp9 
readRDS("samplestrimodal_experiment 10 .rds") %>% drop_na() ->exp10 
#bind_rows(exp1,exp2,exp3,exp4,exp5)->sampleshcc

# readRDS("samplestrimodal_4chains_parallel_lj.rds") -> sampleshcc
readRDS("trimodal_temperedmcmc_lj.rds") ->tmcmc_samples
# sampleshcc %>%  drop_na() ->sampleshcc

llik<- function(x){
  log(0.2*dnorm(x,sd=2)+0.5*dnorm(x,mean=50)+0.3*dnorm(x,mean=-50))
}

 
#true distribution
true <- function(x){
  exp(llik(x))
}

true_values <- tibble(x=seq(-60,60,by=0.001)) %>% mutate(density=sapply(x,true))
true_values %>% mutate(cdf=cumsum(density)*0.001)->true_values


sink(file="trimodal_expt_stats.txt")
runlength = 200000*4
peak0 <- function(x){
  expe<- experiments[[x]]
  expe%>% filter(between(y1,-30,30))%>% 
    group_by(scale) %>%
    summarise(meany=mean(y1),
              sd = sd(y1),
              n=n()/runlength,
              expt=mean(expt))%>% round(3)
}

experiments <- list(exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8,exp9,exp10)
lapply(1:length(experiments),FUN=peak0)->peak1
peak1 %>% bind_rows() -> peak1

print("stats at x=0")
peak1 %>%
  group_by(scale)%>% 
#  filter(scale==0.25) %>% 
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n))%>% as.data.frame

peak50 <- function(x){
  expe<- experiments[[x]]
  expe%>% filter(y1>30)%>% 
    group_by(scale) %>%
    summarise(meany=mean(y1),sd = sd(y1),n=n()/runlength,expt=mean(expt))%>% round(5)
}

lapply(1:length(experiments),FUN=peak50)->peak2
peak2 %>% bind_rows() -> peak2

print("stats at x=50")
peak2 %>%
  group_by(scale) %>% 
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n))%>% as.data.frame

peakn50 <- function(x){
  expe<- experiments[[x]]
  expe%>% filter(y1< -30)%>% 
    group_by(scale) %>%
    summarise(meany=mean(y1),sd = sd(y1),n=n()/runlength,expt=mean(expt))%>% round(5)
}

lapply(1:length(experiments),FUN=peakn50)->peak3
peak3 %>% bind_rows() -> peak3

print("stats at x=-50")
peak3 %>%
  group_by(scale) %>% 
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n))%>% as.data.frame


peakoverall_scale <-function(x){
  expe<- experiments[[x]]
  expe%>% group_by(scale)%>%
    mutate(displacement=c(1,diff(y1)),accept=c(1,1*(abs(diff(y1)) > 1e-7)))%>%
    summarise(meany=mean(y1),
              sd = sd(y1),
              n=n()/runlength,
              acceptancerate=sum(accept)/n(),
              expt=mean(expt))
}

lapply(1:length(experiments),FUN=peakoverall_scale)->peakovscale
peakovscale %>% bind_rows() -> peakovscale

print("stats per scale")
peakovscale %>%
  group_by(scale) %>% 
  summarise(mean=mean(meany),
            bias=mean(meany)-10,
            se=sd(meany),
            meansd=mean(sd),
            sdsd=sd(sd),
            meanaccratio=mean(acceptancerate),
            sdar=sd(acceptancerate))%>%
  mutate(t=(mean-10)/se,p=2*pt(abs(t),df=9,lower.tail=F)) %>%
  as.data.frame()

print("stats all scales combined")
peakoverall <-function(x){
  expe<- experiments[[x]]
  expe%>% 
    summarise(meany=mean(y1),sd = sd(y1),n=n()/runlength,expt=mean(expt))%>% round(3)
}

lapply(1:length(experiments),FUN=peakoverall)->peako
peako %>% bind_rows() -> peako

peako  %>% 
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n))

sink(file=NULL)

#Experiment 1 test==============
##overall mean
exp1%>% mutate(accept=c(1,1*(abs(diff(y1)) > 1e-7)))%>%
  summarise(mean=mean(y1),sd = sd(y1),n=n(),acceptance=sum(accept)/n())%>% round(3)
#separated by scale
exp1 %>% mutate(accept=c(1,1*(abs(diff(y1)) > 1e-7)))%>%
  group_by(scale) %>%
  summarise(mean=mean(y1),bias = mean(y1)-10,sd = sd(y1),acceptance=sum(accept)/n())%>% round(3)



#peak at 0
exp1 %>% filter(between(y1,-30,30))%>% 
  group_by(scale) %>%
  summarise(mean=mean(y1),sd = sd(y1),n=n()/runlength)%>% round(3)
##overall

exp1 %>% filter(between(y1,-30,30))%>% 
  summarise(mean=mean(y1),sd = sd(y1),n=n()/nrow(sampleshcc))%>% round(3)

#traceplot

exp1 %>% filter(scale%in% c(0.05,1,5)) %>% group_by(scale) %>%
  mutate(time=row_number()) %>%
  ggplot(aes(x=time,y=y1)) + geom_line()+
  ylab(expression(theta)) + xlab("Time")+
  facet_grid(vars(as.factor(scale)))->traceplot
ggsave(traceplot,filename="traceplot trimodal.png",
       width=8,height = 8,units="in",dpi=300)

#==================



#convergence check for exp 1 with 0.05

sink("trimodal diagnostics exp1.txt")
#at zero

nchains=4
grstat0 <- function(s)
  {exp1 %>% filter(between(y1,-30,30))%>% 
    filter(scale==s) ->zeropeak_005
  
  zeropeak_005 %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(zeropeak_005$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

ess0 <- function(s){
  exp1 %>% filter(between(y1,-30,30))%>% 
    filter(scale==s) ->zeropeak_005
  ess(zeropeak_005$y1)
}
print("Gelman-Rubin Statistic at 0")
exp1 %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat0))
print("ess at 0")
exp1 %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=ess0))
#must be above 1000


#at 50

nchains=4
grstat50 <- function(s)
{exp1 %>% filter(y1>30)%>% 
    filter(scale==s) ->zeropeak_005
  
  zeropeak_005 %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(zeropeak_005$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

ess50 <- function(s){
  exp1 %>% filter(y1>30)%>% 
    filter(scale==s) ->zeropeak_005
  ess(zeropeak_005$y1)
}
print("Gelman-Rubin Statistic at 50")
exp1 %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat50))
print("ess at 50")
exp1 %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=ess50))
#must be above 1000



#at -50

nchains=4
grstatn50 <- function(s)
{exp1 %>% filter(y1< -30)%>% 
    filter(scale==s) ->zeropeak_005
  
  zeropeak_005 %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(zeropeak_005$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

essn50 <- function(s){
  exp1 %>% filter(y1< -30)%>% 
    filter(scale==s) ->zeropeak_005
  ess(zeropeak_005$y1)
}

print("Gelman-Rubin Statistic at -50")
exp1 %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstatn50))

print("ess at -50")
exp1 %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=essn50)) #must be above 1000
sink(file=NULL)
