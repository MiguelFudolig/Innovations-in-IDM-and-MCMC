library(tidyverse)
library(LaplacesDemon)
library(mcmcse)
readRDS("samplesbimodal_experiment 1 .rds") %>% drop_na() ->exp1 
readRDS("samplesbimodal_experiment 2 .rds")%>% drop_na() ->exp2 
readRDS("samplesbimodal_experiment 3 .rds")%>% drop_na() ->exp3 
readRDS("samplesbimodal_experiment 4 .rds")%>% drop_na() ->exp4 
readRDS("samplesbimodal_experiment 5 .rds") %>% drop_na()->exp5 
readRDS("samplesbimodal_experiment 6 .rds") %>% drop_na() ->exp6 
readRDS("samplesbimodal_experiment 7 .rds") %>% drop_na() ->exp7 
readRDS("samplesbimodal_experiment 8 .rds") %>% drop_na() ->exp8 
readRDS("samplesbimodal_experiment 9 .rds") %>% drop_na() ->exp9 
readRDS("samplesbimodal_experiment 10 .rds") %>% drop_na() ->exp10 
#bind_rows(exp1,exp2,exp3,exp4,exp5)->sampleshcc

# readRDS("samplestrimodal_4chains_parallel_lj.rds") -> sampleshcc
readRDS("bimodal_temperedmcmc.rds") ->tmcmc_samples
# sampleshcc %>%  drop_na() ->sampleshcc

llik<- function(x){
  log(0.5*dnorm(x,sd=1)+0.5*dnorm(x,mean=10))
}

 
#true distribution
true <- function(x){
  exp(llik(x))
}

true_values <- tibble(x=seq(-5,15,by=0.001)) %>% mutate(density=sapply(x,true))
true_values %>% mutate(cdf=cumsum(density)*0.001)->true_values

options(tibble.pillar.sigfig=7)
sink(file="bimodal_expt_stats.txt")
runlength = 150000*4
peak0 <- function(x){
  expe<- experiments[[x]]
  expe%>% filter(y1<5)%>% 
    group_by(scale) %>%
    summarise(meany=mean(y1),sd = sd(y1),n=n()/runlength,expt=mean(expt))%>% round(3)
}

experiments <- list(exp1,exp2,exp3,exp4,exp5,exp6,exp7,exp8,exp9,exp10)
lapply(1:length(experiments),FUN=peak0)->peak1
peak1 %>% bind_rows() -> peak1

print("stats at x=0")
peak1 %>%
  group_by(scale)%>% 
#  filter(scale==0.25) %>% 
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n))%>%
  as.data.frame()

peak50 <- function(x){
  expe<- experiments[[x]]
  expe%>% filter(y1>5)%>% 
    group_by(scale) %>%
    summarise(meany=mean(y1),sd = sd(y1),n=n()/runlength,expt=mean(expt))%>% round(5)
}

lapply(1:length(experiments),FUN=peak50)->peak2
peak2 %>% bind_rows() -> peak2

print("stats at x=10")
peak2 %>%
  group_by(scale) %>% 
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n)) %>%
  as.data.frame

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
            bias=mean(meany)-5,
            se=sd(meany),
            meansd=mean(sd),
            sdsd=sd(sd),
            meanaccratio=mean(acceptancerate),
            sdar=sd(acceptancerate))%>%
  mutate(t=(mean-5)/se,p=2*pt(abs(t),df=9,lower.tail=F)) %>%
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
  summarise(mean=mean(meany),se=sd(meany),meansd=mean(sd),sdsd=sd(sd),meann=mean(n),sdn=sd(n))%>%
  as.data.frame()

print("tempered mean stats")

#temperedmcmc
print("overall")
tmcmc_samples %>% summarise(mean=mean(tempered),sd=sd(tempered))

print("tempered 0")
#peak at 0
tmcmc_samples %>% filter(tempered < 5)%>% 
  summarise(mean=mean(tempered),sd = sd(tempered),n=n()/(150000*4))

print("tempered50")
#peak at 50
tmcmc_samples %>% filter(tempered > 5)%>% 
  summarise(mean=mean(tempered),sd = sd(tempered),n=n()/(150000*4))
sink(file=NULL)

#Experiment 1 test==============
##overall mean
exp1%>% mutate(accept=c(1,1*(abs(diff(y1)) > 1e-7)))%>%
  summarise(mean=mean(y1),sd = sd(y1),n=n(),acceptance=sum(accept)/n())%>% round(3)
#separated by scale
exp1 %>% mutate(accept=c(1,1*(abs(diff(y1)) > 1e-7)))%>%
  group_by(scale) %>%
  summarise(mean=mean(y1),bias = mean(y1)-5,sd = sd(y1),acceptance=sum(accept)/n())%>% round(3)




runlength<- 150000*4

#peak at 0
exp1 %>% 
  filter(y1 < 5)%>% 
  group_by(scale) %>%
  summarise(mean=mean(y1),sd = sd(y1),n=n()/runlength)%>% round(3)

#peak at 50

exp1 %>% filter(y1 > 5)%>% 
  group_by(scale) %>%
  summarise(mean=mean(y1),sd = sd(y1),n=n()/runlength)%>% round(3)

exp1 %>% filter(scale%in% c(0.1))-> ex
hist(ex$y1,plot=F,breaks=25)->gg
ggez <- data.frame(mids=gg$mids,density=gg$density)


exp1 %>% filter(scale%in% c(0.05,1,5)) %>%
  ggplot(aes(x=y1)) + 
  # geom_histogram(aes(y=stat(..count.. / sum(..count..)),group=factor(scale)),
  #                bins=25,
  #                # binwidth=0.75,
  #                fill="white",
  #                size=1.1)+
  # geom_histogram(mapping=aes(y=stat(density)),
  #                bins=25,
  #                # binwidth=0.8,
  #                fill="white",
  #                size=1.1)+
  geom_density(size=0.8,linetype=5)+
  # geom_histogram(aes(y=..density..),binwidth=0.02,alpha=0.7)+
 # geom_point(data=ggez,aes(x=mids,y=density))+
  geom_line(data=true_values,
            aes(x=x,y=density),
            color="blue", 
            size=1.1)+
  # geom_density(data=tmcmc_samples,aes(x=tempered),
  #              color="red",
  #              linetype=3,
  #              size=0.8) + 
  xlab(expression(theta)) + ylab("Density")+
  xlim(c(-5,15)) + 
  # scale_linetype_discrete(name="Scale",
  #                         limits=factor(c(1,3,5)),
  #                         labels=c("True","TMCMC","CFLMC"))+
  scale_linetype_discrete(name="Scale",
                          limits=factor(c(1,5)),
                          labels=c("True","CFLMC"))+
  facet_grid(vars(as.factor(scale)))->pdfplot
pdfplot
ggsave(pdfplot,file="bimodal_dfplotexample_notmcmc.png",width=8,height=6,units="in")


exp1 %>% filter(scale%in% c(0.05,1,5)) %>%
  ggplot(aes(x=y1,color=as.factor(scale))) + 
  stat_ecdf(geom="step")+
  stat_ecdf(data=tmcmc_samples,
            mapping=aes(x=tempered),
            color="red",
            linetype=2)+
  geom_line(data=true_values,aes(x=x,y=cdf),color="blue", size=0.8,linetype=5)+
  # xlab(expression(theta)) + ylab("Density")+
  xlim(c(-5,15)) + ylim(c(0.3,0.7))+
  scale_color_discrete(name="Scale")->cdfplot
cdfplot
ggsave(cdfplot,file="bimodal_cdfplotexample.png",width=8,height=6,units="in")

#traceplots for exp 1


#==================


exp1 %>% filter(scale%in% c(0.05,1,5)) %>% group_by(scale) %>%
  mutate(time=row_number()) %>%
  ggplot(aes(x=time,y=y1)) + geom_line()+
  ylab(expression(theta)) + xlab("Time")+
  facet_grid(vars(as.factor(scale)))->traceplot
ggsave(traceplot,filename="biexptrace.png", width=8, height=8, units="in")

#convergence check for exp 1 with 0.05

sink("bimodal diagnostics exp1.txt")
#at zero

nchains=4
grstat0 <- function(s)
{exp1 %>% filter(y1<5)%>% 
    filter(scale==s) ->zeropeak_005
  
  zeropeak_005 %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(zeropeak_005$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

ess0 <- function(s){
  exp1 %>% filter(y1<5)%>% 
    filter(scale==s) ->zeropeak_005
  ess(zeropeak_005$y1)
}
print("Gelman-Rubin Statistic at 0")
exp1 %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat0))
print("ess at 0")
exp1 %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=ess0))

#at ten

nchains=4
grstat10 <- function(s)
{exp1 %>% filter(y1>5)%>% 
    filter(scale==s) ->zeropeak_005
  
  zeropeak_005 %>% group_by(chain) %>% summarize(mean=mean(y1),n=n(),var=sd(y1)^2) ->means
  mean(zeropeak_005$y1)-> grandmean
  means <- means %>% mutate(B=(mean-grandmean)^2)
  means %>% mutate(W=sum(var)/nchains,sumW=sum(var*((n-1)/n))/nchains,sumB=sum(B))->means
  
  Rhat=sqrt((unique(means$sumW)+unique(means$sumB))/unique(means$W))
  Rhat
}

ess10 <- function(s){
  exp1 %>% filter(y1>5)%>% 
    filter(scale==s) ->zeropeak_005
  ess(zeropeak_005$y1)
}
print("Gelman-Rubin Statistic at 10")
exp1 %>%summarise(scale=unique(scale),Rhat=sapply(unique(scale),FUN=grstat10))
print("ess at 10")
exp1 %>%summarise(scale=unique(scale),ess=sapply(unique(scale),FUN=ess10))
sink(file=NULL)
