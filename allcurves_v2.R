library(tidyverse)
library(deSolve)

flu2008 <- read_csv("2008data.csv")
flu2008%>% mutate(others=a_h1 +a_h3 + a_subtyping_not_performed + a_unable_to_subtype + b + h3n2v)->flu2008
#where nebraska is
# flu2008 <- flu2008[flu2008$region=="Region 7",(names(flu2008) %in% c("week", "total_specimens","a_2009_h1n1","others"))]
fpe2008 <- flu2008 %>% filter(region=="Region 7") %>% select(week, total_specimens,a_2009_h1n1,others) %>% filter(others > 0)
flu2008 <- flu2008 %>% filter(region=="Region 7") %>% select(week, total_specimens,a_2009_h1n1,others) %>% filter(a_2009_h1n1 > 0)
fpe2008 <- fpe2008 %>% anti_join(flu2008,by="week")
fpe2008$time <- 1:nrow(fpe2008)
flu2008$time <- 1:nrow(flu2008)
time <- flu2008$time


ms <- lapply(1:4,function(x) read.csv(paste("multistrain_samples_",x,".csv",sep=""),skip=26))%>%
        bind_rows %>% drop_na
os <- lapply(1:4,function(x) read.csv(paste("origstrain_samples_",x,".csv",sep=""),skip=26))%>%
  bind_rows %>% drop_na

es <-  lapply(1:4,function(x) read.csv(paste("singlestrainh1n1_samples_",x,".csv",sep=""),skip=26))%>%
  bind_rows %>% drop_na


ms %>%select("theta.1",
             "theta.2",
             "theta.3",
             "theta.4",
             "S0",
             "q",
             "rho")%>%
      drop_na()-> ms

os %>%select("theta.1",
             "theta.2",
             "S0",
             "rho")%>%
  drop_na()-> os

es %>%select("theta.1",
             "theta.2",
             "S0",
             "rho")%>%
  drop_na()-> es

ms_SIR <- function(t,y,parms){
  with(as.list(c(y,parms)),{
    ds  <- (-1*beta*s*i - beta2*s*i_n)
    di <- beta*s*i - gamma*i
    dr <- gamma*i - beta2*r*i_n
    dv <- (-1*beta2*v*i_n)
    di_n <- beta2*(s+v+r)*i_n - gamma2*i_n
    df <- gamma2*i_n #final removal 
    
    list(c(ds,di,dr,dv,di_n,df))
  }
  )
}

os_SIR <- function(t,y,parms){
  with(as.list(c(y,parms)),{
  ds  <- (-1*beta*s*i)
  di <- beta*s*i - gamma*i
  dr <- gamma*i
  dv <- 0
  list(c(ds,di,dr,dv))
  }
  )
}

es_SIR <- function(t,y,parms){
  with(as.list(c(y,parms)),{
  ds  <- (-1*beta*s*i)
  di <- beta*s*i - gamma*i
  dr <- gamma*i
  list(c(ds,di,dr))
  }
  )
}

time <- flu2008$time
p <- 0.5
N=13.6E6


infections <- function(x){
  s0 = x["S0"] %>% as.numeric
  emergerate=x["q"] %>% as.numeric
  t1 <- x["theta.1"] %>% as.numeric
  t2 <- x["theta.2"] %>% as.numeric
  t3 <- x["theta.3"] %>% as.numeric
  t4 <- x["theta.4"] %>% as.numeric
  rho = x["rho"] %>% as.numeric
  initial <- c(s = s0,
             i=(1-s0-p)*(1-emergerate),
             r = 0,
             v =p,
             i_n = (1-s0-p)*emergerate,
             f=0)

  parameters <- c(beta=t1,
                gamma=t2,
                beta2=t3,
                gamma2=t4)
  # print(initial)
  # print(parameters)
  ode(y=initial,times=time,func=ms_SIR,parms=parameters,method = "ode45")-> surveillance
  surveillance %>% as.data.frame %>%
    select(c("i","i_n")) %>%
    mutate(origcases=N/rho*i,emergentcases=N/rho*i_n) -> cases
  cases %>% mutate(time=row_number()) -> cases
  cases
}

ms %>%  apply(MARGIN=1,FUN=infections)%>%
  bind_rows ->curves


ss_infections <- function(x){
  s0 = x["S0"] %>% as.numeric
  t1 <- x["theta.1"] %>% as.numeric
  t2 <- x["theta.2"] %>% as.numeric
  rho = x["rho"] %>% as.numeric
  initial <- c(s = s0,
               i=(1-s0),
               r = 0
  )
  
  parameters <- c(beta=t1,
                  gamma=t2
                  )
  # print(initial)
  # print(parameters)
  ode(y=initial,times=time,func=es_SIR,parms=parameters,method = "ode45")-> surveillance
  surveillance %>% as.data.frame %>%
    select(c("i")) %>%
    mutate(ecases=N/rho*i) -> cases
  cases %>% mutate(time=row_number()) -> cases
  cases
}

es %>% apply(MARGIN=1,FUN=ss_infections)%>%
  bind_rows ->curves_ss


os_infections <- function(x){
  s0 = x["S0"] %>% as.numeric
  t1 <- x["theta.1"] %>% as.numeric
  t2 <- x["theta.2"] %>% as.numeric
  rho = x["rho"] %>% as.numeric
  initial <- c(s = s0,
               i=(1-s0-p),
               r = 0,
               v =p
  )
  
  parameters <- c(beta=t1,
                  gamma=t2
  )
  # print(initial)
  # print(parameters)
  ode(y=initial,times=time,func=os_SIR,parms=parameters,method = "ode45")-> surveillance
  surveillance %>% as.data.frame %>%
    select(c("i")) %>%
    mutate(ocases=N/rho*i) -> cases
  cases %>% mutate(time=row_number()) -> cases
  cases
}

os %>% apply(MARGIN=1,FUN=os_infections)%>%
  bind_rows ->curves_os


curves %>% 
  mutate(olq = qpois(0.025,origcases),
            ouq= qpois(0.975,origcases),
            omedian = qpois(0.5,origcases),
            elq = qpois(0.025,emergentcases),
            euq= qpois(0.975,emergentcases),
            emedian = qpois(0.5,emergentcases)) %>%
  group_by(time) %>%
  summarise(olq = mean(olq),
            ouq= mean(ouq),
            omedian = mean(omedian),
            elq = mean(elq),
            euq= mean(euq),
            emedian = mean(emedian)) ->qlimits


curves_ss %>% 
  mutate(elq = qpois(0.025,ecases),
         euq= qpois(0.975,ecases),
         emedian = qpois(0.5,ecases)) %>%
  group_by(time) %>%
  summarise(elq = mean(elq),
            euq= mean(euq),
            emedian = mean(emedian)) ->ssqlimits

curves_os %>% 
  mutate(olq = qpois(0.025,ocases),
         ouq= qpois(0.975,ocases),
         omedian = qpois(0.5,ocases)) %>%
  group_by(time) %>%
  summarise(olq = mean(olq),
            ouq= mean(ouq),
            omedian = mean(omedian)) ->osqlimits
  

flu2008 %>% select(time,others,a_2009_h1n1)%>%
  left_join(qlimits,by="time")%>% 
  left_join(ssqlimits,by="time",suffix=c(".ms",".ss"))%>%
  left_join(osqlimits,by="time",suffix=c(".ms",".os")) ->mastermultistrain

mastermultistrain %>% 
  pivot_longer(cols=c("a_2009_h1n1","emedian.ms","emedian.ss"),names_to = "Strain",values_to = "Cases")%>%
  # head %>% View
  ggplot(aes(x=time,y=Cases,group=Strain,color=Strain,shape=Strain))+
  geom_jitter(size=1.5)+
  geom_ribbon(aes(ymin=elq.ms,ymax=euq.ms),color="green",alpha=0,size=1)+
  geom_ribbon(aes(ymin=elq.ss,ymax=euq.ss),color="blue",alpha=0,linetype=3,size=1)+

  theme_classic() + 
  xlab("Week") + 
  ylab("Reported Cases of A(H1N1)")+
 scale_color_discrete(labels=c("CDC data","Multi-strain Model","Single-strain model"),name="")+
 scale_shape_discrete(labels=c("CDC data","Multi-strain Model","Single-strain model"),name="")->h1n1plot

h1n1plot
ggsave(h1n1plot,filename="poisbasedquantiles_h1n1_withSSmodeljitter.png",
       width=8, height=6,dpi=300)

mastermultistrain %>% rename(a_others=others)%>%
  pivot_longer(cols=c("a_others","omedian.ms","omedian.os"),names_to = "Strain",values_to = "Cases")%>%
  ggplot(aes(x=time,y=Cases,group=Strain,color=Strain,shape=Strain))+
  geom_jitter(size=1.5)+
  geom_ribbon(aes(ymin=olq.ms,ymax=ouq.ms),color="green",alpha=0,size=1)+
  geom_ribbon(aes(ymin=olq.os,ymax=ouq.os),color="blue",alpha=0,linetype=3,size=1)+
  theme_classic() + 
  xlab("Week") + 
  ylab("Reported Cases of non-A(H1N1)")+
  scale_color_discrete(labels=c("CDC data","Multi-strain Model","Single-strain model"),name="")+
  scale_shape_discrete(labels=c("CDC data","Multi-strain Model","Single-strain model"),name="")->nonh1n1plot

nonh1n1plot

ggsave(nonh1n1plot,filename="poisbasedquantiles_others_withSSmodeljitter.png",
       width=8, height=6,dpi=300)
