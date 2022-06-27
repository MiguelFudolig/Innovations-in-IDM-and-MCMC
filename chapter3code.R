library(odin)
library(tidyverse)
library(ggplot2)
library(ggforce)
setwd("C:/Users/migsf/Google Drive/Migs Files/UNL life/Research/nonlinear dynamics/POMP/ODIN")
sir_generator <- odin::odin({
  ## Core equations for transitions between compartments:
  # update(S) <- S - n_SI
  # update(I) <- I + n_SI - n_IR
  # update(R) <- R + n_IR
  update(S) <-   S - dN_SI - dN_SInew;
  update(I) <- I + dN_SI - dN_IR;
  update(R) <-R + dN_IR - dN_RInew;
  update(V) <-  V - dN_VInew;
  update(Inew) <- Inew + dN_SInew + dN_VInew +dN_RInew - dN_IRnew;
  update(Rnew) <- Rnew + dN_IRnew;
  
  
  ## Individual probabilities of transition:
  p_SI <- (1 - exp(-beta * I / N))*(beta*I)/(beta*I+betanew*Inew) # S to I
  p_IR <- 1 - exp(-gamma) # I to R
  p_InS <- (1 - exp(-betanew * Inew / N))*(betanew*Inew)/(beta*I+betanew*Inew) # S to I 
  p_In <- 1 - exp(-betanew * Inew / N)
  p_IRn <- 1- exp(-gammanew)
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  
  dN_S[] <- rmultinom(S,pS) #multinomial to choose between infection of original and emergent strain
  pS[1] <- p_SI
  pS[2] <- p_InS
  pS[3] <- 1-p_SI-p_InS
  dim(pS) <- 3
  dim(dN_S) <- 3
  dN_SI <- dN_S[1]
  dN_IR <- rbinom(I, p_IR)
  dN_VInew <- rbinom(V, p_In)
  dN_SInew <- dN_S[2]
  dN_RInew <- rbinom(R,p_In)
  dN_IRnew <- rbinom(Inew,p_IRn)
  
  
  ## Total population size
  N <- S + I + R + V + Inew + Rnew
  
  ## Initial states:
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- R_ini
  initial(V) <- V_ini
  initial(Inew) <- Inew_ini
  initial(Rnew) <- N-S_ini - I_ini - R_ini - V_ini - Inew_ini
  
  
  
    
  ## User defined parameters - default in parentheses:
  S_ini <- user(49989)
  I_ini <- user(1)
  R_ini <- user(0)
  V_ini <- user(50000)
  Inew_ini <- user(10)
  beta <- user(1.3)
  gamma <- user(1)
  betanew <- user(2)
  gammanew<- user(1)
  
}, verbose = FALSE)

emergence<- odin::odin({
  ## Core equations for transitions between compartments:
  # update(S) <- S - n_SI
  # update(I) <- I + n_SI - n_IR
  # update(R) <- R + n_IR
  update(S) <-   S - dN_SI - dN_SInew;
  update(I) <- I + dN_SI - dN_IR;
  update(R) <-R + dN_IR - dN_RInew;
  update(V) <-  V - dN_VInew;
  update(Inew) <- Inew + dN_SInew + dN_VInew +dN_RInew - dN_IRnew;
  update(Rnew) <- Rnew + dN_IRnew;
  
  
  ## Individual probabilities of transition:
  p_SI <- 1 - exp(-beta * I / N) # S to I
  p_IR <- 1 - exp(-gamma) # I to R
  p_In <- 1 - exp(-betanew * Inew / N) # S to I 
  p_IRn <- 1- exp(-gammanew)
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  
  dN_S[] <- rmultinom(S,pS) #multinomial to choose between infection of original and emergent strain
  pS[1] <- p_SI
  pS[2] <- p_In
  pS[3] <- 1-p_SI-p_In
  dim(pS) <- 3
  dim(dN_S) <- 3
  dN_SI <- dN_S[1]
  dN_IR <- rbinom(I, p_IR)
  dN_VInew <- rbinom(V, p_In)
  dN_SInew <- dN_S[2]
  dN_RInew <- rbinom(R,p_In)
  dN_IRnew <- rbinom(Inew,p_IRn)
  
  
  ## Total population size
  N <- S + I + R + V + Inew + Rnew
  
  ## Initial states:
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- R_ini
  initial(V) <- V_ini
  initial(Inew) <- Inew_ini
  initial(Rnew) <- N-S_ini - I_ini - R_ini - V_ini - Inew_ini
  
  
  
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(49989)
  I_ini <- user(1)
  R_ini <- user(0)
  V_ini <- user(50000)
  Inew_ini <- user(10)
  beta <- user(1.3)
  gamma <- user(1)
  betanew <- user(2)
  gammanew<- user(1)
  
}, verbose = FALSE)

surv<- function(N,p,q,beta,gamma,betanew,gammanew,t=20)
{
S_ini=N*(1-p)*(1-q)
I_ini = N*(1-p)*q
V_ini =N*p
Inew_ini = 0
sir <- sir_generator(S_ini=N-I_ini-V_ini - Inew_ini, 
                     I_ini = I_ini,
                     V_ini =V_ini,
                     R_ini = 0,
                     Inew_ini=Inew_ini,
                     beta=beta,
                     gamma=gamma,
                     betanew=betanew,
                     gammanew=gammanew)


#sir


set.seed(20210429)
sir_col <- c("black", "blue", "pink","green","red","purple")

res <- sir$run(0:t)
res


}


survemerge<- function(N,p,pre,Inew_prop=0.01,beta,gamma,betanew,gammanew,t=20)
{
  pre[nrow(pre),]-> init
  Inew_ini = round(Inew_prop*as.numeric(init["S"]))
  S_ini=as.numeric(init["S"])-Inew_ini
  I_ini = as.numeric(init["I"])
  V_ini =as.numeric(init["V"])
  R_ini=as.numeric(init["R"])
  sir <- sir_generator(S_ini=S_ini, 
                       I_ini = I_ini,
                       V_ini =V_ini,
                       R_ini = R_ini,
                       Inew_ini=Inew_ini,
                       beta=beta,
                       gamma=gamma,
                       betanew=betanew,
                       gammanew=gammanew)
  
  # sir <- sir_generator(I_ini = 10)
  sir
  
  
  set.seed(20210429)
  sir_col <- c("black", "blue", "pink","green","red","purple")
  
  res <- sir$run((init["step"]):(init["step"]+t))
  postemerge <- bind_rows(as_tibble(pre[-nrow(pre),]),as_tibble(res))
  postemerge
}

N <- 10000000
p <- 0.5
q <- 0.01
beta <- 0.66#1.3 or 0.3
gamma <- 0.3
betanew <- 0.4
gammanew <- 0.3

rm(eko,emeko)
emergetime=10
surv(N,p,q,beta,gamma,betanew,gammanew,t=emergetime)->eko
survemerge(N,p,pre=eko,Inew_prop=0.01,beta,gamma,betanew,gammanew)->emeko
emeko %>%mutate(S=S/N,I=I/N,V=V/N,R=R/N,Inew=Inew/N,Rnew=Rnew/N)%>%pivot_longer(c("S", "I", "R","V","Inew","Rnew"),names_to = "Compartments",values_to = "SData") %>%
  ggplot(aes(x=step,y=SData,color=Compartments))+geom_path(size=1)+ geom_point(aes(shape=Compartments),size=2)+
  geom_vline(xintercept=emergetime,linetype="dashed",color="black") +
  xlab("Day") + ylab("Proportion of Population")+theme_classic()#+facet_zoom(xlim=c(10,30),ylim=c(0,0.1))

ggsave(filename=paste("beta,gamma,betanew,gammanew=(",beta,gamma,betanew,gammanew,") 10mil.png",sep=" "),
       dpi=300, width=6, height=6,units = "in")

emeko %>% mutate(emergeprop= Inew/(I+Inew)) %>% 
  ggplot(aes(x=step,y=emergeprop))+geom_step(size=2) + 
  xlab("Day") + ylab("Proportion of Emergent Strain")+
  ylim(c(0,1))+
  geom_hline(yintercept=0.5,linetype="dashed")+
  theme_classic()
ggsave(filename=paste("beta,gamma,betanew,gammanew=(",beta,gamma,betanew,gammanew,") infection proportion 10mil.png",sep=" "),
       dpi=300, width=6, height=6,units = "in")


##NEW EMERGENT SIZES

N <- 10000000
p <- 0.5
q <- 0.01
beta <- 0.8#1.3 or 0.3
gamma <- 0.3
betanew <- 0.33
gammanew <- 0.3

rm(emeko,eko)
surv(N,p,q,beta,gamma,betanew,gammanew,t=emergetime)->eko
survemerge(N,p,pre=eko,Inew_prop=0.01,beta,gamma,betanew,gammanew)->emeko
emeko$emerge<-0.01
survemerge(N,p,pre=eko,Inew_prop=0.05,beta,gamma,betanew,gammanew)->emeko5
emeko5$emerge<-0.05
survemerge(N,p,pre=eko,Inew_prop=0.1,beta,gamma,betanew,gammanew)->emeko10
emeko10$emerge<-0.1
survemerge(N,p,pre=eko,Inew_prop=0.005,beta,gamma,betanew,gammanew)->emekohalf
emekohalf$emerge<-0.005
survemerge(N,p,pre=eko,Inew_prop=0.025,beta,gamma,betanew,gammanew)->emeko25
emeko25$emerge<-0.025



bind_rows(emekohalf,emeko,emeko25,emeko5,emeko10) ->emergesize
emergesize %>% mutate(emergeprop=Inew/(I+Inew)) %>%
  ggplot(aes(x=step,y=emergeprop,color=as.factor(emerge)))+ 
  geom_point(aes(shape=as.factor(emerge)),size=3)+geom_step(size=1) +
  xlab("Day") + ylab("Proportion of Emergent Strain")+
  labs(color="Initial Emergence Proportion",shape="Initial Emergence Proportion")+
  ylim(c(0,1))+xlim(c(9,30))+
  geom_hline(yintercept=0.5,linetype="dashed")+
  theme_classic()

ggsave(filename=paste("beta,gamma,betanew,gammanew=(",beta,gamma,betanew,gammanew,") emergence sizes 10mil.png",sep=" "),
       dpi=300, width=6, height=6,units = "in")


N <- 10000000
p <- 0.5
q <- 0.01
beta <- 0.66#1.3 or 0.3
gamma <- 0.3
betanew <- 0.4
gammanew <- 0.3

rm(emeko,emeko5,emeko10,emekohalf,emeko25,eko)
surv(N,p,q,beta,gamma,betanew,gammanew,t=emergetime)->eko
survemerge(N,p,pre=eko,Inew_prop=0.01,beta,gamma,betanew,gammanew)->emeko
emeko$emerge<-0.01
survemerge(N,p,pre=eko,Inew_prop=0.05,beta,gamma,betanew,gammanew)->emeko5
emeko5$emerge<-0.05
survemerge(N,p,pre=eko,Inew_prop=0.1,beta,gamma,betanew,gammanew)->emeko10
emeko10$emerge<-0.1
survemerge(N,p,pre=eko,Inew_prop=0.005,beta,gamma,betanew,gammanew)->emekohalf
emekohalf$emerge<-0.005
survemerge(N,p,pre=eko,Inew_prop=0.025,beta,gamma,betanew,gammanew)->emeko25
emeko25$emerge<-0.025



bind_rows(emekohalf,emeko,emeko25,emeko5,emeko10) ->emergesize
emergesize %>% mutate(emergeprop=Inew/(I+Inew)) %>%
  ggplot(aes(x=step,y=emergeprop,color=as.factor(emerge)))+ 
  geom_point(aes(shape=as.factor(emerge)),size=3)+geom_step(size=1) +
  xlab("Day") + ylab("Proportion of Emergent Strain")+
  labs(color="Initial Emergence Proportion",shape="Initial Emergence Proportion")+
  ylim(c(0,1))+xlim(c(9,30))+
  geom_hline(yintercept=0.5,linetype="dashed")+
  theme_classic()

ggsave(filename=paste("beta,gamma,betanew,gammanew=(",beta,gamma,betanew,gammanew,") emergence sizes 10mil.png",sep=" "),
       dpi=300, width=6, height=6,units = "in")


infectious<- function(N,p,beta,gamma,betanew,gammanew,mu,filename)
{
  S_ini=round((gamma+mu)/beta*N,0)
  I_ini = max(0,round(N*mu*(beta*(1-p)-gamma-mu)/(beta*(gamma+mu)),0))
  V_ini =N*p
  Inew_ini = 1000
  sir <- sir_generator(S_ini=N-I_ini-V_ini - Inew_ini, 
                       I_ini = I_ini,
                       V_ini =V_ini,
                       R_ini = 0,
                       Inew_ini=Inew_ini,
                       beta=beta,
                       gamma=gamma,
                       betanew=betanew,
                       gammanew=gammanew)
  
  # sir <- sir_generator(I_ini = 10)
  sir
  

  set.seed(20210429)
  sir_col <- c("black", "blue", "pink","green","red","purple")

  res <- sir$run(0:20)
  head(res,20)

  png(file=filename,
      width=6, height=6,units="in",res=300)
  # par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  matplot(res[, 3], res[, 6], xlab = "Original infections", ylab = "Emergent infections",
          type = "l", col = sir_col, lty = 1,lwd=3,main="Surveillance data for N=100000")
  legend(x=12,y=70000, inset=c(-20,0),lwd = 1, col = sir_col, legend = c("S", "I", "R","V","I2","R2"), bty = "n",ncol=2)
  dev.off()
}



N <- 100000
p <- 0.3
beta <- 3 #1.3 or 0.3
gamma <- 0.6
betanew <- 1.3
gammanew <- 0.6
mu<-0.2

infectious(N,p,beta,gamma,betanew,gammanew,mu,filename="infectious3.png")