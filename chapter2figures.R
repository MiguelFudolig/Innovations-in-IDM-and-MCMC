library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggforce)
#Initializing population and birth/death rate mu andp
N = 1
mu = 0.5#0.2
p=0.5



#This function produces the plots with the values for beta, beta prime, gammma, and gamma prime as a PNG with filename fn
for_plots <- function(b0,bp0,gp,g,fn){
  S = c(round((1-p)*N,2))-0.01
  V = c(round((p*N),2))
  I1 = c(.01*N)
  I2 <- c(0)
  R2 <- c(0)
  R = c(N-S-I1-V)
  
  t=1
  
  #mutation
  
  #tau <- rgeom(1,0.01)
  #tau <- 50
  tau <- 100
  
  timeend <-200 
  while (t < timeend){
    b<- max(0,b0+ 0*rnorm(1,0,1))
    bp<- max(0,bp0+ 0*rnorm(1,0,1))
    
    if (t < tau){
      S[t+1] <- max(0,S[t] + (1-p)*mu -b*S[t]*I1[t] -mu*S[t])
      I1[t+1] <- max(0,I1[t] + b*S[t]*I1[t]- g*I1[t] - mu*I1[t])
      V[t+1] <- max(0,V[t] + p*mu- mu*V[t])
      R[t+1] <- 1-S[t+1]-I1[t+1]-V[t+1]
      I2[t+1] <- 0
      R2[t+1]<-0
      t <- t+1
    }
    if (t==tau){
      S[t+1] <- max(0,S[t]+ (1-p)*mu -b*S[t]*I1[t] -mu*S[t]-.01)
      I2[t+1] <- I2[t]+0.01 #1% of the people get infected
      I1[t+1] <- max(0,I1[t] + b*S[t]*I1[t]- g*I1[t] - mu*I1[t])
      R2[t+1] <- R2[t]
      V[t+1]<- max(0,V[t] + p*mu- mu*V[t])
      R[t+1]<-1-S[t+1]-I1[t+1]-V[t+1] - I2[t+1]-R2[t+1]
      t<- t+1
    }
    if(t>tau){
      S[t+1] <- max(0,S[t] + (1-p)*mu -b*S[t]*I1[t] -bp*S[t]*I2[t] -mu*S[t])
      I1[t+1] <- max(0,I1[t] + b*S[t]*I1[t]- g*I1[t] - mu*I1[t])
      V[t+1] <- max(0,V[t] + p*mu- bp*V[t]*I2[t]- mu*V[t])
      R[t+1] <- max(0,R[t]+g*I1[t] - bp*R[t]*I2[t] - mu*R[t])
      I2[t+1] <- min(max(0,I2[t]+bp*I2[t]*(S[t]+R[t]+V[t])-(gp+mu)*I2[t]),1)
      R2[t+1]<-max(0,1-S[t+1]-I1[t+1]-V[t+1] - I2[t+1]-R[t+1])
      #R2[t+1]<-max(0,R2[t]+(gp)*I2[t]-mu*R2[t])
      
      t<- t+1
    }
  }
  
  
  pop<-round(cbind(S,I1,R,V,I2,R2),2)
  sums <- apply(pop,1,sum)
  
  
  par(mfrow=c(2,3),oma=c(0,0,1,0))
  #"R0 = 6, R0' = 1.67, p =0.7: Endemic Equilibrium"
  plot(1:timeend,S, main=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =0.7.",sep=" "))
  plot(1:timeend,I1, main=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =0.7.",sep=" "))
  
  plot(1:timeend,V,main=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =0.7.",sep=" "))
  plot(1:timeend,R, main=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =0.7.",sep=" "))
  
  plot(1:timeend,I2, main=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =0.7.",sep=" "))
  plot(1:timeend,R2, main=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =0.7.",sep=" "))
  
  
  par(mfrow=c(1,1))
  head(R)
  S[1]
  V[1]
  
  R[1]
  head(V)
  
  
  surveillance <- as.data.frame(cbind(S,I1,V,R,I2))
  surveillance$time <- 1:nrow(surveillance)-1
  s_long <- gather(surveillance,Compartment,count,names(surveillance)[-6],factor_key=T)
  splots <- ggplot(s_long,aes(x=time,y=count,group=Compartment,color=Compartment)) + geom_line(aes(linetype=Compartment,color=Compartment))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
  
  splots <- splots + labs(x="Time",y="Proportion of Population",title=paste("R0 =", round(b/(g+mu),2), "R0' =", round(bp/(gp+mu),2), "p =",round(p,1),sep=" "))+scale_y_continuous(breaks=seq(0,0.8,by=0.1))
  splots <- splots + scale_color_manual(values=c("Black","Red","Blue","Green","Pink"),labels=c("s",expression('i'[1]),"v",expression('r'[1]),expression('i'[2])))+ scale_linetype_manual(values=c(1,2,1,8,12),labels=c("s",expression('i'[1]),"v",expression('r'[1]),expression('i'[2])))
  # splots <- splots + facet_zoom(xlim=c(0,5),ylim=c(0,0.015))
  splots <- splots + facet_zoom(xlim=c(100,105),ylim=c(0,0.015))
  ggsave(fn,height=6,width=8,units="in",dpi=300)
  # splots
}


#Running the simulation plots with each part having its own values for the
#transmission coefficients and removal coefficients.

#ENDEMIC EQUILIBRIUM
b0<-2
bp0 <-2
gp = 1
g=.1
for_plots(b0,bp0,gp,g,"newEE_test.png")

#EQ2EE

b0<-2
bp0 <-2
gp = .1
g=2/3-0.5

for_plots(b0,bp0,gp,g,"EQ2EE_test.png")
#EQ2DFE

b0<-2
bp0 <-2
gp = .1
g=2
for_plots(b0,bp0,gp,g,"EQ2DFE_test.png")

#EQ1

b0<-2 #3
bp0 <-2
gp = 1.2
g=.1



for_plots(b0,bp0,gp,g,"EE1_test.png")
#DFE

b0<-2
bp0 <-2
gp = 3
g=2

for_plots(b0,bp0,gp,g,"DFE_test.png")

b0<-2
bp0 <-2
gp = 3
g=2

for_plots(b0,bp0,gp,g,"DFE_test_zoom0.png")

b0<-2
bp0 <-2
gp = 3
g=2

for_plots(b0,bp0,gp,g,"DFE_test_zoom100.png")

#Producing the colormap for different values of the reproduction numbers for both original and emergent strain.

surv <- function(b0,bp0,gp,g){
  S = c(round((1-p)*N,2))-0.01
  V = c(round((p*N),2))
  I1 = c(.01*N)
  I2 <- c(0)
  R2 <- c(0)
  R = c(N-S-I1-V)
  
  t=1
  
  tau <- 100
  
  timeend <-300 
  while (t < timeend){
    b<- max(0,b0+ 0*rnorm(1,0,1))
    bp<- max(0,bp0+ 0*rnorm(1,0,1))
    
    if (t < tau){
      S[t+1] <- max(0,S[t] + (1-p)*mu -b*S[t]*I1[t] -mu*S[t])
      I1[t+1] <- max(0,I1[t] + b*S[t]*I1[t]- g*I1[t] - mu*I1[t])
      V[t+1] <- max(0,V[t] + p*mu- mu*V[t])
      R[t+1] <- 1-S[t+1]-I1[t+1]-V[t+1]
      I2[t+1] <- 0
      R2[t+1]<-0
      t <- t+1
    }
    if (t==tau){
      S[t+1] <- max(0,S[t]+ (1-p)*mu -b*S[t]*I1[t] -mu*S[t]-.01)
      I2[t+1] <- I2[t]+0.01 #1% of the people get infected
      I1[t+1] <- max(0,I1[t] + b*S[t]*I1[t]- g*I1[t] - mu*I1[t])
      R2[t+1] <- R2[t]
      V[t+1]<- max(0,V[t] + p*mu- mu*V[t])
      R[t+1]<-1-S[t+1]-I1[t+1]-V[t+1] - I2[t+1]-R2[t+1]
      t<- t+1
    }
    if(t>tau){
      S[t+1] <- max(0,S[t] + (1-p)*mu -b*S[t]*I1[t] -bp*S[t]*I2[t] -mu*S[t])
      I1[t+1] <- max(0,I1[t] + b*S[t]*I1[t]- g*I1[t] - mu*I1[t])
      V[t+1] <- max(0,V[t] + p*mu- bp*V[t]*I2[t]- mu*V[t])
      R[t+1] <- max(0,R[t]+g*I1[t] - bp*R[t]*I2[t] - mu*R[t])
      I2[t+1] <- min(max(0,I2[t]+bp*I2[t]*(S[t]+R[t]+V[t])-(gp+mu)*I2[t]),1)
      R2[t+1]<-max(0,1-S[t+1]-I1[t+1]-V[t+1] - I2[t+1]-R[t+1])
      #R2[t+1]<-max(0,R2[t]+(gp)*I2[t]-mu*R2[t])
      
      t<- t+1
    }
  }
  
  
  pop<-round(cbind(S,I1,R,V,I2,R2),2)
  sums <- apply(pop,1,sum)
  
  
  surveillance <- as.data.frame(cbind(S,I1,V,R,I2))
  surveillance$time <- 1:nrow(surveillance)-1
  surveillance[nrow(surveillance),]
}

epi_map <- function(len,b,bp,rlow,rhigh){
  # gg<- surv(b0,bp0,gp,g)
  # lol <- (gg$I2 > 0) + (gg$I1 > 0)
  # len <- 1000
  beta <- rep(b,len)
  rep<- seq(rlow,rhigh,length.out =len)
  betaprime <- rep(bp,len)
  reppr<- seq(rlow,rhigh,length.out =len)
  
  gamma<- beta/rep - mu
  gammaprime <- betaprime/reppr - mu
  
  states <- matrix(NA,nrow=length(gamma),ncol=length(gammaprime))
  
  for (j in 1:length(gamma)){
    for(k in 1:length(gammaprime)){
      endvalues <- surv(beta[j],betaprime[k],gammaprime[k],gamma[j])
      states[j,k] <- (endvalues$I1>1E-3) + 2*(endvalues$I2>1E-3)
    }
    
  }
  
  gg <- expand.grid(rep, reppr)
  names(gg)<- c("r0","rp")
  states_vec <- as.vector(states)
  states_vec <- as.factor(states_vec)
  forheatmap <- as.data.frame(cbind(gg,states_vec))
  forheatmap$rv <- forheatmap$r0*(1-p)
  heat <-ggplot(forheatmap,aes(r0,rp,fill=states_vec))+geom_tile()+scale_x_continuous(expand=c(0,0),name=expression(R[0]),breaks=0:10)+scale_y_continuous(expand=c(0,0),name=bquote(paste(R[0], "'")),breaks=0:10)+
    # scale_fill_discrete(colors=c("red","blue","green","yellow"))
    # geom_abline(slope=(1-p),intercept=0)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5),aspect.ratio = 1)+
    scale_fill_discrete(name = "Equilibrium point",labels = c("DFE","Original Strain","Emergent Strain","Endemic Equilibrium")) +
    ggtitle("Epidemic State at Equilibrium")
  heat
}

#simulating with beta = betaprime = 2, which was used in the manuscript.
p1 <- epi_map(1000,2,2,0.1,6)

ggsave(p1,height=5,width=5,dpi=300,file="heat_b2bp2.png")