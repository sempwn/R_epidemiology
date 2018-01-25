############################################################################################
###########    Simple SIR model dynamics                                  ##################
###########                                                               ##################
###########                                                               ##################
###########                                                               ##################
###########    Written by: Mike A Irvine                                  ##################
###########                                                               ##################
###########    written: 11/04/17                                          ##################
############################################################################################
#If running for first time and neccessary, uncomment following lines to install packages 
#install library to handle differential equations
#install.packages("deSolve")
#install.packages("ggplot2")

# resets R to fresh
rm(list = ls(all = TRUE))  

# load in library for solving differential equations
library("deSolve")

# load in libraries for plotting
library("ggplot2")
library("reshape2")
library("scales")

SIR_func=function(t, x, vparameters){

  S = x[1]  
  I = x[2]  
  R = x[3]  
  
  with(as.list(vparameters),{
    
    npop = S+I+R   
    
    dS = -beta*S*I/npop +delta*R         
    dI = +beta*S*I/npop - gamma*I  
    dR = +gamma*I -delta*R                
    
    vout = c(dS,dI,dR)
    list(vout)
  })
}


npop = 100000
I_0 = 1
R_0 = 0
S_0 = npop-I_0-R_0

tbegin = 0
tend   = 365
vt = seq(tbegin,tend,1)  


gamma = 1/3 #recovery rate (days)        
R0    = 2.0 #Basic reproduction number    
beta  = R0*gamma #infectivity: product of R0 and recovery rate.
delta = 0.

vparameters = c(gamma=gamma,beta=beta,delta=delta)
inits = c(S=S_0,I=I_0,R=R_0)

model.results = as.data.frame(lsoda(inits, vt, SIR_func, vparameters))
meltdf <- melt(model.results,id="time")

p<- ggplot(meltdf, aes(x=time ,y=value,colour=variable,group=variable)) + geom_line() +
    xlab("time") + scale_y_continuous(labels = scales::comma) +
    ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45))

show(p)

#################################################
# Second model: Include vaccination and treatment
#################################################

SIR_VT_func=function(t, x, vparameters){
  #define the states of the model (susceptible, infected and recovered)
  S = x[1]  
  I = x[2]  
  R = x[3]  
  V = x[4]
  T = x[5]
  
  with(as.list(vparameters),{
    #total population size.
    npop = S+I+R+V+T  
    
    #Instantaneous change in states
    dS = -beta*S*I/npop +delta*R - vr*S       
    dI = +beta*S*I/npop - gamma*I -tr*I
    dR = +gamma*I -delta*R + te*T
    dV = +vr*S
    dT = +tr*I - te*T
    
    #output instananeous change in states as list. 
    vout = c(dS,dI,dR,dV,dT)
    list(vout)
  })
}

npop = 100000
I_0 = 1
R_0 = 0
V_0 = 0
T_0 = 0
S_0 = npop - I_0 - R_0 - V_0 - T_0

tbegin = 0
tend   = 365
vt = seq(tbegin,tend,1)  


gamma = 1/3 #recovery rate (1/days)        
R0    = 2.1 #Basic reproduction number (no dimensions)   
beta  = R0*gamma #infectivity: product of R0 and recovery rate.
delta = 0. #waning immunity
te = 1/7 #treatment efficacy (1/days)
tr = 0.0 #rate of treatment for infected (1/days)
vr = 1/365 #rate of vaccination per person (1/days)

vparameters = c(gamma=gamma,beta=beta,delta=delta,te=te,tr=tr,vr=vr)
inits = c(S=S_0,I=I_0,R=R_0,V=V_0,T=T_0)

model.results = as.data.frame(lsoda(inits, vt, SIR_VT_func, vparameters))
meltdf <- melt(model.results,id="time")

p<- ggplot(meltdf, aes(x=time ,y=value,colour=variable,group=variable)) + geom_line() +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45))

show(p)

#################################################
# Calculate critical vaccination rate ###########
#################################################
vrs <- seq(1.,365.,length=100)
vrs <- 1/vrs
tot_Is <- array(0,c(100))
for(i in 1:100){
  vr = vrs[i]
  vparameters = c(gamma=gamma,beta=beta,delta=delta,te=te,tr=tr,vr=vr)
  model.results = as.data.frame(lsoda(inits, vt, SIR_VT_func, vparameters))
  tot_Is[i] <- sum(model.results['I'])
}
vaccination.results <- data.frame(vaccination=1/vrs, infections=tot_Is)
p<- ggplot(vaccination.results, aes(x=vaccination ,y=infections)) + geom_line() +
  xlab("Expected time to vaccinate (days)") + scale_y_continuous(labels = scales::comma) +
  ylab("Total number of individuals infected") + theme(axis.text.y  = element_text(angle=45))

show(p)

