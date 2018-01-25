############################################################################################
###########    Frequency versus Density                                   ##################
###########      population dynamics                                      ##################
###########                                                               ##################
###########                                                               ##################
###########    Written by: Mike A Irvine                                  ##################
###########                                                               ##################
###########    written: 17/01/25                                          ##################
############################################################################################
#If running for first time and neccessary, uncomment following lines to install packages 
#install library to handle differential equations
#install.packages("deSolve")
#install library to handle graphics
#install.packages("ggplot2")

# resets R to fresh
rm(list = ls(all = TRUE))  

# load in library for solving differential equations
library("deSolve")

# load in libraries for plotting
library("ggplot2")
library("reshape2")
library("scales")

density_func=function(t, x, vparameters){
  
  S = x[1]  
  I = x[2]  
  
  with(as.list(vparameters),{
    
     
    beta = kappa*nu
    dS = -beta*S*I/A      
    dI = +beta*S*I/A  
    
    vout = c(dS,dI)
    list(vout)
  })
}


npop = 2463431 #population of Metro Vancouver
I_0 = 1
S_0 = npop-I_0

tbegin = 0
tend   = 10
vt = seq(tbegin,tend,0.1)  


kappa = 0.1 #rate of contact per person per square km per day.
nu = 0.008 #probability that a contact leads to an infection
A = 2878.52 #area metro Vancouver (km^2)

vparameters = c(kappa=kappa,nu=nu,A=A) #set parameters
inits = c(S=S_0,I=I_0) #set initial conditions

model.results = as.data.frame(lsoda(inits, vt, density_func, vparameters))
meltdf <- melt(model.results,id="time",variable.name = "City")
meltdf <- meltdf[ which(meltdf$City=='I'), ]
meltdf$City = 'Vancouver'

p<- ggplot(meltdf, aes(x=time ,y=value,colour=City,group=City)) + geom_line() +
  xlab("Time (days)") + scale_y_continuous(labels = scales::comma) +
  ylab("Number of infected individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('Density-dependent transmission')

show(p)

##################################################################################
#####                                                                #############
#####            Compare density-dependent transmission              #############
#####                   in different cities                          #############
#####                                                                #############
##################################################################################

density_results_from_city = function(N,A,city_name){
  #N = population
  #A = Area in km^2
  npop = N
  I_0 = 1
  S_0 = npop-I_0
  
  tbegin = 0
  tend   = 10
  vt = seq(tbegin,tend,0.1)  
  
  
  kappa = 0.1 #rate of contact per person per square km per day.
  nu = 0.008 #probability that a contact leads to an infection
  
  vparameters = c(kappa=kappa,nu=nu,A=A) #set parameters
  inits = c(S=S_0,I=I_0) #set initial conditions
  
  model.results = as.data.frame(lsoda(inits, vt, density_func, vparameters))
  meltdf <- melt(model.results,id="time",variable.name = "City")
  meltdf <- meltdf[ which(meltdf$City=='I'), ]
  meltdf$City = city_name
  return(meltdf)
}

city_populations <- c(24256800,8874724,8787892,8443675,5928040,2463431,367770)
city_areas <- c(6340.5,1485.49,1572.15,783.84,5905.71,2878.52,696.15)
city_names <- c('Shanghai','Mexico City','London','New York','Toronto','Vancouver','Victoria')
cities_df <- data.frame(city_names,city_populations,city_areas)
colnames(cities_df)<-c('name','population','area')

p<- ggplot(cities_df, aes(x=area, y=population, fill=population)) +
  geom_label(aes(label=name), color="white", size=3,alpha=0.8) +
  scale_fill_continuous(guide = 'none') + 
  scale_y_continuous(label=comma_format(),
                     limits=c(0,25000000)) + 
  scale_x_continuous(label=comma_format(),limits=c(0,7000)) +
  theme(axis.text.y  = element_text(angle=45)) +
  xlab(bquote("Area ("*km^2*')')) + 
  ylab("Population size")
show(p)


df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df)<-c('time','City','value')

for(i in 1:nrow(cities_df)){
  pop <- cities_df$population[[i]]
  name <- cities_df$name[[i]]
  area <- cities_df$area[[i]]

  mdf <- density_results_from_city(pop,area,name)
  df <- rbind(df,mdf)
}

#####
###  Plot infections per city
#####


p<- ggplot(df, aes(x=time ,y=value,colour=City,group=City)) + geom_line() +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time (days)") + 
  ylab("Number of infected individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('Density-dependent transmission by city')
show(p)

#####
###  Plot infections per city in log scale.
#####

p<- ggplot(df, aes(x=time ,y=value,colour=City,group=City)) + geom_line() +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab("Time (days)") + 
  ylab("Number of infected individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('Density-dependent transmission by city')

show(p)

############################################################################################
###########                                                               ##################
###########                                                               ##################
###########                                                               ##################
###########                   Frequency transmission                      ##################
###########                                                               ##################
###########                                                               ##################
###########                                                               ##################
############################################################################################

frequency_func=function(t, x, vparameters){
  
  S = x[1]  
  I = x[2]  
  
  with(as.list(vparameters),{
    
    
    beta = eta*nu
    dS = -beta*S*I/N    
    dI = +beta*S*I/N
    
    vout = c(dS,dI)
    list(vout)
  })
}


N = 2463431 #population of Metro Vancouver
I_0 = 1
S_0 = N-I_0

tbegin = 0
tend   = 365
vt = seq(tbegin,tend,0.1)  

eta = 1 #rate of contact per person per day
nu = 0.008 #probability that a contact leads to an infection


vparameters = c(eta=eta,nu=nu,N=N) #set parameters
inits = c(S=S_0,I=I_0) #set initial conditions

model.results = as.data.frame(lsoda(inits, vt, frequency_func, vparameters))
meltdf <- melt(model.results,id="time",variable.name = "City")
meltdf <- meltdf[ which(meltdf$City=='I'), ]
meltdf$City = 'Vancouver'

p<- ggplot(meltdf, aes(x=time ,y=value,colour=City,group=City)) + geom_line() +
  xlab("Time (days)") + scale_y_continuous(labels = scales::comma) +
  ylab("Number of infected individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('Frequency-dependent transmission')
  
show(p)

#####
##
##      look at frequency for all cities
##
#####

frequency_results_from_city = function(N,A,city_name){
  #N = population
  #A = Area in km^2
  I_0 = 1
  S_0 = N-I_0
  
  tbegin = 0
  tend   = 365
  vt = seq(tbegin,tend,5)  
  
  eta = 1 #rate of contact per person per day
  nu = 0.008 #probability that a contact leads to an infection
  
  
  vparameters = c(eta=eta,nu=nu,N=N) #set parameters
  inits = c(S=S_0,I=I_0) #set initial conditions
  
  model.results = as.data.frame(lsoda(inits, vt, frequency_func, vparameters))
  meltdf <- melt(model.results,id="time",variable.name = "City")
  meltdf <- meltdf[ which(meltdf$City=='I'), ]
  meltdf$City = city_name
  return(meltdf)
}

df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df)<-c('time','City','value')

for(i in 1:nrow(cities_df)){
  pop <- cities_df$population[[i]]
  name <- cities_df$name[[i]]
  area <- cities_df$area[[i]]
  
  mdf <- frequency_results_from_city(pop,area,name)
  df <- rbind(df,mdf)
}

#####
###  Plot infections per city
#####


p<- ggplot(df, aes(x=time ,y=value,colour=City,group=City)) + geom_line() +
  scale_y_continuous(labels = scales::comma) +
  xlab("Time (days)") + 
  ylab("Number of infected individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('Frequency-dependent transmission by city')

show(p)


