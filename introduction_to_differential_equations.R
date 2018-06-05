############################################################################################
###########    Introduction to Differential                               ##################
###########       Equations                                               ##################
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
#install.packages("diagram")

# resets R to fresh
rm(list = ls(all = TRUE)) 

#change working directory to model file path.
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir) 

library(deSolve) # package for solving ode (state) models
library(reshape) # manipulate dataframe
library(ggplot2) # data visualisation
library(plyr) #data frame manipulation
library(diagram)

#draw SIR diagram
openplotmat()
pos <- coordinates(c(3))
mid1 <- straightarrow(from = pos[1, ], to = pos[2, ])
mid2 <- straightarrow(from = pos[2, ], to = pos[3, ])
curvedarrow (from = pos[2, ], to = .5*(pos[1, ] +pos[2,]),curve = -2.0, lty = 6, lcol = 1)
text(mid1[1]-0.00, mid1[2]+0.1, expression(paste(beta, 'I')), cex = 1.5)
text(mid2[1]-0.00, mid2[2]+0.1, expression(gamma), cex = 1.5)
labels <- c('S','I','R')
text_size <- 1.3
edge_length <- 0.08
for(i in 1:length(labels)){
  textrect(mid = pos[i,], radx = edge_length, rady = edge_length, lab = labels[i], cex = text_size, box.col = "#0072B2")
}

#load data
data <- read.csv("./data/data_sir.csv")

# set initial conditions
init <- c(S = 0, I = 0, R = 0)

# set the time at which to get output (every year for forty years).
times <- seq(0, 40, by = 1)

# create SIR function.
sir <- function(time, state, parameters) {
  # {{Curly brackets indicate the beginning and end of functions}} This is the
  # meat of the function - what happens to the parameters. Here it needs to be
  # a list}}
  with(as.list(c(state, parameters)), {
    dS = -beta * S * I  # The change in S 
    dI = beta * S * I - gamma * I
    dR = gamma * I
    return(list(c(dS, dI, dR)))  # {{'Return' as in 'spit out' the values that you want}}
  })
}

# set list of parameters
parameters <- c(beta = 2/1e+05, gamma = 0.5)


# run model
init <- c(S = 1e+05 - 1, I = 1, R = 0)
out <- ode(y = init, func = sir, times = times, parms = parameters)
print(head(out))

melt.out <- melt(as.data.frame(out),id="time")
p<- ggplot(melt.out, aes(x=time ,y=value,colour=variable,group=variable)) + geom_line() +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45))

show(p)


# Does this match up well to data?

#change data into a proportion.
data$value = data$proportion * 1e+05
ddata <- as.data.frame(data)

p<- ggplot() + geom_line(data=melt.out, aes(x=time ,y=value,colour=variable,group=variable)) +
  geom_point(data=ddata, aes(x=time,y=value)) +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45))

show(p)


# Run for many values of gamma
gamma_values <- seq(0.2, 0.8, 0.1)  # vector of possible gamma values

run_results <- data.frame(time=double(),value=double(),gamma=double())
# Run through and solve with each gamma value Length = how many entries in
# vector What gamma value in this run?
for (i in 1:length(gamma_values)) {
  gamma_new <- gamma_values[i]
  # New parameter set
  parameters <- c(beta = 2/1e+05, gamma = gamma_new)
  # Run with this new parameter set
  out <- ode(y = init, times = times, func = sir, parms = parameters)
  # Store the number of infecteds this produces
  stored <- data.frame("time"=out[, "time"], "value" = out[, "I"], "gamma" = gamma_new)  # Bind the columns with time, the infected and the new gamma
  run_results <- rbind(run_results,stored)
}

# Have a look at what you have...(again wouldn't usually include this in the
# saved R code)
print(head(run_results))


# Plot
p<- ggplot() + geom_line(data=run_results, aes(x=time ,y=value,colour=gamma,group=gamma)) +
  geom_point(data=ddata, aes(x=time,y=value)) +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45))

show(p)


# Create SIR function with births and deaths
# This part of the script extends the original model to include births 
# in the susceptible population and deaths in the S,I,R population

#draw SIR diagram with births and deaths
openplotmat()
pos <- coordinates(c(4,4))
mid0 <- straightarrow(from = pos[1, ], to = pos[2, ])
mid1 <- straightarrow(from = pos[2, ], to = pos[3, ])
mid2 <- straightarrow(from = pos[3, ], to = pos[4, ])

curvedarrow (from = pos[3, ], to = .5*(pos[2, ] +pos[3,]),curve = 2.0, lty = 6, lcol = 1)

text(mid0[1]-0.00, mid0[2]+0.1, expression(paste(nu, 'N')), cex = 1.5)
text(mid1[1]-0.00, mid1[2]-0.1, expression(paste(beta, 'I')), cex = 1.5)
text(mid2[1]-0.00, mid2[2]+0.1, expression(gamma), cex = 1.5)

mid0 <- straightarrow(from = pos[2, ], to = pos[6, ])
mid1 <- straightarrow(from = pos[3, ], to = pos[7, ])
mid2 <- straightarrow(from = pos[4, ], to = pos[8, ])

text(mid0[1]-0.02, mid0[2]+0.0, expression(nu), cex = 1.5)
text(mid1[1]-0.02, mid1[2]-0.0, expression(nu), cex = 1.5)
text(mid2[1]-0.02, mid2[2]+0.0, expression(nu), cex = 1.5)

labels <- c('S','I','R')
text_size <- 1.3
edge_length <- 0.08
for(i in 1:length(labels)){
  textrect(mid = pos[i+1,], radx = edge_length, rady = edge_length, lab = labels[i], cex = text_size, box.col = "#0072B2")
}

sirbd <- function(time, state, parameters) {
  # {{Curly brackets indicate the beginning and end of functions}} This is the
  # meat of the function - what happens to the parameters. Here it needs to be
  # a list}}
  with(as.list(c(state, parameters)), {
    N = S + I + R
    dS = -beta * S * I + nu*N - nu*S  # The change in S 
    dI = beta * S * I - gamma * I - nu*I
    dR = gamma * I - nu*R
    return(list(c(dS, dI, dR)))  # {{'Return' as in 'spit out' the values that you want}}
  })
}

# run model
init <- c(S = 1e+05 - 1, I = 1, R = 0)
parameters <- c(beta = 2/1e+05, gamma = 0.5, nu = 0.02)
times <- seq(0, 100, by = 1)
out <- ode(y = init, func = sirbd, times = times, parms = parameters)
print(head(out))

melt.out <- melt(as.data.frame(out),id="time")
p<- ggplot(melt.out, aes(x=time ,y=value,colour=variable,group=variable)) + geom_line() +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('SIR with births and deaths')
  

show(p)


# Create SIRS function with births and deaths
# Here there is waning immunity where individuals who are recovered can end up becoming susceptible again.
# This part of the script extends the original model to include births 
# in the susceptible population and deaths in the S,I,R population

sirsbd <- function(time, state, parameters) {
  # {{Curly brackets indicate the beginning and end of functions}} This is the
  # meat of the function - what happens to the parameters. Here it needs to be
  # a list}}
  with(as.list(c(state, parameters)), {
    N = S + I + R
    dS = nu*N + delta * R -beta * S * I  - nu*S  # The change in S 
    dI = beta * S * I - gamma * I - nu * I
    dR = gamma * I - nu*R - delta * R
    return(list(c(dS, dI, dR)))  # {{'Return' as in 'spit out' the values that you want}}
  })
}

# run model
init <- c(S = 1e+05 - 1, I = 1, R = 0)
parameters <- c(beta = 2/1e+05, gamma = 0.5, nu = 0.02,delta=0.1)
times <- seq(0, 100, by = 1)
out <- ode(y = init, func = sirsbd, times = times, parms = parameters)
print(head(out))

melt.out <- melt(as.data.frame(out),id="time")
p<- ggplot(melt.out, aes(x=time ,y=value,colour=variable,group=variable)) + geom_line() +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45)) +
  ggtitle('SIRS with births and deaths')


show(p)

