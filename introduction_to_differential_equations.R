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

# resets R to fresh
rm(list = ls(all = TRUE)) 

#change working directory to model file path.
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir) 

library(deSolve)
library(reshape)
library(plyr)

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
head(run_results)


# Plot
p<- ggplot() + geom_line(data=run_results, aes(x=time ,y=value,colour=gamma,group=gamma)) +
  geom_point(data=ddata, aes(x=time,y=value)) +
  xlab("time") + scale_y_continuous(labels = scales::comma) +
  ylab("number of individuals") + theme(axis.text.y  = element_text(angle=45))

show(p)