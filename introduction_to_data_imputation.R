############################################################################################
###########           Introduction to data imputation                     ##################
###########                                                               ##################
###########                                                               ##################
###########             Written by: Mike A Irvine                         ##################
###########                                                               ##################
###########                written: 2018-06-04                            ##################
############################################################################################
#install.packages("GGally")
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally)


#Construct data
construct.data <- function(coefficient=c(5,2,3)){
  n <- 100
  x1 <- rnorm(100)
  x2 <- rnorm(100)
  x3 <- -0.1*x1 + rnorm(100)
  y <- 5*x1 + 2*x2 + 3*x3 + rnorm(100)
  return(data.frame(x1,x2,x3,y))
}

impute <- function (a, a.impute){ 
  ifelse (is.na(a), a.impute, a)
}

plot.df.histogram <-function(df){
  p <- ggplot(data=df, aes(y)) + geom_histogram(bins=30)
  show(p)
}

plot.compare.histogram <-function(y1,y2,labels=c("true","reconstructed")){
  df1 <- data.frame(y=y1)
  df2 <- data.frame(y=y2)
  df1$category <- labels[1]
  df2$category <- labels[2]
  
  df<-rbind(df1,df2)
  p<- ggplot(df, aes(y, fill = category))   + 
  geom_histogram(alpha = 0.2, position="identity")
  #geom_density(alpha = 0.2)

  show(p)
}

#construct data
df <- construct.data()

ggpairs(df)

# fit a simple linear regression to complete data
model.full = lm(y ~ x1 + x2 + x3, data = df)

summary(model.full)

###################################################
##                                               ##
##                                               ##
##               SECTION I                       ##
##           Missing at random                   ##
##                                               ##
###################################################

# Now let's replace 20 randomly selected outcome measures with NA.
x1 <- df$x1 #save full x1 values for later
df[sample(seq_len(nrow(df)), 20, replace = FALSE),1] <- NA


###################################################
##                Solution one                   ##
##                ignore all NA                  ##
###################################################

df.ignore.na <- df[complete.cases(df),]

plot.compare.histogram(x1,df.ignore.na$x1,labels=c("true","reconstructed from \nnot NA"))

model.ignore.na = lm(y ~ x1 + x2 + x3, data = df.ignore.na)

summary(model.ignore.na)



###################################################
##                Solution two                   ##
##               replace with mean               ##
###################################################

df$x1.mean <- df$x1
df$x1.mean[is.na(df$x1.mean)] <- mean(na.omit(df$x1))

plot.compare.histogram(x1,df$x1.mean,labels=c("true","reconstructed from \nmean"))

model.mean = lm(y ~ x1.mean + x2 + x3, data = df)

summary(model.mean)

###################################################
##              Solution three                   ##
##                imputation                     ##
###################################################

lm.imp <- lm(x1 ~ y + x2 + x3,data=df)

pred <- predict (lm.imp, df)
df$x1.imp <- impute (x1, pred)

plot.compare.histogram(x1,df$x1.imp,labels=c("true","reconstructed from \ndeterministic imputation"))

model.imp = lm(y ~ x1.imp + x2 + x3, data = df)

summary(model.imp)


###################################################
##                                               ##
##                                               ##
##              SECTION II                       ##
##       Missing not at random                   ##
##                                               ##
###################################################

df <- construct.data()

# probability that x1 value is missing dependent on y
p <- 1/(1+exp(-3*df$x1))
excluded <- runif(length(p)) < p

x1 <- df$x1

df[which(excluded),1] <- NA

###################################################
##                Solution one                   ##
##                ignore all NA                  ##
###################################################

df.ignore.na <- df[complete.cases(df),]

plot.compare.histogram(x1,df.ignore.na$x1,labels=c("true","reconstructed from \nnot NA"))

model.ignore.na = lm(y ~ x1 + x2 + x3, data = df.ignore.na)

summary(model.ignore.na)



###################################################
##                Solution two                   ##
##               replace with mean               ##
###################################################

df$x1.mean <- df$x1
df$x1.mean[is.na(df$x1.mean)] <- mean(na.omit(df$x1))

plot.compare.histogram(x1,df$x1.mean, labels=c("true","reconstructed from \nmean"))

model.mean = lm(y ~ x1.mean + x2 + x3, data = df)

summary(model.mean)

###################################################
##              Solution three                   ##
##                imputation                     ##
###################################################

lm.imp <- lm(x1 ~ y + x2 + x3,data=df)

pred <- predict (lm.imp, df)
df$x1.imp <- impute (x1, pred)

plot.compare.histogram(x1,df$x1.imp,labels=c("true","reconstructed from \ndeterministic imputation"))

model.imp = lm(y ~ x1.imp + x2 + x3, data = df)

summary(model.imp)




