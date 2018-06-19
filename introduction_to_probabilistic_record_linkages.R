############################################################################################
###########             Introduction to probabilistic                     ##################
###########                    record linkages                            ##################
###########                                                               ##################
###########             Author: Mike Irvine                               ##################
###########             Contact: m.irvine@math.ubc.ca                     ##################
###########                Date: 2018-06-18                               ##################
############################################################################################

#install.packages("RecordLinkage")
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RecordLinkage)

#create theme for table plots
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))

# create function to plot tables
plotTable <- function(df){
  frame()
  
  
  myt <- gridExtra::tableGrob(df, theme = mytheme)
  
  grid.draw(myt)
} 

# create example datasets

firstname <- c("Adrian","ASH","Yoav","Fiona")
lastname <- c("Sayers","Blom","BenShlomo","Steel")
sex <- c("m","m","m","f")
df1 = data.frame(firstname,lastname,sex, stringsAsFactors=FALSE)

firstname <- c("Fiona","Yoav","Ashley","Adrian")
lastname <- c("Steele","Ben-Shlomo","Blom","Sayers")
edu <- c("phd","phd","phd","msc")
df2 = data.frame(firstname,lastname,edu, stringsAsFactors=FALSE)

frame()
grid.table(df1)

frame()
grid.table(df2)

###################################################
##                                               ##
##                                               ##
##               SECTION I                       ##
##           Deterministic linkage               ##
##                                               ##
###################################################


#join on first and last name. Compare the inner join and outer join
df.innerjoin <- merge(x=df1, y=df2, by = c("firstname","lastname"))
df.outerjoin <- merge(x=df1, y=df2, by = c("firstname","lastname"), all=TRUE)

# cartesian product
df.cp <- merge(x = df1, y = df2, by = NULL)
df.cp <- transform(df.cp, firstname.match= firstname.x==firstname.y)
df.cp <- transform(df.cp, lastname.match= lastname.x==lastname.y)

plotTable(df.cp)

# calculate log ratio

calc.R <-function(c1,c2,m=0.95,u=0.25){
  c1 <- as.numeric(c1)
  c2 <- as.numeric(c2)
  M <- log2(m)*(c1+c2) +log2(1-m)*(1-c1) + log2(1-m)*(1-c2)
  U <- log2(u)*(c1+c2) +log2(1-u)*(1-c1) + log2(1-u)*(1-c2)
  return (M-U)
}

df.cp$R <- calc.R(df.cp$firstname.match,df.cp$lastname.match)


plotTable(df.cp[,c(1,2,3,4,5,6,7,8,9)])

###################################################
##                                               ##
##                                               ##
##               SECTION II                      ##
##           Probabilistic linkage               ##
##                                               ##
###################################################

#Use the Levenshtein (edit) distance on first and last names of each dataframe
df.cp <- merge(x = df1, y = df2, by = NULL)
df.cp$firstname.match <- RecordLinkage::levenshteinSim(df.cp$firstname.x,df.cp$firstname.y)
df.cp$lastname.match <- RecordLinkage::levenshteinSim(df.cp$lastname.x,df.cp$lastname.y)

head(df.cp)

calc.R <-function(c1,c2,m=0.95,u=0.25){
  c1 <- as.numeric(c1)
  c2 <- as.numeric(c2)
  R1 <- log2(m/u - ( m/u - (1-m)/(1-u) )*(1-c1))
  R2 <- log2(m/u - (m/u - (1-m)/(1-u))*(1-c2))
  return (R1+R2)
}

df.cp$R <- calc.R(df.cp$firstname.match,df.cp$lastname.match)


plotTable(df.cp[,c(1,2,3,4,5,6,9)])


#convert into a posterior probability

N.expected.matches <- 4
N.1 <- nrow(df1) 
N.2 <- nrow(df2)

p.M.prior <- N.expected.matches/(N.1*N.2)
prior.ratio <- p.M.prior/(1-p.M.prior)

df.cp$pM <- exp(log(2)*(df.cp$R + log2(prior.ratio)))
df.cp$pM <- df.cp$pM/(1+df.cp$pM)

plotTable(df.cp[,c(1,2,3,4,5,6,10)])


###################################################
##                                               ##
##                                               ##
##               SECTION III                     ##
##             Implemention in                   ##
##           RecordLinkage package               ##
##                                               ##
###################################################

#perform string-based comparison for columns one and two
rpairs <- compare.linkage(df1,df2,strcmp = c(1,2))

#calculate weights based on E-M
rpairs <- epiWeights(rpairs)

# summarise fitted weights including providing histogram
summary(rpairs)

# classify based on threshold
result <- epiClassify(rpairs, 0.55)

# find all matched pairs
matched <- result$pairs[result$prediction=='L',]

# generate resulting dataframe
result.df <- data.frame(cbind(df1[matched$id1,],df2[matched$id2,]))

plotTable(result.df)

