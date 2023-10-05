
library(tidyverse)
library(epitools)

set.seed(5)

expit <- function(x) {
  1 / (1 + exp(-x))
}

d <- tibble(x = rbinom(100, size = 1, p = 0.5))
d$z <- rbinom(100, size = 1, p = 0.5)
d$y <- rbinom(100, size = 1, p = expit(0.5 + 0.5 * d$x -  2 * d$z))

# contingency table
ftable(d[, c("x", "y")])

#calculate odds ratio
od <- oddsratio(ftable(d[, c("x", "y")]))
print("Total odds")
print(od$measure[2, 1])

od <- oddsratio(ftable(d[d$z == 1, c("x", "y")]))
or_z <- od$measure[2, 1]

od <- oddsratio(ftable(d[d$z == 0, c("x", "y")]))
or_nz <- od$measure[2, 1]

# marginal sum
print("marginal sum")
print((or_nz * table(d$z)[[1]] + or_z * table(d$z)[[2]]) / 100)



# unadjusted OR
m1 <- glm(y ~ x, data = d, family = "binomial")
print("Unadjusted OR")
print(exp(coef(m1))["x"])

# adjusted OR
m2 <- glm(y ~ x + z, data = d, family = "binomial")
print("adjusted OR")
print(exp(coef(m2))["x"])
