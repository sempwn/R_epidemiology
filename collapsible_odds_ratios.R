
library(tidyverse)
library(epitools)


expit <- function(x) {
  1 / (1 + exp(-x))
}

d <- tibble(x = rbernoulli(100, p = 0.5))
d$z <- rbinom(100, size = 1, p = expit(0.1 + 0.0 * d$x))
d$y <- rbinom(100, size = 1, p = expit(0.5 + 0.5 * d$x -  2 * d$z))

ftable(d[, c("x", "y")])

od <- oddsratio(ftable(d[d$z == 1, c("x", "y")]))
or_z <- od$measure[2, 1]

od <- oddsratio(ftable(d[d$z == 0, c("x", "y")]))
or_nz <- od$measure[2, 1]


# marginal sum
(or_nz * table(d$z)[[1]] + or_z * table(d$z)[[2]]) / 100

#calculate odds ratio
od <- oddsratio(ftable(d[, c("x", "y")]))
od$measure[2, 1]

# unadjusted OR
m1 <- glm(y ~ x, data = d, family = "binomial")
exp(coef(m1))["xTRUE"]

# adjusted OR
m2 <- glm(y ~ x + z, data = d, family = "binomial")
exp(coef(m2))["xTRUE"]
