## 665 project simulation
library(dplyr)

#simulate dataset for one species
data <- data.frame("year" = sample(seq(1990, 2000, by = 1), 100, replace = T),
                   "route" = sample(seq(1,10, by = 1), 100, replace = T),
                   "count" = rpois(100, 2), 
                   "observer" = sample(seq(1,20, by = 1), 100, replace = T))
  
#observer-route ID
data$obsroute <- paste0(data$observer, data$route, sep = "")

#add dummy variable for first time observers
data <- arrange(data, year)

firstyr <- c(1)
for(i in 2:nrow(data)) {
  d <- data[i, ]
  d2 <- data[1:(i-1), ]
  uniq <- distinct(d, obsroute)
  uniq2 <- distinct(d2, obsroute)
  if (uniq$obsroute %in% uniq2$obsroute) {
    firstyr <- c(firstyr, 0)
  } else firstyr <- c(firstyr, 1)
}

data$firstyr <- firstyr

#add stratum ID
data <- data %>%
  mutate(strata = case_when(route >= 1 & route <= 3 ~ 1,
                       route > 3 & route < 8 ~ 2,
                       route >= 8 ~ 3))

#JAGS model of counts with random effects for year and obsroute
i <- data$strata
t <- data$year
j <- data$obsroute

#Model format
#log(counts) = Si + Bt(t - t*) + wj + nI(j,t) + Yit + Eijt
#poisson distribution

## Hyperparameters
#w ~ N(0, sigmaW) observer
#Y ~ N(0, sigmaYi) year
#E ~ N(0, sigmaE) overdispersion

#sigmaW ~ IG(0, 1000)
#sigmaYi ~ IG(0, 1000) *allowed to vary among strata
#sigmaE ~ IG(0, 1000)

## Parameters
#S ~ N(0, 1000) - strata specific intercept
#B ~ N(0, 1000) - strata specific slope
#n ~ N(0, 1000) - first year observer

counts <- model.frame(count ~ strata + year + obsroute + firstyr, data)
X <- model.matrix(count ~ strata + year + obsroute + firstyr, counts)
y <- counts$count
n <- length(y)


