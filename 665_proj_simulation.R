## 665 project simulation
library(dplyr)
source("clarkFunctions2018.r")

#simulate dataset for one species
data <- data.frame("year" = sample(seq(1990, 2000, by = 1), 1000, replace = T),
                   "route" = sample(seq(1,10, by = 1), 1000, replace = T),
                   "count" = rpois(100, 2), 
                   "observer" = sample(seq(1,20, by = 1), 1000, replace = T))
  
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
j <- data$obsroute
t <- data$year - 1989
data$t <- t

#Model format
#log(counts) = Si + Bt(t - t*) + wj + nI(j,t)
#poisson distribution
#t* is the first year of observations

## Hyperparameters
#w ~ N(0, sigmaW) observer

#sigmaW ~ IG(0, 1000)

## Parameters
#S ~ N(0, 1000) - strata specific intercept
#B ~ N(0, 1000) - strata specific slope
#n ~ N(0, 1000) - first year observer

#ID observers
z <- as.character(data$obsroute)
zall <- sort(unique(z))
ix <- match(z, zall)
data$ix <- ix

counts <- model.frame(count ~ strata + t + ix + firstyr, data)
X <- model.matrix(count ~ strata + t + ix + firstyr, counts)
y <- counts$count
n <- length(y)

#JAGS fixed effects
cat("model{
    for(i in 1:n){
      y[i] ~ dpois(lambda[i])
      lambda[i] <- exp(b1*X[i, 2] + b2*X[i, 3] + b3*X[i, 5])
    }

    b1 ~ dnorm(0.0, 1000)
    b2 ~ dnorm(0.0, 1000)
    b3 ~ dnorm(0.0, 1000)
    }
    ", fill = T, file = "countsFixed.txt"
    )

library(R2jags)
countData <- list(y = y, X = X, n = n)
parNames <- c("b1", "b2", "b3")
countFit <- jags(data = countData, param = parNames,
                 n.iter = 5000, n.burnin = 1000, model.file = "countsFixed.txt")
print(countFit)

par(mar=c(3,3,1,1))
plot( as.mcmc(countFit) )

#JAGS hierarchical
nobs <- max(ix)            # no. observers

cat("model{
    for(i in 1:n){
    y[i] ~ dpois(lambda[i])
    lambda[i] <- exp(b1*X[i, 2] + b2*X[i, 3] + b3*X[i, 5] + a1[ix[i]]*X[i, 4])
    }
    
    for(j in 1:nobs){
      a1[j] ~ dnorm(0.0, sigma)
    }
    
    sigma ~ dgamma(0.001, 1000)
    b1 ~ dnorm(0.0, 1000)
    b2 ~ dnorm(0.0, 1000)
    b3 ~ dnorm(0.0, 1000)
    }
    ", fill = T, file = "countsRE.txt"
)

countData   <- list(y = y, X = X, ix = ix, nobs = nobs, n = n)

parNamesRE <- c('b1','b2','b3','a1','sigma')

parInit <- function(){ 
  list(b1 = rnorm(1,0,1), b2 = rnorm(1,0,1), b3 = rnorm(1,0,1), sigma = runif(1,0,1),
       a1 = rnorm(nobs,0,1)) 
}

countRE <- jags(data=countData, inits=parInit, param=parNamesRE,
                 n.iter=2000, n.burnin=100, model.file="countsRE.txt")
rePars <- getJagsPars(countRE)
fixed  <- rePars$fixed
beta   <- fixed[1:4,1] #get very very small with random effect for simulated dataset
alpha  <- rePars$mean

