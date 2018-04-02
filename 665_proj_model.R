## 665 final project
## Grace Di Cecco
## JAGS models for count data from BBS

library(dplyr)
library(geosphere)
library(maps)
library(rgdal)
library(rgeos)

source("clarkFunctions2018.r")

#Read in BBS data
##Note: on mac once connected to Bioark server, path is "/Volumes/hurlbertlab/Databases/BBS/2017/filename.csv"
routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")
bcrshp <- readOGR("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\bcr_terrestrial_shape\\BCR_Terrestrial_master.shp") #BCRs
##Typo in bcrshp - has TENNESSE

#Subset species - from Huang et al. 2017
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\BBS-centroids\\huang-2017-bbs-species.csv", header = TRUE)

#Remove routes weather runtype = 0
routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year", "obsn"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, year >= 1969, select = c("statenum","stateroute", "year", "obsn", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

library(R2jags)
results <- matrix(0, ncol = 27)
colnames(results) <- colnames(counts.merge)

for(aou in species$aou) {
#Subset counts by species
counts.short <- counts %>%
  filter(year >= 1969) %>%
  filter(aou == aou)
counts.merge <- merge(routes.short, counts.short, by = c("stateroute", "year"))

#Add obs-route ID and dummy variable for first year observers
counts.merge$obsroute <- paste0(counts.merge$observer, counts.merge$route, sep = "")

counts.merge <- arrange(counts.merge, year)

#Loop takes some time
init.time = Sys.time()
firstyr <- c(1)
for(i in 2:nrow(counts.merge)) {
  d <- counts.merge[i, ]
  d2 <- counts.merge[1:(i-1), ]
  uniq <- distinct(d, obsroute)
  uniq2 <- distinct(d2, obsroute)
  if (uniq$obsroute %in% uniq2$obsroute) {
    firstyr <- c(firstyr, 0)
  } else firstyr <- c(firstyr, 1)
  curr.time = Sys.time()
  elapsed = curr.time - init.time
  percelltime = elapsed/i
  estimated.end = (nrow(counts.merge) - i)*percelltime + curr.time
  print(paste(i, "out of", nrow(counts.merge), "; current time:", curr.time,
              "; estimated end time:", estimated.end))
}
#1 spp - about 10 minutes
counts.merge$firstyr <- firstyr

Sys.time()

#JAGS model of counts with random effects for year and obsroute
s <- as.character(paste0(counts.merge$statenum.x, counts.merge$bcr, sep = ""))
sall <- sort(unique(s))
i <- match(s, sall)
counts.merge$strata <- i

j <- counts.merge$obsroute
t <- counts.merge$year - 1968
counts.merge$t <- t

#ID observers
z <- as.character(counts.merge$obsroute)
zall <- sort(unique(z))
ix <- match(z, zall)
counts.merge$ix <- ix

counts <- model.frame(stoptotal ~ strata + t + ix + firstyr, counts.merge)
X <- model.matrix(stoptotal ~ strata + t + ix + firstyr, counts)
y <- counts$stoptotal
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

countData <- list(y = y, X = X, n = n)
parNames <- c("b1", "b2", "b3")
countFit <- jags(data = countData, param = parNames,
                 n.iter = 10000, n.burnin = 2500, model.file = "countsFixed.txt")
# 1 spp - < 2 hours
print(countFit)
countFit.out <- countFit$BUGSoutput$summary
write.csv(countFit.out, paste0(aou, "_fixed_jagssummary.csv", sep = ""))

Sys.time()

par(mar=c(3,3,1,1))
plot( as.mcmc(countFit) )

#JAGS random effects

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

parInit <- function() {
  list(b1 = rnorm(1,0,1), b2 = rnorm(1,0,1), b3 = rnorm(1,0,1)) 
}

Sys.time()
countRE <- jags(data=countData, inits=parInit, param=parNamesRE,
                n.iter=20000, n.burnin=5000, model.file="countsRE.txt")
#3.5 hours for bobwhite
Sys.time()
rePars <- getJagsPars(countRE)
fixed  <- rePars$fixed
beta   <- fixed[1:4,1] 
alpha  <- rePars$mean
write.csv(fixed, paste0(aou, "_beta_jagsRE.csv", sep = ""), row.names = T)
write.csv(alpha, paste0(aou, "_alpha_jagsRE.csv", sep = ""))

#calculate strata specific abundance indices
counts.merge$abundind <- exp(beta[1]*counts.merge$strata + beta[2]*counts.merge$t) #model with RE
fixedb1 <- countFit$BUGSoutput$mean$b1
fixedb2 <- countFit$BUGSoutput$mean$b2
counts.merge$fixedabundind <- exp(fixedb1*counts.merge$strata + fixedb2*counts.merge$t)

results <- rbind(results, counts.merge)
}

write.csv(results, "counts_w_modeloutput.csv", row.names = F)
