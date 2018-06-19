# Script for running fixed effects and hierarchical models on Longleaf

library(dplyr)
library(R2jags)

source("/proj/hurlbertlab/gdicecco/665_proj/clarkFunctions2018.r")

#Read in BBS data
routes <- read.csv("/proj/hurlbertlab/bbs/bbs_routes_20170712.csv")
bbscounts <- read.csv("/proj/hurlbertlab/bbs/bbs_counts_20170712.csv")
weather <- read.csv("/proj/hurlbertlab/bbs/bbs_weather_20170712.csv")

#Subset species - from Huang et al. 2017
species <- read.csv("/proj/hurlbertlab/gdicecco/665_proj/huang-2017-bbs-species.csv", header = TRUE)

#Remove routes weather runtype = 0
routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year", "obsn"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, year >= 1969, select = c("statenum","stateroute", "year", "obsn", "latitude", "longitude", "bcr"))
bbscounts$stateroute <- bbscounts$statenum*1000 + bbscounts$route

results <- matrix(0, ncol = 27)
colnames(results) <- c("stateroute", "year", "statenum.x", "obsn", "latitude", "longitude",
                       "bcr", "record_id", "countrynum", "statenum.y", "route", "rpid", "aou",
                       "count10", "count20", "count30", "count40", "count50", "stoptotal",
                       "speciestotal", "obsroute", "firstyr", "strata", "t", "ix", "abundind", "fixedabundind")
dics <- matrix(0, ncol = 3)
colnames(dics) <- c("aou", "jags", "dic")

pars.results <- matrix(0, ncol = 6)
colnames(pars.results) <- c("aou", "model", "param", "mean", "cilo", "cihi")

setwd("/proj/hurlbertlab/gdicecco/665_proj/")
for(i in 1:length(species$aou)) {
  tryCatch({
    AOU <- species$aou[i]
    
    #Subset counts by species
    counts.short <- bbscounts %>%
      filter(year >= 1969) %>%
      filter(aou == AOU)
    counts.merge <- merge(routes.short, counts.short, by = c("stateroute", "year"))
    
    #Add obs-route ID and dummy variable for first year observers
    counts.merge$obsroute <- paste0(counts.merge$obsn, counts.merge$stateroute, sep = "")
    
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
    
    counts.merge$firstyr <- firstyr
    
    #JAGS models of counts
    s <- as.character(paste0(counts.merge$statenum.x, counts.merge$bcr, sep = ""))
    sall <- sort(unique(s))
    is <- match(s, sall)
    counts.merge$strata <- is
    nstrata <- max(is)
    
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
        lambda[i] <- exp(b1[is[i]] + b2[is[i]]*X[i, 3] + b3*X[i, 5])
        }
        
        for(k in 1:nstrata) {
        b1[k] ~ dnorm(0.0, 1000)
        b2[k] ~ dnorm(0.0, 1000)
        }
        
        b3 ~ dnorm(0.0, 1000)
  }
        ", fill = T, file = "countsFixed.txt"
    )
    
    countData <- list(y = y, X = X, n = n, is = is, nstrata = nstrata)
    parNames <- c("b1", "b2", "b3")
    countFit <- jags(data = countData, param = parNames,
                     n.iter = 10000, n.burnin = 2500, model.file = "countsFixed.txt")
    
    countFit.out <- countFit$BUGSoutput$summary
    names <- rownames(countFit.out)
    countFit.df <- as.data.frame(countFit.out)
    parstmp <- data.frame(aou = AOU, 
                          model = "fixedEffects", 
                          param = names,
                          mean = countFit.df$mean,
                          cilo = countFit.df$'2.5%',
                          cihi = countFit.df$'97.5%')
    
    pars.results <- rbind(pars.results, parstmp)
    
    
    dictmp <- c(AOU, "fixed", countFit$BUGSoutput$DIC)
    dics <- rbind(dics, dictmp)
    
    #JAGS hierarchical
    
    nobs <- max(ix)            # no. observers
    
    cat("model{
        for(i in 1:n){
        y[i] ~ dpois(lambda[i])
        lambda[i] <- exp(b1[is[i]] + b2[is[i]]*X[i, 3] + b3*X[i, 5] + a1[ix[i]]*X[i, 4])
        }
        
        for(k in 1:nstrata) {
        b1[k] ~ dnorm(0.0, 1000)
        b2[k] ~ dnorm(0.0, 1000)
        }
        
        for(m in 1:nobs){
        a1[m] ~ dnorm(0.0, sigma)
        }
        
        sigma ~ dgamma(1, 1000)
        
        b3 ~ dnorm(0.0, 1000)
        }
        ", fill = T, file = "countsRE.txt"
    )
    
    countData   <- list(y = y, X = X, ix = ix, nobs = nobs, is = is, nstrata = nstrata, n = n)
    
    parNamesRE <- c('b1','b2','b3','a1','sigma')
    
    parInit <- function() {
      list(b1 = rnorm(nstrata,0,1), b2 = rnorm(nstrata,0,1), b3 = rnorm(1,0,1)) 
    }
    
    countRE <- jags(data=countData, inits=parInit, param=parNamesRE,
                    n.iter=40000, n.burnin=20000, model.file="countsRE.txt")
    
    rePars <- getJagsPars(countRE)
    
    beta   <- na.omit(rePars$mean[, 2:3])
    
    fixed.df <- as.data.frame(rePars$fixed)
    pars.df <- data.frame(aou = AOU, model = "randomEffects", 
                          param = c("b3", "deviance", "sigma", 
                                    rep("a1", length(rePars$mean[, 1])), rep("b1", length(rePars$mean[, 2])),
                                    rep("b2", length(rePars$mean[, 3]))),
                          mean = c(fixed.df$mean, rePars$mean[, 1], rePars$mean[, 2], rePars$mean[, 3]),
                          cilo = c(fixed.df$'2.5%', rePars$ciLo[, 1], rePars$ciLo[, 2], rePars$ciLo[, 3]),
                          cihi = c(fixed.df$'97.5%', rePars$ciHi[, 1], rePars$ciHi[, 2], rePars$ciHi[, 3]))
    
    pars.results <- rbind(pars.results, pars.df)
    
    dicRE <- c(AOU, "randomefx", countRE$BUGSoutput$DIC)
    dics <- rbind(dics, dicRE)
    
    #calculate strata specific abundance indices
    counts.merge$abundind <- exp(beta[counts.merge$strata, 1] + beta[counts.merge$strata, 2]*counts.merge$t) #model with RE
    fixedb1 <- countFit$BUGSoutput$mean$b1
    fixedb2 <- countFit$BUGSoutput$mean$b2
    counts.merge$fixedabundind <- exp(fixedb1[counts.merge$strata] + fixedb2[counts.merge$strata]*counts.merge$t)
    
    results <- rbind(results, counts.merge)
  }, error = function(i) {print("ERROR : SPECIES # ", "i", sep = "")})
  }

write.csv(dics[-1, ], "/proj/hurlbertlab/gdicecco/665_proj/DIC_all_spp_models.csv", row.names = F)
write.csv(results[-1, ], "/proj/hurlbertlab/gdicecco/665_proj/counts_w_modeloutput.csv", row.names = F)
write.csv(pars.results[-1, ], "/proj/hurlbertlab/gdicecco/665_proj/model_parameters.csv", row.names = F)
