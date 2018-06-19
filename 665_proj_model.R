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
bbscounts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
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
bbscounts$stateroute <- bbscounts$statenum*1000 + bbscounts$route

library(R2jags)
results <- matrix(0, ncol = 27)
colnames(results) <- c("stateroute", "year", "statenum.x", "obsn", "latitude", "longitude",
                       "bcr", "record_id", "countrynum", "statenum.y", "route", "rpid", "aou",
                       "count10", "count20", "count30", "count40", "count50", "stoptotal",
                       "speciestotal", "obsroute", "firstyr", "strata", "t", "ix", "abundind", "fixedabundind")
dics <- matrix(0, ncol = 3)
colnames(dics) <- c("aou", "jags", "dic")

pars.results <- matrix(0, ncol = 6)
colnames(pars.results) <- c("aou", "model", "param", "mean", "cilo", "cihi")

spp15 <- species[1:15, ]

for(i in 1) {
  AOU <- spp15$aou[i]

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

counts <- model.frame(speciestotal ~ strata + t + ix + firstyr, counts.merge)
X <- model.matrix(speciestotal ~ strata + t + ix + firstyr, counts)
y <- counts$speciestotal
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
}

write.csv(dics[-1, ], "DIC_all_spp_models.csv", row.names = F)
write.csv(results[-1, ], "counts_w_modeloutput.csv", row.names = F)
write.csv(pars.results[-1, ], "model_parameters.csv", row.names = F)

#### Figures for presentation #####
# delta DIC

delta <- as.data.frame(dics[-1, ]) %>%
  group_by(aou) %>% 
  summarize(delta = diff(as.numeric(as.character(dic))))

library(ggplot2)
theme_set(theme_bw())

ggplot(delta, aes(x = aou, y = log10(abs(delta)))) + geom_bar(stat = "identity", fill = "grey") + 
  xlab("AOU") +  ylab("log10(Fixed effects DIC - Random effects DIC)")

# distribution of coefficients
## AOU == 3850

# strata specific intercept
# 34 strata
betas <- as.data.frame(exp(beta), row.names = F)

library(cowplot)

# strata specific fixed effects (intercept, slope)
b1 <- ggplot(betas, aes(x = b1)) + geom_histogram(binwidth = 0.0001, col = "white") +
  ylab("Number of strata") + xlab("Intercept")
b2 <- ggplot(betas, aes(x = b2)) + geom_histogram(binwidth = 0.005, col = "white") +
  ylab("Number of strata") + xlab("Slope")
plot_grid(b1, b2, nrow = 2)

# compare these from RE model to FE model w/ cred int
fixed.df <- data.frame(param = row.names(countFit.out), countFit.out, row.names = NULL)
params <- strsplit(as.character(fixed.df$param), "\\[|\\]")
params.df <- data.frame(matrix(unlist(params[1:68]), nrow = 68, byrow = T))
fixedbetas <- data.frame(param = params.df$X1, strata = params.df$X2, mean = exp(fixed.df$mean[1:68]), 
                         cilow = exp(fixed.df$X2.5.[1:68]), cihi = exp(fixed.df$X97.5.[1:68]))

paramnames <- c("b1" = "Intercept", "b2" = "Slope")
fixplot <- ggplot(fixedbetas, aes(x = as.numeric(strata), y = mean)) + geom_point() + ylim(0.9, 1.3) + 
  geom_errorbar(aes(ymin = cilow, ymax = cihi)) + ylab("Effect size") + xlab(" ") + 
  facet_grid(~param, labeller = as_labeller(paramnames))

random.df <- data.frame(param = row.names(countRE$BUGSoutput$summary), countRE$BUGSoutput$summary, row.names = NULL)
paramre <- strsplit(as.character(random.df$param), "\\[|\\]")
paramre.df <- data.frame(matrix(unlist(paramre[1:1288]), nrow = 1288, byrow = T))
paramre.betas <- filter(paramre.df, X1 == "b1" | X1 == "b2")
re.betas12 <- cbind(random.df[1:1288, ], paramre.df)
randombetas <- filter(re.betas12, X1 == "b1" | X1 == "b2")

replot <- ggplot(randombetas, aes(x = as.numeric(as.character(X2)), y = exp(mean))) + geom_point() + ylim(0.9, 1.3) +
  geom_errorbar(aes(ymin = exp(X2.5.), ymax = exp(X97.5.))) + ylab("Effect size") + xlab("Strata") + 
  facet_grid(~X1, labeller = as_labeller(paramnames))

plot_grid(fixplot, replot, nrow = 2, labels = c("(a)", "(b)"))

# random effects: observer
# 1220 obsroute, 747 observers
alphas <- exp(alpha)
alph.df <- data.frame(alpha = alphas)

#histogram of random effects
a1 <- ggplot(alph.df, aes(x = alpha)) + geom_histogram(binwidth = 0.05, col = "white") + ylab("Number of observers") + 
  xlab("Random effect")
a1
# 502 routes, 3790 observations (route x year)

#point plot of random effects
alphapoints <- ggplot(random.df[1:1220, ], aes(x = param, y = exp(mean))) + geom_point() + 
  geom_errorbar(aes(ymin = exp(X2.5.), ymax = exp(X97.5.))) + xlab("Observer x Route ID") + ylab("Effect size") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
alphapoints

#years of observations per observer
obsyrs <- counts.merge %>% group_by(obsroute) %>% distinct(year) %>% summarize(nyears = n())

ggplot(obsyrs, aes(x = nyears)) + geom_histogram(binwidth = 1, col = "white") +
  xlab("Number of years") + ylab("Frequency")

# First year observer 
# for which species is there a nonzero overlapping CI 
setwd("C:/Users/gdicecco/Desktop/git/665_proj/model_output/")
files <- list.files()

library(stringr)
sums <- files[str_detect(files, "fixed_")]

firsts <- data.frame(X = NA, mean = NA, X2.5. = NA, X97.5. = NA, aou = NA, model = NA)
for(i in 1:length(sums)){
  df <- read.csv(sums[i])
  char <- str_split(sums[i], "_")
  subs <- df[df$X == "b3", c("X", "mean", "X2.5.", "X97.5.")]
  subs$aou <- char[[1]][1]
  subs$model <- char[[1]][3]
  firsts <- rbind(firsts, subs)
}

firsts <- na.omit(firsts)
firsts$model <- factor(firsts$model, levels = c("jagssummary.csv", "jagsRE.csv"))

modnames <- c("jagssummary.csv" = "Fixed effects", "jagsRE.csv" = "Hierarchical")
ggplot(firsts, aes(x = aou, y = exp(mean))) + geom_point() + ylab("First year observer effect size") + xlab("AOU") + 
  geom_errorbar(aes(ymin = exp(X2.5.), ymax = exp(X97.5.))) +
  facet_grid(~model, labeller = as_labeller(modnames)) + theme_bw()

# Plot empirical and predicted abundance on a map

data <- read.csv("counts_w_modeloutput.csv")

# Greater roadrunner

roadrun <- data %>% 
  filter(aou == 3850) %>%
  group_by(latitude, longitude, stateroute) %>%
  summarize(meancount = mean(speciestotal), meanre = mean(abundind), meanfixed = mean(fixedabundind))

par(mfrow = c(1,1))
longs = c(-125,-85)
lats = c(26,40)

# Raw counts
rbPal <- colorRampPalette(c("blue","red"))
roadrun$col <- rbPal(10)[as.numeric(cut(roadrun$meancount, breaks = 10))]

plot(bcrshp[bcrshp$WATER == 3 & bcrshp$COUNTRY == "USA",], xlim = longs, ylim = lats, border = "gray73", col = "gray95")
points(roadrun$longitude, roadrun$latitude, pch = 16, col = roadrun$col)
legend("bottomleft",title="Abund.",legend=c(1:10),col =rbPal(10),pch=20)

# Fixed effects
rbPal <- colorRampPalette(c("blue","red"))
roadrun$col <- rbPal(10)[as.numeric(cut(roadrun$meanfixed, breaks = 10))]

plot(bcrshp[bcrshp$WATER == 3 & bcrshp$COUNTRY == "USA",], xlim = longs, ylim = lats, border = "gray73", col = "gray95")
points(roadrun$longitude, roadrun$latitude, pch = 16, col = roadrun$col)
legend("bottomleft",title="Abund.",legend=c(1:10),col =rbPal(10),pch=20)

# Random effects 

rbPal <- colorRampPalette(c("blue","red"))
roadrun$col <- rbPal(10)[as.numeric(cut(roadrun$meanre, breaks = 10))]

plot(bcrshp[bcrshp$WATER == 3 & bcrshp$COUNTRY == "USA",], xlim = longs, ylim = lats, border = "gray73", col = "gray95")
points(roadrun$longitude, roadrun$latitude, pch = 16, col = roadrun$col)
legend("bottomleft",title="Abund.",legend=c(1:10),col =rbPal(10),pch=20)
