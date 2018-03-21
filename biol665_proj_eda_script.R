## BIOL 665 project EDA
  
library(dplyr)
library(geosphere)
library(maps)
library(rgdal)
library(rgeos)


#Read in BBS data
##Note: on mac once connected to Bioark server, path is "/Volumes/hurlbertlab/Databases/BBS/2017/filename.csv"
routes <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_routes_20170712.csv")
counts <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_counts_20170712.csv")
species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_species_20170712.csv")
weather <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\Databases\\BBS\\2017\\bbs_weather_20170712.csv")
bcrshp <- readOGR("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\bcr_terrestrial_shape\\BCR_Terrestrial_master.shp") #BCRs
##Typo in bcrshp - has TENNESSE

#species used in Huang 2017 GCB
huang_species <- read.csv("\\\\Bioark.bio.unc.edu\\hurlbertlab\\DiCecco\\BBS-centroids\\huang-2017-bbs-species.csv", header = TRUE)

#Remove routes weather runtype = 0
routes$stateroute <- routes$statenum*1000 + routes$route
weather$stateroute <-weather$statenum*1000 + weather$route
RT1 <- subset(weather, runtype == 1, select = c("stateroute", "year"))
RT1.routes <- merge(RT1, routes[ , c("statenum", "stateroute", "latitude", "longitude","bcr")], by = "stateroute", all.x = TRUE)
routes.short <- subset(RT1.routes, year >= 1969, select = c("statenum","stateroute", "year", "latitude", "longitude", "bcr"))
counts$stateroute <- counts$statenum*1000 + counts$route

longs2 = c(-175,-60)
lats2 = c(24,69)

#all routes from 1969 - 2016, runtype = 1
map(database = "world", xlim = longs2, ylim = lats2)
map(database = "state", add = TRUE)
points(routes.short$longitude, routes.short$latitude, pch = 16, cex = 0.1)

#all routes in strata with strata centers
#Calculate BCR/state strata centroids
regioncodes <- read.table("C:/Users/gdicecco/Desktop/git/bbs-centroid/regioncodes.txt", sep = "\t")
colnames(regioncodes) <- c("countrynum", "statenum", "PROVINCE_S")

polys.df <- data.frame(BCR = c(0), BCRNAME = c(0), PROVINCE_S = c(0), COUNTRY = c(0), REGION = c(0),
                       WATER = c(0), Shape_Leng = c(0), Shape_Area = c(0), x = c(0), y = c(0))
for(state in regioncodes$PROVINCE_S) {
  shape <- subset(bcrshp, PROVINCE_S == state)
  centers <- gCentroid(shape, byid = TRUE)
  df.temp <- data.frame(cbind(shape@data, centers@coords))
  polys.df <- rbind(polys.df, df.temp)
}
polys.df <- polys.df[-1,] #remove zero row
polys.merged <- merge(regioncodes, polys.df, by = "PROVINCE_S", all.x = TRUE) #merge PROVINCE_S with statenum
polys.merged.land <- subset(polys.merged, WATER == 3) #remove centroids for bodies of water

counts.short <- counts %>%
  filter(year >= 1969) %>%
  filter(aou %in% huang_species$aou) %>%
  select(year, aou, speciestotal, stateroute)

#merge strata centroids with routes
routes.short$statebcr <- routes.short$statenum*1000 + routes.short$bcr
polys.merged.land$statebcr <- polys.merged.land$statenum*1000 + polys.merged.land$BCR
routes.short.centers <- merge(routes.short, polys.merged.land, by = "statebcr")
counts.merged.centers <- merge(routes.short.centers, counts.short, by = c("stateroute","year"))

plot(bcrshp[bcrshp$WATER == 3,], ylim = lats2, xlim = longs2, border = "white", col = "white") #plot centroids on BCR map
points(counts.merged.centers$longitude, counts.merged.centers$latitude, cex = 0.1, pch = 16, col = "gray40")
plot(bcrshp[bcrshp$WATER == 3,], ylim = lats2, xlim = longs2, border = "black", add = TRUE)
points(counts.merged.centers$x, counts.merged.centers$y, cex = 0.5, pch = 16, col = "red")
legend(-174,40, legend = c("Routes","Geographic centers"), pch = 16, col = c("gray40", "red"), bty= "n")

#### Plots about observers?
obs <- weather %>%
  group_by(obsn) %>%
  summarize(total = n(), mean.spp = mean(totalspp), sd.spp = sd(totalspp))
            
plot(log(obs$total), log(obs$mean.spp), xlab = "Log total routes run", ylab = "Log mean number of species identified per route")

exp <- obs %>%
  select(total, mean.spp) %>%
  group_by(total) %>%
  summarize(mean = mean(mean.spp))

barplot(height = exp$mean, names.arg = exp$total, xlab = "Number of routes surveyed", ylab = "Mean number of species seen per route")
cor(exp)

#observer experience vs. total birds
observers <- left_join(counts, weather, by = c("countrynum", "statenum", "route", "year"))
observers$routeyear <- paste(observers$countrynum, observers$statenum, observers$route, observers$year, sep = "")
abunds <- observers %>%
  group_by(obsn) %>%
  summarize(total = length(unique(routeyear)), abund = mean(stoptotal))
plot(log(abunds$total), log(abunds$abund), xlab = "Log total routes run", ylab = "Log number of individuals recorded per route")
cor(abunds$total, abunds$abund)
