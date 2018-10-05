# This script file includes the code necessary to construct a delta model for estimating spatial overlap between Pacific Halibut and Arrowtooth Flounder in the Gulf of Alaska. Methods modified from Hunsicker et al. (2013) and Shelton et al. (2017).

# We used standardized survey data procured from the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]). Bottom trawl survey data (1990 to 2017) were collected by the Resource Assessment and Conservation Engineering (RACE) Division and are publicly accessible here: https://www.afsc.noaa.gov/RACE/groundfish/survey_data/data.htm. 

# We used generalized additive models to quantify and predict the probability of occurrence and relative abundance of Pacific Halibut and Arrowtooth Flounder across a uniform grid spanning the study area. We multiplied probabilities of occurrence and relative abundances to estimate overall abundance in each survey year-grid cell combination. We then multiplied standardized abundances of each species to estimate spatial overlap in each combination of survey year and grid cell. 

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu

# References:
# Hunsicker ME, Ciannelli L, Bailey KM, Zador S, Stige L. Climate and demography dictate the strength of predator-prey overlap in a subarctic marine ecosystem. PLOS ONE. 2013;8(6):e66025. doi:10.1371/journal.pone.006602
# Livingston PA, Aydin K, Buckley TW, Lang GM, Yang MS, Miller BS. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes. 2017;100(4):443–470.
# Shelton AO, Hunsicker ME, Ward EJ, Feist BE, Blake R, Ward CL, et al. Spatio-temporal models reveal subtle changes to demersal communities following the Exxon Valdez oil spill. ICES J Mar Sci. 2017. doi: 10.1093/icesjms/fsx079
# von Szalay PG, Raring NW. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle (WA): National Oceanic and Atmospheric Administration; 2016. Technical Memorandum: NMFS-AFSC-325. Sponsored by the US Department of Commerce.

rm(list=ls())
graphics.off()

setwd("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/ResourcePartitioning/ResourcePartitioning")
################################################################
### DATA PREPARATION

# Format trawl data:
trawl = read.csv("AFSC_TrawlData_1984_2017.csv") # These data include both positive and zero CPUE values!
trawl$STATIONID = as.factor(trawl$STATIONID)
trawl$SPECIES_CODE = as.factor(trawl$SPECIES_CODE)

trawl$Species = trawl$SPECIES_CODE
levels(trawl$Species)
levels(trawl$Species) = list(ATF="10110", PH="10120", SBL="20510", PC="21720", WEP="21740")

# Assign Management Areas (values) based upon Survey Stratum (index):
index = c(unique(trawl$STRATUM))

values = c("650", "650", "650", "650", "650", "640", "640", "640", "640", "640", "640", "640", "640", "620", "640", "620",
           "640", "610", "650", "610", "650", "610", "610", "610", "610", "610", "610", "610", "610", "610", "620", "620",
           "620", "620", "620", "630", "630", "630", "630", "630", "630", "630", "630", "640", "640", "630", "630", "630", 
           "630", "630", "620", "620", "620", "620", "630", "630", "630", "650", "650") 

trawl$MgmtArea = values[match(trawl$STRATUM, index)]
table(trawl$STRATUM, trawl$MgmtArea)

# Exclude 1984 and 1987 data:
trawl = subset(trawl, YEAR!=1984)
trawl = subset(trawl, YEAR!=1987)

ATF = subset(trawl, Species=="ATF")
PH = subset(trawl, Species=="PH")
PH_ATF = rbind(PH, ATF)

require(dplyr)
meanCPUE = PH_ATF %>% 
  group_by(Species, YEAR) %>% 
  summarise(aveCatch = mean(NUMCPUE))

require(ggplot2)
PH_ATF_NUMCPUE = ggplot(data=meanCPUE, aes(x=YEAR, y=aveCatch, group=Species)) +
  geom_line() +
  facet_wrap(~ Species, scale="free_y") +
  scale_color_manual(values=c("blue", "red")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.0375, 0.945)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=12)) +
  theme(axis.text.x = element_text(family="Arial", size=11)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=12)) +
  theme(axis.title.y = element_text(vjust=1.1, size=12)) +
  theme(strip.text.x = element_text(family="Arial", size=12)) +
  theme(strip.text = element_text(hjust=0.5)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) 

ggsave(filename="PH_ATF_NUMCPUE.png", plot=PH_ATF_NUMCPUE, dpi=500, width=12, height=8, units="in")

### Set depth bins (100 m incremements: < 100, 100-199, 200-299, 300-399, 400-499, >=500):
trawl$DepthBin = cut(trawl$GEAR_DEPTH, breaks = c(0, 99, 199, 299, 399, 499, 600))
levels(trawl$DepthBin) = c("<100", "100-199", "200-299", "300-399", "400-499", ">=500")
# Check bins:
table(trawl$NUMBER_FISH, trawl$DepthBin)

require(reshape2)
trawl_wide_all = dcast(trawl, YEAR + HAULJOIN + VESSEL + CRUISE + STRATUM + DISTANCE_FISHED + NET_WIDTH + STATIONID + START_LATITUDE + START_LONGITUDE + GEAR_DEPTH + GEAR_TEMPERATURE + MgmtArea ~ Species, value.var="NUMCPUE")

# Sum NUMBER_FISH by species, management area, and year:

# Exclude incomplete entries not used in modeling (i.e., rows with missing depth or BT data):
trawl_comp1 = subset(trawl, GEAR_DEPTH >=0)
trawl_comp = subset(trawl_comp1, GEAR_TEMPERATURE >=0)

trawlz = trawl[, c(2,4,20,21,23)]
trawlz = unique(trawlz)

depthN1 = trawlz %>% 
  group_by(MgmtArea, YEAR) %>% 
  summarise(CountDepth1 = length(GEAR_DEPTH))
depthN1 = as.data.frame(depthN1)

depthN2 = trawlz %>% 
  group_by(MgmtArea, YEAR) %>% 
  summarise(CountDepth2 = sum(!is.na(GEAR_DEPTH)))
depthN2 = as.data.frame(depthN2)

depthN = cbind(depthN1, depthN2)
depthN$propNA = 1 - (depthN$CountDepth2 / depthN$CountDepth1)

require(tidyr)
trawl_red0 = trawl[, c(2,4,20,21,23)]
trawl_red0 = unique(trawl_red0)
trawl_red1 = trawl_red0 %>% drop_na(GEAR_DEPTH)
trawl_red2 = trawl_red0 %>% drop_na(GEAR_TEMPERATURE)
trawl_red3 = trawl_red0 %>% drop_na(GEAR_DEPTH, GEAR_TEMPERATURE)

table(trawl_red0$MgmtArea, trawl_red0$YEAR)
table(trawl_red3$MgmtArea, trawl_red3$YEAR)

tempN1 = trawlz %>% 
  group_by(MgmtArea, YEAR) %>% 
  summarise(CountTemp1 = length(GEAR_TEMPERATURE))
tempN1 = as.data.frame(tempN1)

tempN2 = trawlz %>% 
  group_by(MgmtArea, YEAR) %>% 
  summarise(CountTemp2 = sum(!is.na(GEAR_TEMPERATURE)))
tempN2 = as.data.frame(tempN2)

tempN = cbind(tempN1, tempN2)
tempN$propNA = 1 - (tempN$CountTemp2 / tempN$CountTemp1)

# Treat year as a factor for all analyses:
trawl_comp$YEAR = as.factor(trawl_comp$YEAR)

# Create own "Haul_Join" column by concatenating Vessel, Cruise, and Haul:
trawl_comp$Haul_Join = paste(trawl_comp$VESSEL, trawl_comp$CRUISE, trawl_comp$HAUL, sep="")

cor.test(trawl_red0$GEAR_DEPTH, trawl_red0$GEAR_TEMPERATURE)
################################################################
# Subset data by species:
PH = subset(trawl_comp, Species=="PH")
ATF = subset(trawl_comp, Species=="ATF")
################################################################
# Convert data from long to wide format:
trawl_wide = dcast(trawl_comp, YEAR + HAULJOIN + VESSEL + CRUISE + Haul_Join + STRATUM + DISTANCE_FISHED + NET_WIDTH + STATIONID + START_LATITUDE + START_LONGITUDE + GEAR_DEPTH + GEAR_TEMPERATURE + MgmtArea ~ Species, value.var="NUMCPUE")

# Set NA values to zero for when a particular species was not caught during an individual haul:
trawl_wide[is.na(trawl_wide)] = 0

# Rename variables to match previously written code:
colnames(trawl_wide) = c("year", "HAULJOIN", "Vessel", "Cruise", "Haul_Join", "Stratum", "Distance_Fished", "Net_Width", "Station", "lat", "lon", "depth", "BT", "MgmtArea", "ATF", "PH", "SBL", "PC", "WEP")

# Calculate number of tows conducted:
No_tows_YR = trawl_wide %>% 
  group_by(MgmtArea, year) %>% 
  summarise(CountTows = length(Haul_Join))
# print(No_tows_YR, n=400)

trawl_wide %>% 
  group_by(year) %>% 
  summarise(CountTows = length(Haul_Join))

# Calculate number of positive tows by species:
PH_n = subset(trawl_wide, PH > 0)
ATF_n = subset(trawl_wide, ATF > 0)

# Calculate number of positive tows by species:
PH_tows_YR = PH_n %>% 
  group_by(year) %>% 
  summarise(CountTows = length(Haul_Join))

ATF_tows_YR = ATF_n %>% 
  group_by(year) %>% 
  summarise(CountTows = length(Haul_Join))

PH_tows = PH_n %>% 
  group_by(MgmtArea, year) %>% 
  summarise(CountTows = length(Haul_Join))

ATF_tows = ATF_n %>% 
  group_by(MgmtArea, year) %>% 
  summarise(CountTows = length(Haul_Join))
################################################################
# Input proportional length information (30-69 cm fish):
PHpropLength = read.csv("~/Documents/UAF/Dissertation/Analyses/Bartolino_Methods/HaulData_PH_30_69cm.csv")

# Set NA values to zero:
PHpropLength[is.na(PHpropLength)] = 0

# Remove unnecessary columns:
PHpropLength = PHpropLength[,c(2,4:6)]

ATFpropLength = read.csv("~/Documents/UAF/Dissertation/Analyses/Bartolino_Methods/HaulData_ATF_30_69cm.csv")

# Set NA values to zero:
ATFpropLength[is.na(ATFpropLength)] = 0

# Remove unnecessary columns:
ATFpropLength = ATFpropLength[,c(2,4:6)]

# Remove other species from data frame:
trawl_wide = trawl_wide[,1:16]

# Join proportional length data frame and other trawl data:
trawl_wide$Haul_Join = as.character(trawl_wide$Haul_Join)
PHpropLength$Haul_Join = as.character(PHpropLength$Haul_Join)
ATFpropLength$Haul_Join = as.character(ATFpropLength$Haul_Join)

trawl_wide_30_69 = trawl_wide %>% left_join(ATFpropLength)
trawl_wide_30_69 = trawl_wide_30_69 %>% left_join(PHpropLength)

# Set NA values to zero:
trawl_wide_30_69[is.na(trawl_wide_30_69)] = 0
##################################################################
# Calculate new CPUE (by number) based on proportion of catch that measured 30 to 69 cm fork length!
trawl_wide_30_69$CPUE_ATF = trawl_wide_30_69$ATF * trawl_wide_30_69$ATFToIncl

trawl_wide_30_69$CPUE_PH = trawl_wide_30_69$PH * trawl_wide_30_69$PHToIncl
##################################################################
### MODEL FITTING
require(mgcv)
require(lme4)
require(MuMIn)
options(na.action = "na.fail") 

# Model P/A for Pacific Halibut:
trawl_wide_30_69$PHpa = as.numeric(trawl_wide_30_69$CPUE_PH > 0)
nrow(trawl_wide_30_69)

# Remove super deep PH:
trawl_wide_30_69sub = subset(trawl_wide_30_69, Haul_Join!="148201101201")
nrow(trawl_wide_30_69sub)
table(trawl_wide_30_69sub$year)
table(trawl_wide_30_69sub$year)

# Sample sizes:
trawl_wide_30_69subPHpa = subset(trawl_wide_30_69sub, PHpa=="1")
nrow(trawl_wide_30_69subPHpa)
table(trawl_wide_30_69subPHpa$year)
table(trawl_wide_30_69subPHpa$MgmtArea)
table(trawl_wide_30_69subPHpa$MgmtArea, trawl_wide_30_69subPHpa$year)

levels(trawl_wide_30_69sub$year)
# levels(trawl_wide_30_69sub$year) = c("1990", "93", "96", "99", "01", "03", "05", "07", "09", "11", "13", "15", "2017")

# Full (global) model:
PH.pa.gam_full = gam(PHpa ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = trawl_wide_30_69sub, family = binomial(link=logit), method="GCV.Cp")
summary(PH.pa.gam_full)

# Generate all possible alternative models and select best-fit based on delta AIC (year not as factor = necessary):
PH.pa.gam_select = dredge(PH.pa.gam_full, beta=FALSE, evaluate=TRUE, rank="AIC", trace=FALSE)

print(PH.pa.gam_select, abbrev.names=FALSE, warnings=TRUE)
summary(PH.pa.gam_select)

PH.pa.gam_1 = gam(PHpa ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = trawl_wide_30_69sub, family = binomial(link=logit), method="GCV.Cp")
summary(PH.pa.gam_1)
PH.pa.gam_2 = gam(PHpa ~ year + s(lon, lat) + s(depth), data = trawl_wide_30_69sub, family = binomial(link=logit), method="GCV.Cp")
summary(PH.pa.gam_2)
PH.pa.gam_3 = gam(PHpa ~ s(lon, lat) + s(depth) + s(BT, k=4), data = trawl_wide_30_69sub, family = binomial(link=logit), method="GCV.Cp")
summary(PH.pa.gam_3)

PH.pa.gam_best = PH.pa.gam_full
summary(PH.pa.gam_best)
gam.check(PH.pa.gam_best)
gam.fit(PH.pa.gam_best)
# plot(PH.pa.gam_best, scale=0, lwd=2)
# gam.check(PH.pa.gam_best)
# PH.pa.gam_best_fit = fitted(PH.pa.gam_best)

# Variogram for assessment of spatial autocorrelation:
require(geoR)
pred<-PH.pa.gam_best$fitted

#Convert coordinates from Lat & Lon to UTM (km)
subdataLL<-trawl_wide_30_69[,c("lon","lat")]
names(subdataLL)<-c('X','Y')
attr(subdataLL,'projection')<-'LL'
subdataUTM<-convUL(subdataLL)

subdataUTM$resd<-trawl_wide_30_69$PHpa-pred
subdataUTM$PHpa<-trawl_wide_30_69$PHpa

#Estimate variogram of data and residuals using geoR library
subset.geor.PHpa<-as.geodata(subdataUTM, coords.col=c(1,2),data.col=4)
subset.geor.res<-as.geodata(subdataUTM, coords.col=c(1,2),data.col=3)
vario.PHpa<-variog(subset.geor.PHpa,uvec=seq(20,500,by=25),pairs.min=10)

vario.res<-variog(subset.geor.res,uvec=seq(20,500,by=25),pairs.min=10)

plot(vario.PHpa,main='Variogram',xlab='Distance
     (km)',ylab='Semivariance',type='b',scaled=T, ylim=c(0.5,1.5))
lines(vario.res,lty=2,scaled=T, col="red")
abline(v=100)
# abline(h=.88)
legend(200,0.7,c('PHpa','residuals'),lty=c(1,2))
#################################################################
require(sp)

# Plot GAM results:
require(maps)
require(mapdata)
require(visreg)
require(ggplot2)
require(PBSmapping)
require(extrafont)

setEPS()
postscript("some_graph.eps")

lonmin = -172
lonmax = -130
latmin = 52
latmax = 62

data(worldHiresMapEnv)

# Use vis.gam for lat lon data:
pdf("vis.gam_PHpa_30_69.pdf", width=10, height=6)
vis.gam(PH.pa.gam_best, c("lon", "lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Pacific Halibut, Partial Effect on Presence or Absence", too.far=0.025, n.grid=250, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))
maps::map('worldHires', fill=T, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax), add=T, col="lightgrey")
dev.off()

PH.pa.gam_best$gcv.ubre

PHpaYRgam = visreg(fit=PH.pa.gam_best, xvar="year", band=TRUE, partial=FALSE, rug=FALSE, line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), trans=binomial()$linkinv, type="conditional", gg=TRUE, labels=c("1990", "", "96", "", "01", "", "05", "", "09", "", "13", "", "2017")) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on Presence (1) or Absence (0)") +
  scale_y_continuous(limits=c(0.30,1.0), breaks=c(0.30, 0.4, 0.50, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme(legend.background = element_rect(fill="transparent")) 
PHpaYRgam

ggsave(filename="PHpaYRgam_30_69.png", plot=PHpaYRgam, dpi=500, width=12, height=8, units="in")

PHpaDepthGAM = visreg(PH.pa.gam_best, xvar="depth",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on Presence (1) or Absence (0)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000), breaks=c(0,150,300,450,600,750,900)) +
  theme(legend.background = element_rect(fill="transparent")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
PHpaDepthGAM

ggsave(filename="PHpaDepthGAM_30_69.png", plot=PHpaDepthGAM, dpi=500, width=12, height=8, units="in")

PHpaBTgam = visreg(PH.pa.gam_best, xvar="BT",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x=expression(paste("Bottom Temperature (",degree,"C)")), y="Partial Effect on Presence (1) or Absence (0)") +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25, 0.50, 0.75, 1)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
PHpaBTgam

ggsave(filename="PHpaBTgam_30_69.png", plot=PHpaBTgam, dpi=500, width=12, height=8, units="in")

# Can visualize by multiple factors at once (but won't):
# visreg2d(PH.pa.gam_full, "year", "BT", plot.type="image")

plot(trawl_wide_30_69$year, resid(PH.pa.gam_full)); abline(h=0, lty=2)
plot(trawl_wide_30_69$depth, resid(PH.pa.gam_full)); abline(h=0, lty=2)
plot(trawl_wide_30_69$BT, resid(PH.pa.gam_full)); abline(h=0, lty=2)
plot(fitted(PH.pa.gam_full),resid(PH.pa.gam_full), # Evidence of heteroscedasticity?
     xlab="Fitted values", ylab="Residuals")
abline(h=0,lty=2); grid()
################################################################
# Model CPUE (where present) for Pacific Halibut:
PH = subset(trawl_wide_30_69, CPUE_PH > 0)
PH = PH[, c(1:14, 16, 20:22, 24:25)]
PH$logPH = log(PH$CPUE_PH)
max(PH$depth)
nrow(PH)
table(PH$Haul_Join, PH$year)

# Exclude extreme outlier by depth:
PH = subset(PH, Haul_Join!="148201101201")
max(PH$depth)
nrow(PH)
table(PH$year)

# Full (global) model:
PH.cpue.gam_full = gam(logPH ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = PH, family = gaussian(link=identity), method="GCV.Cp")
summary(PH.cpue.gam_full)

# Generate all possible alternative models and select best-fit based on delta AIC (GRAYED OUT TO AVOID RE-RUNNING):
PH.cpue.gam_select = dredge(PH.cpue.gam_full, beta=FALSE, evaluate=TRUE, rank="AIC", trace=FALSE)

print(PH.cpue.gam_select, abbrev.names=FALSE, warnings=TRUE)

PH.cpue.gam_1 = gam(logPH ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = PH, family = gaussian(link=identity), method="GCV.Cp")
summary(PH.cpue.gam_1)
PH.cpue.gam_2 = gam(logPH ~ year + s(lon, lat) + s(depth), data = PH, family = gaussian(link=identity), method="GCV.Cp")
summary(PH.cpue.gam_2)
PH.cpue.gam_3 = gam(logPH ~ s(lon, lat) + s(depth) + s(BT, k=4), data = PH, family = gaussian(link=identity), method="GCV.Cp")
summary(PH.cpue.gam_3)

PH.cpue.gam_best = PH.cpue.gam_full
summary(PH.cpue.gam_best)
# plot(PH.cpue.gam_best, scale=0, lwd=2)
# gam.check(PH.cpue.gam_best)
# PH.cpue.gam_best_fit = fitted(PH.cpue.gam_best)

# Variogram for assessment of spatial autocorrelation:
pred2<-PH.cpue.gam_best$fitted

#Convert coordinates from Lat & Lon to UTM (km)
subdataLL2<-PH[,c("lon","lat")]
names(subdataLL2)<-c('X','Y')
attr(subdataLL2,'projection')<-'LL'
subdataUTM2<-convUL(subdataLL2)

subdataUTM2$resd<-PH$logPH-pred2
subdataUTM2$PHcpue<-PH$logPH

#Estimate variogram of data and residuals using geoR library
subset.geor.PHcpue<-as.geodata(subdataUTM2, coords.col=c(1,2),data.col=4)
subset.geor.res2<-as.geodata(subdataUTM2, coords.col=c(1,2),data.col=3)
vario.PHcpue<-variog(subset.geor.PHcpue,uvec=seq(20,500,by=25),pairs.min=10)

vario.res2<-variog(subset.geor.res2,uvec=seq(20,500,by=25),pairs.min=10)

plot(vario.PHcpue,main='Variogram',xlab='Distance
     (km)',ylab='Semivariance',type='b',scaled=T, ylim=c(0.5,1.5))
lines(vario.res2,lty=2,scaled=T, col="red")
abline(v=100)
# abline(h=1.06)
legend(200,1.34,c('PHcpue','residuals'),lty=c(1,2))
#################################################################
# Use vis.gam for lat lon data:
pdf("vis.gam_PHcpue_30_69.pdf", width=10, height=6)
vis.gam(PH.cpue.gam_best, c("lon", "lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Pacific Halibut, Partial Effect on log-CPUE (number per hectare)", too.far=0.025, n.grid=500, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))
maps::map('worldHires', fill=T, add=T, col="lightgrey")
dev.off()

# visreg provides better visualizations than plot(gam), but it splits the spatial smooth into separate fits for longitudinal effects (at mean latitude?) and latitudinal effects (at mean longitude?). Also, I'm not sure how the confidence intervals are constructed (F. Mueter)!
PHcpueYRgam = visreg(PH.cpue.gam_best, xvar="year", line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), band=TRUE, partial=FALSE, rug=FALSE, trans=gaussian()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(size=18, hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="", y="Partial Effect on log-CPUE (number per hectare)") +
  scale_y_continuous(limits=c(4.5, 6.0), breaks=c(4.5, 5.0, 5.5, 6.0)) +
  theme(legend.background = element_rect(fill="transparent"))
PHcpueYRgam

ggsave(filename="PHcpueYRgam_30_69.png", plot=PHcpueYRgam, dpi=500, width=12, height=8, units="in")

PHcpueDepthGAM = visreg(fit=PH.cpue.gam_best, xvar="depth",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on log-CPUE (number per hectare)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000), breaks=c(0,150,300,450,600,750,900)) +
  scale_y_continuous(limits=c(2,8), breaks=c(2,4,6,8)) +
  theme(legend.background = element_rect(fill="transparent")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
PHcpueDepthGAM

ggsave(filename="PHcpueDepthGAM_30_69.png", plot=PHcpueDepthGAM, dpi=500, width=12, height=8, units="in")

PHcpueBTgam = visreg(PH.cpue.gam_best, xvar="BT",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Pacific Halibut") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  scale_x_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +
  scale_y_continuous(limits=c(4.5,7), breaks=c(5,6,7)) +
  labs(x=expression(paste("Bottom Temperature (",degree,"C)")), y="Partial Effect on log-CPUE (number per hectare)") +
  theme(legend.background = element_rect(fill="transparent")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
PHcpueBTgam

ggsave(filename="PHcpueBTgam_30_69.png", plot=PHcpueBTgam, dpi=500, width=12, height=8, units="in")

# Can visualize by multiple factors at once (but won't):
# visreg2d(PH.pa.gam_best, "year", "BT", plot.type="image")

# Diagnostic Plots:
plot(PH$year, resid(PH.cpue.gam_best)); abline(h=0, lty=2)
plot(PH$depth, resid(PH.cpue.gam_best)); abline(h=0, lty=2)
# plot(trawl_wide$BT, resid(PH.cpue.gam_best)); abline(h=0, lty=2)
plot(fitted(PH.cpue.gam_best), resid(PH.cpue.gam_best), # Evidence of heteroscedasticity?
     xlab="Fitted values", ylab="Residuals")
abline(h=0,lty=2); grid()

# Residuals by depth and over time:
r <- resid(PH.cpue.gam_best)
j <- r > 0
symbols(PH$year[j], PH$depth[j], circles = r[j], bg = 4, inches=0.05) 
symbols(PH$year[!j], PH$depth[!j], circles = abs(r[!j]), bg = 2, inches=0.05, add=T) 

# Residuals by BT and over time:
symbols(PH$year[j], PH$BT[j], circles = r[j], bg = 4, inches=0.05) 
symbols(PH$year[!j], PH$BT[!j], circles = abs(r[!j]), bg = 2, inches=0.05, add=T) 
################################################################
################################################################
### MODEL FITTING

# Model P/A for Arrowtooth Flounder:
trawl_wide_30_69$ATFpa = as.numeric(trawl_wide_30_69$CPUE_ATF > 0)

# Sample sizes:
trawl_wide_30_69ATFpa = subset(trawl_wide_30_69, ATFpa=="1")
nrow(trawl_wide_30_69ATFpa)
table(trawl_wide_30_69ATFpa$year)
table(trawl_wide_30_69ATFpa$MgmtArea)
table(trawl_wide_30_69ATFpa$MgmtArea, trawl_wide_30_69ATFpa$year)

# Full (global) model:
ATF.pa.gam_full = gam(ATFpa ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = trawl_wide_30_69, family = binomial(link=logit), method="GCV.Cp")
summary(ATF.pa.gam_full)

# Generate all possible alternative models and select best-fit based on delta AIC (GRAYED OUT TO AVOID RE-RUNNING):
ATF.pa.gam_select = dredge(ATF.pa.gam_full, beta=FALSE, evaluate=TRUE, rank="AIC", trace=FALSE)

print(ATF.pa.gam_select, abbrev.names=FALSE, warnings=TRUE)

ATF.pa.gam_1 = gam(ATFpa ~ year + s(lon, lat) + s(depth), data = trawl_wide_30_69, family = binomial(link=logit), method="GCV.Cp")
summary(ATF.pa.gam_1)
ATF.pa.gam_2 = gam(ATFpa ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = trawl_wide_30_69, family = binomial(link=logit), method="GCV.Cp")
summary(ATF.pa.gam_2)
ATF.pa.gam_3 = gam(ATFpa ~ s(lon, lat) + s(depth) + s(BT, k=4), data = trawl_wide_30_69, family = binomial(link=logit), method="GCV.Cp")
summary(ATF.pa.gam_3)

ATF.pa.gam_best = ATF.pa.gam_full
summary(ATF.pa.gam_best)
# plot(ATF.pa.gam_best, scale=0, lwd=2)
# gam.check(ATF.pa.gam_best)
# ATF.pa.gam_best_fit = fitted(ATF.pa.gam_best)

# Variogram for assessment of spatial autocorrelation:
pred3<-ATF.pa.gam_best$fitted

#Convert coordinates from Lat & Lon to UTM (km)
subdataLL3<-trawl_wide_30_69[,c("lon","lat")]
names(subdataLL3)<-c('X','Y')
attr(subdataLL3,'projection')<-'LL'
subdataUTM3<-convUL(subdataLL3)

subdataUTM3$resd<-trawl_wide_30_69$ATFpa-pred3
subdataUTM3$ATFpa<-trawl_wide_30_69$ATFpa

#Estimate variogram of data and residuals using geoR library
subset.geor.ATFpa<-as.geodata(subdataUTM3, coords.col=c(1,2),data.col=4)
subset.geor.res3<-as.geodata(subdataUTM3, coords.col=c(1,2),data.col=3)
vario.ATFpa<-variog(subset.geor.ATFpa,uvec=seq(20,500,by=25),pairs.min=10)

vario.res3<-variog(subset.geor.res3,uvec=seq(20,500,by=25),pairs.min=10)

plot(vario.ATFpa,main='Variogram',xlab='Distance
     (km)',ylab='Semivariance',type='b',scaled=T, ylim=c(0.5,1.5))
lines(vario.res3,lty=2,scaled=T, col="red")
abline(v=100)
# abline(h=0.98)
legend(200,1.34,c('ATFpa','residuals'),lty=c(1,2))
################################################################
# Plot GAM results:
table(trawl_wide_30_69$year)

# Use vis.gam for lat lon data:
pdf("vis.gam_ATFpa_30_69.pdf", width=10, height=6)
vis.gam(ATF.pa.gam_full, c("lon", "lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Arrowtooth Flounder, Partial Effect on Presence or Absence", too.far=0.025, n.grid=500, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))
maps::map('worldHires', fill=T, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax), add=T, col="lightgrey")
dev.off()

ATFpaYRgam = visreg(fit=ATF.pa.gam_full, xvar="year", band=TRUE, partial=FALSE, rug=FALSE, line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), trans=binomial()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Arrowtooth Flounder") +
  theme(plot.title = element_text(size=18, hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="", y="Partial effect on Presence (1) or Absence (0)") +
  scale_y_continuous(limits=c(0.30,1.0), breaks=c(0.30, 0.4, 0.50, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme(legend.background = element_rect(fill="transparent"))
ATFpaYRgam

ggsave(filename="ATFpaYRgam_30_69.png", plot=ATFpaYRgam, dpi=500, width=12, height=8, units="in")

ATFpaDepthGAM = visreg(ATF.pa.gam_full, xvar="depth",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Arrowtooth Flounder") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on Presence (1) or Absence (0)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000), breaks=c(0,150,300,450,600,750,900)) +
  theme(legend.background = element_rect(fill="transparent")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ATFpaDepthGAM

ggsave(filename="ATFpaDepthGAM_30_69.png", plot=ATFpaDepthGAM, dpi=500, width=12, height=8, units="in")

ATFpaBTgam = visreg(ATF.pa.gam_full, xvar="BT",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=binomial()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Arrowtooth Flounder") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x=expression(paste("Bottom Temperature (",degree,"C)")), y="Partial Effect on Presence (1) or Absence (0)") +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.25, 0.50, 0.75, 1)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ATFpaBTgam

ggsave(filename="ATFpaBTgam_30_69.png", plot=ATFpaBTgam, dpi=500, width=12, height=8, units="in")

# Can visualize by multiple factors at once (but won't):
# visreg2d(ATF.pa.gam_best, "year", "BT", plot.type="image")

# Diagnostic Plots:
plot(trawl_wide_30_69$year, resid(ATF.pa.gam_full)); abline(h=0, lty=2)
plot(trawl_wide_30_69$depth, resid(ATF.pa.gam_full)); abline(h=0, lty=2)
plot(trawl_wide_30_69$BT, resid(ATF.pa.gam_full)); abline(h=0, lty=2)
plot(fitted(ATF.pa.gam_full),resid(ATF.pa.gam_full), # Evidence of heteroscedasticity?
     xlab="Fitted values", ylab="Residuals")
abline(h=0,lty=2); grid()

# Residuals by depth and over time:
r <- resid(ATF.pa.gam_best)
j <- r > 0
symbols(trawl_wide_30_69$year[j], trawl_wide_30_69$depth[j], circles = r[j], bg = 4, inches=0.05) 
symbols(trawl_wide_30_69$year[!j], trawl_wide_30_69$depth[!j], circles = abs(r[!j]), bg = 2, inches=0.05, add=T) 

################################################################
# Model CPUE (where present) for Arrowtooth Flounder:
ATF = subset(trawl_wide_30_69, CPUE_ATF > 0)
ATF = ATF[, c(1:15, 17:19, 23, 26)]
ATF$logATF = log(ATF$CPUE_ATF)

nrow(ATF)
table(ATF$year)

# Full (global) model:
ATF.cpue.gam_full = gam(logATF ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = ATF, family = gaussian(link=identity), method="GCV.Cp")
summary(ATF.cpue.gam_full)

# Generate all possible alternative models and select best-fit based on delta AIC (GRAYED OUT TO AVOID RE-RUNNING):
ATF.cpue.gam_select = dredge(ATF.cpue.gam_full, beta=FALSE, evaluate=TRUE, rank="AIC", trace=FALSE)

print(ATF.cpue.gam_select, abbrev.names=FALSE, warnings=TRUE)

ATF.cpue.gam_1 = gam(logATF ~ year + s(lon, lat) + s(depth) + s(BT, k=4), data = ATF, family = gaussian(link=identity), method="GCV.Cp")
summary(ATF.cpue.gam_1)
ATF.cpue.gam_2 = gam(logATF ~ year + s(lon, lat) + s(depth), data = ATF, family = gaussian(link=identity), method="GCV.Cp")
summary(ATF.cpue.gam_2)
ATF.cpue.gam_3 = gam(logATF ~ s(lon, lat) + s(depth) + s(BT, k=4), data = ATF, family = gaussian(link=identity), method="GCV.Cp")
summary(ATF.cpue.gam_3)

ATF.cpue.gam_best = ATF.cpue.gam_full
summary(ATF.cpue.gam_best)
# plot(ATF.pa.gam_best, scale=0, lwd=2)
# gam.check(ATF.cpue.gam_best)
# ATF.cpue.gam_best_fit = fitted(ATF.cpue.gam_best)

# Variogram for assessment of spatial autocorrelation:
pred4<-ATF.cpue.gam_best$fitted

#Convert coordinates from Lat & Lon to UTM (km)
subdataLL4<-ATF[,c("lon","lat")]
names(subdataLL4)<-c('X','Y')
attr(subdataLL4,'projection')<-'LL'
subdataUTM4<-convUL(subdataLL4)

subdataUTM4$resd<-ATF$logATF-pred4
subdataUTM4$ATFcpue<-ATF$logATF

#Estimate variogram of data and residuals using geoR library
subset.geor.ATFcpue<-as.geodata(subdataUTM4, coords.col=c(1,2),data.col=4)
subset.geor.res4<-as.geodata(subdataUTM4, coords.col=c(1,2),data.col=3)
vario.ATFcpue<-variog(subset.geor.ATFcpue,uvec=seq(20,500,by=25),pairs.min=10)

vario.res4<-variog(subset.geor.res4,uvec=seq(20,500,by=25),pairs.min=10)

plot(vario.ATFcpue,main='Variogram',xlab='Distance
     (km)',ylab='Semivariance',type='b',scaled=T, ylim=c(0.5,1.5))
lines(vario.res2,lty=2,scaled=T, col="red")
abline(v=100)
# abline(h=1.06)
legend(200,1.34,c('ATFcpue','residuals'),lty=c(1,2))
################################################################
# Plot GAM results:

# Use vis.gam for lat lon data:
pdf("vis.gam_ATFcpue_30_69.pdf", width=10, height=6)
vis.gam(ATF.cpue.gam_best, c("lon", "lat"), plot.type = "contour", type="response", contour.col="black", color="heat", xlab="Longitude", ylab="Latitude", main="Arrowtooth Flounder, Partial Effect on log-CPUE (number per hectare)", too.far=0.025, n.grid=500, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))
maps::map('worldHires', fill=T, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax), add=T, col="lightgrey")
dev.off()

# visreg provides better visualizations than plot(gam), but it splits the spatial smooth into separate fits for longitudinal effects (at mean latitude?) and latitudinal effects (at mean longitude?). Also, I'm not sure how the confidence intervals are constructed (F. Mueter)!
ATFcpueYRgam = visreg(ATF.cpue.gam_best, xvar="year", line=list(col="red"), fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), band=TRUE, partial=FALSE, rug=FALSE, trans=gaussian()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Arrowtooth Flounder") +
  theme(plot.title = element_text(size=18, hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="", y="Partial Effect on log-CPUE (number per hectare)") +
  theme(legend.background = element_rect(fill="transparent"))
ATFcpueYRgam

ggsave(filename="ATFcpueYRgam_30_69.png", plot=ATFcpueYRgam, dpi=500, width=12, height=8, units="in")

ATFcpueDepthGAM = visreg(ATF.cpue.gam_best, xvar="depth",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Arrowtooth Flounder") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x="Depth (m)", y="Partial Effect on log-CPUE (number per hectare)") +
  scale_x_continuous(expand=c(0,0), limits=c(0,1000), breaks=c(0,150,300,450,600,750,900)) +
  scale_y_continuous(limits=c(2,10), breaks=c(2,4,6,8,10)) +
  theme(legend.background = element_rect(fill="transparent")) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter") 
ATFcpueDepthGAM

ggsave(filename="ATFcpueDepthGAM_30_69.png", plot=ATFcpueDepthGAM, dpi=500, width=12, height=8, units="in")

ATFcpueBTgam = visreg(ATF.cpue.gam_best, "BT",  fill=list(col="gray", alpha=0.2), points=list(col="lightgray", alpha=0.5, cex=0.25), line=list(col="red"), band=TRUE, partial=FALSE, rug=1, trans=gaussian()$linkinv, type="conditional", gg=TRUE) +
  theme_bw() +
  ggtitle("Arrowtooth Flounder") +
  theme(plot.title = element_text(hjust = 0.5, family="Arial", size=18)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.0, size=14)) +
  labs(x=expression(paste("Bottom Temperature (",degree,"C)")), y="Partial Effect on log-CPUE (number per hectare)") +
  theme(legend.background = element_rect(fill="transparent")) +
  scale_x_continuous(expand=c(0,0), limits=c(0,13.5), breaks=c(0,3,6,9,12)) +
  scale_y_continuous(limits=c(7,10.7), breaks=c(7,8,9,10)) +
  geom_rug(sides="b", alpha=0.2, size=0.5, position="jitter")
ATFcpueBTgam

ggsave(filename="ATFcpueBTgam_30_69.png", plot=ATFcpueBTgam, dpi=500, width=12, height=8, units="in")

# Can visualize by multiple factors at once (but won't):
# visreg2d(ATF.pa.gam_best, "year", "BT", plot.type="image")

# Diagnostic Plots:
plot(ATF$year, resid(ATF.cpue.gam_best)); abline(h=0, lty=2)
plot(ATF$depth, resid(ATF.cpue.gam_best)); abline(h=0, lty=2)
# plot(trawl_wide$BT, resid(ATF.cpue.gamm_best_corr)); abline(h=0, lty=2)
plot(fitted(ATF.cpue.gam_best),resid(ATF.cpue.gam_best), # Evidence of heteroscedasticity?
     xlab="Fitted values", ylab="Residuals")
abline(h=0,lty=2); grid()

# Residuals by depth and over time:
r <- resid(ATF.cpue.gam_best)
j <- r > 0
symbols(ATF$year[j], ATF$depth[j], circles = r[j], bg = 4, inches=0.05) 
symbols(ATF$year[!j], ATF$depth[!j], circles = abs(r[!j]), bg = 2, inches=0.05, add=T) 

# Residuals by BT and over time:
symbols(ATF$year[j], ATF$BT[j], circles = r[j], bg = 4, inches=0.05) 
symbols(ATF$year[!j], ATF$BT[!j], circles = abs(r[!j]), bg = 2, inches=0.05, add=T) 

################################################################
##### Make Predictions - on GAMMs
# Read in grid cell centroid data:
grid_cent = read.csv("HaulCentEnviroData_for_Predictions_yr_100KM.csv")
grid_cent = na.omit(grid_cent)

grid_cent_yrGAM = grid_cent
grid_cent_yrGAM = grid_cent_yrGAM[, 2:8]
colnames(grid_cent_yrGAM) = c("year", "id2", "EEZgrid", "depth", "BT", "lon", "lat")
grid_cent_yrGAM$UniqueID = paste(grid_cent_yrGAM$id2, grid_cent_yrGAM$year, sep="_")

# PH presence or absence:
PHpa_predictGAM = predict.gam(PH.pa.gam_best, newdata=grid_cent_yrGAM, type="response", se.fit=TRUE)

PHpa_predGAM = cbind(grid_cent_yrGAM, PHpa_predictGAM)
PHpa_pred_CIgam = within(PHpa_predGAM, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(PHpa_pred_CIgam) = c("year", "id2", "EEZgrid", "depth", "BT", "lon", "lat", "UniqueID", "PHpa_fit", "PHpa_se.fit", "PHpa_lowerCI", "PHpa_upperCI")

PHpa_pred_CIgam %>%
  group_by(year) %>%
  summarize(meanP = mean(PHpa_fit)) 

# PH CPUE (where present):
PHcpue_predictGAM = predict.gam(PH.cpue.gam_best, newdata=grid_cent_yrGAM, type="response", se.fit=TRUE)

PHcpue_predGAM = cbind(grid_cent_yrGAM, PHcpue_predictGAM)
PHcpue_pred_CIgam = within(PHcpue_predGAM, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(PHcpue_pred_CIgam) = c("year", "id2", "EEZgrid", "depth", "BT", "lon", "lat", "UniqueID", "PHcpue_fit", "PHcpue_se.fit", "PHcpue_lowerCI", "PHcpue_upperCI")

PHpredictionsGAM = merge(PHpa_pred_CIgam, PHcpue_pred_CIgam, by="UniqueID", all=TRUE)
PHpredictionsGAM$PHpredBio = PHpredictionsGAM$PHpa_fit * PHpredictionsGAM$PHcpue_fit
PHpredictionsGAM$PHstdBio = PHpredictionsGAM$PHpredBio/max(PHpredictionsGAM$PHpredBio)

# ATF presence or absence:
ATFpa_predictGAM = predict.gam(ATF.pa.gam_best, newdata=grid_cent_yrGAM, type="response", se.fit=TRUE)

ATFpa_predGAM = cbind(grid_cent_yrGAM, ATFpa_predictGAM)
ATFpa_pred_CIgam = within(ATFpa_predGAM, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(ATFpa_pred_CIgam) = c("year", "id2", "EEZgrid", "depth", "BT", "lon", "lat", "UniqueID", "ATFpa_fit", "ATFpa_se.fit", "ATFpa_lowerCI", "ATFpa_upperCI")

# ATF CPUE (where present):
ATFcpue_predictGAM = predict.gam(ATF.cpue.gam_best, newdata=grid_cent_yrGAM, type="response", se.fit=TRUE)

ATFcpue_predGAM = cbind(grid_cent_yrGAM, ATFcpue_predictGAM)
ATFcpue_pred_CIgam = within(ATFcpue_predGAM, {
  lower = fit-1.96*se.fit
  upper = fit+1.96*se.fit
})

colnames(ATFcpue_pred_CIgam) = c("year", "id2", "EEZgrid", "depth", "BT", "lon", "lat", "UniqueID", "ATFcpue_fit", "ATFcpue_se.fit", "ATFcpue_lowerCI", "ATFcpue_upperCI")

ATFpredictionsGAM = merge(ATFpa_pred_CIgam, ATFcpue_pred_CIgam)
ATFpredictionsGAM$ATFpredBio = ATFpredictionsGAM$ATFpa_fit * ATFpredictionsGAM$ATFcpue_fit
ATFpredictionsGAM$ATFstdBio = ATFpredictionsGAM$ATFpredBio/max(ATFpredictionsGAM$ATFpredBio)

################################################################
# PH-ATF Spatial Overlap by grid cell:
PH_ATF_spatialGAM = merge(PHpredictionsGAM, ATFpredictionsGAM)

# Remove grid cells with std abundances < 0.25 for PH and ATF:
PH_ATF_spatialGAMred = PH_ATF_spatialGAM[!(PH_ATF_spatialGAM$PHstdBio < 0.25 & PH_ATF_spatialGAM$ATFstdBio < 0.25), ]

# Calculate Spatial Overlap:
PH_ATF_spatialGAMred$SpatOver = PH_ATF_spatialGAMred$PHstdBio * PH_ATF_spatialGAMred$ATFstdBio
PH_ATF_spatialGAMred$SpatOver = round(PH_ATF_spatialGAMred$SpatOver, digits=3)

PH_ATF_spatialGAMred = PH_ATF_spatialGAMred[ , c(1:12, 20:25, 33:43)]
colnames(PH_ATF_spatialGAMred) = c("UniqueID", "year", "id", "EEZ_grid", "depth", "BT", "Long", "Lat", "PHpa.fit", "PHpa.se", "PHpa.low", "PHpa.hi", "PHcpueFit", "PHcpueSE", "PHcpueLow", "PHcpueHi", "PHpredBio", "PHstdBio", "ATFpa.fit", "ATFpa.se", "ATFpa.low", "ATFpa.hi", "ATFcpueFit", "ATFcpueSE", "ATFcpueLow", "ATFcpueHi", "ATFpredBio", "ATFstdBio", "SpatOver")
PH_ATF_spatialGAMred$EEZ_grid = as.character(PH_ATF_spatialGAMred$EEZ_grid)
write.csv(PH_ATF_spatialGAMred, "PH_ATF_spatialGAMred_yr_max_100KM.csv")

################################################################
# Plot spatial overlap onto grid:
setwd("~/Documents/UAF/Dissertation/Analyses/MappingInfo/")

library(dplyr)
library(tidyr)
library(sp)
library(raster)
library(rgeos)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
require(purrr)
require(mapproj)
require(devtools)
require(stringr)
require(maptools)
require(rgdal)
library(PBSmapping)
require(mapdata)
require(ggplot2)
require(ggmap)

# Read NMFS Area 610-650 Shapefile:
NMFS = readOGR(".", "GOA_Shapes")

# Project layer:
NMFS_Pr = spTransform(NMFS, CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs"))

# plot(NMFS_Pr); axis(3); axis(4)
print(proj4string(NMFS_Pr))

# Dissolve individual polygons prior to making grid:
NMFSdata = data.frame()
NMFSdata = rbind(NMFSdata, NMFS_Pr@data)
NMFSdata$OBJECTID = as.character(NMFSdata$OBJECTID)
NMFSdata$REP_AREA = as.character(NMFSdata$REP_AREA)
NMFSdata$Region = NA
NMFSdata$Region = "GOA"

# Merge into polygons,
NMFS_Pr@data$OBJECTID = as.character(NMFS_Pr@data$OBJECTID)
NMFS_Pr@data = full_join(NMFS_Pr@data, NMFSdata, by = "REP_AREA")

# Ensure shapefile row.names and polygon IDs are sensible
row.names(NMFS_Pr) = row.names(NMFS_Pr@data)
NMFS_Pr = spChFIDs(NMFS_Pr, row.names(NMFS_Pr))

# Now dissolve
NMFS_Pr = gUnaryUnion(NMFS_Pr, id = NMFS_Pr@data$Region)

# Make sure row names match:
row.names(NMFS_Pr) = as.character(1:length(NMFS_Pr))

# Extract the data you want (the larger geography)
NMFSdata = unique(NMFSdata$Region)
NMFSdata = as.data.frame(NMFSdata)
colnames(NMFSdata) = "Region"  # your data will probably have more than 1 row!

# And add the data back in
NMFS_Pr = SpatialPolygonsDataFrame(NMFS_Pr, NMFSdata)

# Check it's all worked
# plot(NMFS_Pr)

#-------------------------------------------------------------
#  Make your grid - INPFC
#-------------------------------------------------------------
# Determine your grid size in degrees, assuming square grids:
my.interval=100 # May need to go smaller (e.g., 0.5)
# AFSC SAMPLING DESIGN:
# 5x5 km grids = individual stations
# tows ~ 1.5 km distance swept

# Specify your coordinates for the grid. Certain grid combinations don't work because they'll fall outside of the range of your EEZ layer or on the other side of -180, which gets wonky. With an interval set at 5 degree squares, pick a range that creates 5 degree increments. -177.5 will add 2.5 degrees in each direction (5 divided by 2), to yield a range starting at -180. If you use an interval of 1, you have to make sure that you don't start with like -180 and end with -124.5. It would need to be -180 and -124 or -125.

lonmin = -875
lonmax = 1975
latmin = 5075
latmax = 7975

mygrd = expand.grid(
  LON = seq(lonmin, lonmax, by=my.interval),
  LAT = seq(latmin, latmax, by=my.interval)) %>% 
  mutate(my.z=1:n()) %>% 
  data.frame
coordinates(mygrd) = ~ LON + LAT

# Convert SpatialPoints object to SpatialPixelsDataFrame:
# I got a funky error message sometimes on this line after I changed the dimensions of the grid cells above. I mostly just rebooted R and it seemed to work again. Magic!!!
mygrd = (as(SpatialPixelsDataFrame(mygrd, mygrd@data, tolerance=.00086), "SpatialPolygonsDataFrame"))

# See if it worked:
# plot(mygrd)

# Project mygrid to be same as NMFS file:
proj4string(mygrd) = proj4string(NMFS_Pr)

# Clip the two. If you use 1 degree squares, this will take a few minutes! Experiment with larger cells first. Once you get this to work on small cells, I'd recommend saving the output as an RData file so you don't have to redo. 
clip2_NMFS = gIntersection(NMFS_Pr, mygrd, byid = TRUE, drop_lower_td = TRUE)

# Project the clipped output:
proj4string(clip2_NMFS) = proj4string(NMFS_Pr)

# Make sure it worked:
# plot(clip2_NMFS, col="grey")
################################################################
# Make a copy while exploring so you don't have to rerun the clip2 intersection, which is slow.
clip3_NMFS = clip2_NMFS
clip3_NMFS = spTransform(clip3_NMFS, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Load data for map of N. Pacific (from PBSmapping):
data(nepacLLhigh)

lonmin = -172.5
lonmax = -129.5
latmin = 49
latmax = 62

# Clip the map to our spatial boundaries:
world = fortify(nepacLLhigh)
world2 = clipPolys(world, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

Canada = raster::getData("GADM", country = "CAN", level = 0)
Canada = fortify(Canada)
names(Canada) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Canada$PID = as.numeric(Canada$PID)
Canada = clipPolys(Canada, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

NMFS_shape = readOGR(".", "GOA_Shapes")
NMFS_shape = spTransform(NMFS_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
NMFS_shape = fortify(NMFS_shape)
names(NMFS_shape) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
NMFS_shape$PID = as.numeric(NMFS_shape$PID)
NMFS_shape = clipPolys(NMFS_shape, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

NMFS_shape1 = readOGR(".", "GOA_Shapes")
NMFS_shape1 = spTransform(NMFS_shape1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
NMFS_shape = fortify(NMFS_shape1)
names(NMFS_shape) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
NMFS_shape$PID = as.numeric(NMFS_shape$PID)
NMFS_shape = clipPolys(NMFS_shape, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))


NMFS_cent = coordinates(NMFS_shape1)
NMFS_cent = NMFS_cent[c(1, 3:6),]
NMFS_cent = as.data.frame(NMFS_cent)
names(NMFS_cent) = c("lon", "lat")
rownames(NMFS_cent) = c()
NMFS_cent$MgmtArea = NA
NMFS_cent[1,3] = "SE"
NMFS_cent[2,3] = "Shumagin"
NMFS_cent[3,3] = "Chirikof"
NMFS_cent[4,3] = "Kodiak"
NMFS_cent[5,3] = "Yakutat"
NMFS_cent$MgmtArea = as.factor(NMFS_cent$MgmtArea)

# Shift labels down to avoid overlap with grid cells:
NMFS_cent$lat[NMFS_cent$MgmtArea=="Kodiak"] = NMFS_cent$lat[NMFS_cent$MgmtArea=="Kodiak"] - 1.5
NMFS_cent$lat[NMFS_cent$MgmtArea=="SE"] = NMFS_cent$lat[NMFS_cent$MgmtArea=="SE"] - 0.8
NMFS_cent$lon[NMFS_cent$MgmtArea=="SE"] = NMFS_cent$lon[NMFS_cent$MgmtArea=="SE"] - 0.2
NMFS_cent$lat[NMFS_cent$MgmtArea=="Chirikof"] = NMFS_cent$lat[NMFS_cent$MgmtArea=="Chirikof"] - 1.15
NMFS_cent$lat[NMFS_cent$MgmtArea=="Yakutat"] = NMFS_cent$lat[NMFS_cent$MgmtArea=="Yakutat"] - 0.1
NMFS_cent$lat[NMFS_cent$MgmtArea=="Shumagin"] = NMFS_cent$lat[NMFS_cent$MgmtArea=="Shumagin"] - 0.4

# Read IPHC Area Shapefile:
IPHC_shape = readOGR(".", "GOA_Den")
IPHC_shape = spTransform(IPHC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
IPHC_2C = subset(IPHC_shape, REG_AREA=="2C")
IPHC_3A = subset(IPHC_shape, REG_AREA=="3A")
IPHC_3B = subset(IPHC_shape, REG_AREA=="3B")
IPHC_4A = subset(IPHC_shape, REG_AREA=="4A")

Area2C = gUnaryUnion(IPHC_2C, IPHC_2C@data$REG_Area)
Area3A = gUnaryUnion(IPHC_3A, IPHC_3A@data$REG_Area)
Area3B = gUnaryUnion(IPHC_3B, IPHC_3B@data$REG_Area)
Area4A = gUnaryUnion(IPHC_4A, IPHC_4A@data$REG_Area)

Area2C_df = fortify(Area2C)
row.names(Area2C_df) = c()
Area2C_df = subset(Area2C_df, group=="1.1")
# plot(Area2C_df$long, Area2C_df$lat)
Area3A_df = fortify(Area3A)
row.names(Area3A_df) = c()
Area3A_df = subset(Area3A_df, group=="1.1")
# plot(Area3A_df$long, Area3A_df$lat)
Area3B_df = fortify(Area3B)
row.names(Area3B_df) = c()
Area3B_df = subset(Area3B_df, group=="1.1")
# plot(Area3B_df$long, Area3B_df$lat)
Area4A_df = fortify(Area4A)

names(Area2C_df) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Area2C_df$PID = as.numeric(Area2C_df$PID)
Area2C_df = clipPolys(Area2C_df, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

names(Area3A_df) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Area3A_df$PID = as.numeric(Area3A_df$PID)
Area3A_df = clipPolys(Area3A_df, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

names(Area3B_df) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Area3B_df$PID = as.numeric(Area3B_df$PID)
Area3B_df = clipPolys(Area3B_df, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))

names(Area4A_df) = c("X", "Y", "POS", "hole", "piece", "id", "PID")
Area4A_df$PID = as.numeric(Area4A_df$PID)
Area4A_df = clipPolys(Area4A_df, xlim=c(lonmin, lonmax), ylim=c(latmin, latmax))


# Convert your dataset to a spatial data frame:
PH_ATF_spatialGAMred = read.csv("PH_ATF_spatialGAMred_yr_max_100KM.csv")
coordinates(PH_ATF_spatialGAMred) = c("Long", "Lat")
proj4string(PH_ATF_spatialGAMred) = proj4string(clip3_NMFS)

# Find where your data set overlaps with the clipped grid layer:
tempdat = data.frame(myrows=names(over(PH_ATF_spatialGAMred, clip3_NMFS)), mygrid=over(PH_ATF_spatialGAMred, clip3_NMFS))
PH_ATF_spatialGAMred$id2 = over(PH_ATF_spatialGAMred, clip3_NMFS)

# Make it a data frame again:
PH_ATF_spatialGAMred = as.data.frame(PH_ATF_spatialGAMred)
PH_ATF_spatialGAMred$EEZ_grid = as.character(PH_ATF_spatialGAMred$EEZ_grid)

PH_ATF_spatialGAMred = PH_ATF_spatialGAMred[,2:30]
PH_ATF_spatialGAMred = PH_ATF_spatialGAMred[!duplicated(PH_ATF_spatialGAMred), ]

### Summary INFORMATION
Overlap_Grid_spatial = PH_ATF_spatialGAMred[,c(1,2,28)]
colnames(Overlap_Grid_spatial) = c("year", "id2", "S")
# Create a data frame that includes the id values and the grid cell centroid lat-lon:
mycenter_NMFS = as.data.frame(gCentroid(clip3_NMFS, byid=TRUE)) %>% 
  mutate(id2=1:n(),
         EEZgrid=rownames(.))

# Join dietary overlap data (Overlap_Grid) with the centroid locations. 
str(mycenter_NMFS)
temp_NMFS = Overlap_Grid_spatial %>% 
  right_join(mycenter_NMFS)
temp_NMFS = as.data.frame(temp_NMFS)

# Convert the shapefile polygons to ggplot polygons
goa.df = fortify(clip3_NMFS, region='id')

# Join the data that contains the grid cell ids with the shapefile.
colnames(goa.df)[colnames(goa.df) == "id"] = "EEZgrid"
newdat_NMFS = goa.df %>% left_join(temp_NMFS)
newdat_NMFS = subset(newdat_NMFS, S >=0)

# Create INPFC polygon labels:
setwd("~/Documents/UAF/Dissertation/Analyses/MappingInfo/")
NMFS_shape1 = readOGR(".", "GOA_Shapes")
NMFS_shape1 = spTransform(NMFS_shape1, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Assign INPFC area to individual grid cells for summary stats:
summ_NMFS = newdat_NMFS
coordinates(summ_NMFS) = c("long", "lat")
proj4string(summ_NMFS) = proj4string(NMFS_shape1)

# tempID = data.frame(myrows=names(over(summ_NMFS, NMFS_shape1)), mygrid=over(summ_NMFS, NMFS_shape1))
summ_NMFS$MgmtArea = over(summ_NMFS, NMFS_shape1)
summ_NMFS = as.data.frame(summ_NMFS)
summ_NMFS = summ_NMFS[complete.cases(summ_NMFS), ]
# write.csv(summ_NMFS, "SpatialOverlapGrid_30_69_100KM.csv")

require(plyr)
summ_NMFS$GOAreg = revalue(summ_NMFS$MgmtArea.REP_AREA, c("610"="western", "620"="central", "630"="central", "640"="eastern", "650"="eastern"))

summ_NMFS2 = summ_NMFS[ , c(6, 8:10, 14:15,18)]
summ_NMFS2 = summ_NMFS2[!duplicated(summ_NMFS2), ]
summ_NMFS2 = subset(summ_NMFS2, GOAreg != "649")
summ_NMFS2 = subset(summ_NMFS2, GOAreg != "659")

require(dplyr)
# Calculate mean overlap by subregion:
meanS_region = summ_NMFS2 %>%
  group_by(GOAreg) %>%
  summarize(meanS = mean(S)) 
meanS_region 

sdS_region = summ_NMFS2 %>%
  group_by(GOAreg) %>%
  summarize(stdevS = sd(S)) 
sdS_region

range(summ_NMFS2$S)
mean(summ_NMFS2$S)
sd(summ_NMFS2$S)

# Calculate mean overlap by INPFC stat area:
meanS_area = summ_NMFS2 %>%
  group_by(MgmtArea.REP_AREA) %>%
  summarize(meanS = mean(S)) 
meanS_area

sdS_area = summ_NMFS2 %>%
  group_by(MgmtArea.REP_AREA) %>%
  summarize(stdevS = sd(S)) 
sdS_area

numS_area = summ_NMFS2 %>%
  group_by(MgmtArea.REP_AREA) %>%
  summarize(number = length(S)) 

meanS_yr = summ_NMFS2 %>%
  group_by(year) %>%
  summarize(meanS = mean(S)) 

sdS_yr = summ_NMFS2 %>%
  group_by(year) %>%
  summarize(stdevS = sd(S)) 

numS_yr = summ_NMFS2 %>%
  group_by(year) %>%
  summarize(number = length(S)) 

# Assign INPFC area to individual grid cells for summary stats:
summ_IPHC = newdat_NMFS
coordinates(summ_IPHC) = c("long", "lat")
proj4string(summ_IPHC) = proj4string(IPHC_shape)

# tempID = data.frame(myrows=names(over(summ_IPHC, IPHC_shape1)), mygrid=over(summ_IPHC, IPHC_shape1))
summ_IPHC$MgmtArea = over(summ_IPHC, IPHC_shape)
summ_IPHC = as.data.frame(summ_IPHC)
summ_IPHC = summ_IPHC[complete.cases(summ_IPHC), ]
# write.csv(summ_IPHC, "SpatialOverlapGrid_30_69_100KM.csv")

# Calculate mean dietary overlap by INPFC stat area:
summ_IPHC2 = summ_IPHC[ , c(6, 8:10, 18)]
summ_IPHC2 = summ_IPHC2[!duplicated(summ_IPHC2), ]

meanS_area = summ_IPHC2 %>%
  group_by(MgmtArea.REG_AREA) %>%
  summarize(meanS = mean(S)) 

sdS_area = summ_IPHC2 %>%
  group_by(MgmtArea.REG_AREA) %>%
  summarize(stdevS = sd(S)) 

numS_area = summ_IPHC2 %>%
  group_by(MgmtArea.REG_AREA) %>%
  summarize(number = length(S)) 

Overlap_Grid_meanYR = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanD = mean(SpatOver))

Overlap_Grid_stdYR = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(stdD = sd(SpatOver))

Overlap_Grid_mean_YR = PH_ATF_spatialGAMred %>%
  group_by(year, id) %>%
  mutate(meanD = mean(SpatOver))

Overlap_Grid_std_YR = PH_ATF_spatialGAMred %>%
  group_by(year, id) %>%
  mutate(stdD = sd(SpatOver))

PHpaGAMfit = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanPA = mean(PHpa.fit))

PHpaGAMse = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(sePA = mean(PHpa.se))

ATFpaGAMfit = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanPA = mean(ATFpa.fit))

ATFpaGAMse = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(sePA = mean(ATFpa.se))

PHcpueGAMfit = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanCPUE = mean(PHcpueFit))

PHcpueGAMse = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(seCPUE = mean(PHcpueSE))

ATFcpueGAMfit = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanCPUE = mean(ATFcpueFit))

ATFcpueGAMse = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(seCPUE = mean(ATFcpueSE))


StdBioPH = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanBio = mean(PHstdBio))

StdBioATF = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(meanBio = mean(ATFstdBio))

StdBioPHse = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(seBio = sd(PHstdBio))

StdBioATFse = PH_ATF_spatialGAMred %>%
  group_by(id) %>%
  mutate(seBio = sd(ATFstdBio))
#-------------------------------------------------------------
# Prepare for plotting, INPFC
#-------------------------------------------------------------
# Create a data frame that includes the id values and the grid cell centroid lat-lon:
mycenter_NMFS = as.data.frame(gCentroid(clip2_NMFS, byid=TRUE)) %>% 
  mutate(id2=1:n(),
         EEZgrid=rownames(.))
colnames(mycenter_NMFS)[3] = "id"
# Join dietary overlap data (Overlap_Grid) with the centroid locations. 
temp6_NMFS = Overlap_Grid_meanYR %>% 
  right_join(mycenter_NMFS)
temp6_NMFS = as.data.frame(temp6_NMFS)

temp7_NMFS = Overlap_Grid_stdYR %>% 
  right_join(mycenter_NMFS)
temp7_NMFS = as.data.frame(temp7_NMFS)

temp8_NMFS = Overlap_Grid_mean_YR %>% 
  right_join(mycenter_NMFS)
temp8_NMFS = as.data.frame(temp8_NMFS)

temp9_NMFS = Overlap_Grid_std_YR %>% 
  right_join(mycenter_NMFS)
temp9_NMFS = as.data.frame(temp9_NMFS)

temp10_NMFS = PHpaGAMfit %>% 
  right_join(mycenter_NMFS)
temp10_NMFS = as.data.frame(temp10_NMFS)

temp11_NMFS = PHpaGAMse %>% 
  right_join(mycenter_NMFS)
temp11_NMFS = as.data.frame(temp11_NMFS)

temp12_NMFS = ATFpaGAMfit %>% 
  right_join(mycenter_NMFS)
temp12_NMFS = as.data.frame(temp12_NMFS)

temp13_NMFS = ATFpaGAMse %>% 
  right_join(mycenter_NMFS)
temp13_NMFS = as.data.frame(temp13_NMFS)

temp14_NMFS = PHcpueGAMfit %>% 
  right_join(mycenter_NMFS)
temp14_NMFS = as.data.frame(temp14_NMFS)

temp15_NMFS = PHcpueGAMse %>% 
  right_join(mycenter_NMFS)
temp15_NMFS = as.data.frame(temp15_NMFS)

temp16_NMFS = ATFcpueGAMfit %>% 
  right_join(mycenter_NMFS)
temp16_NMFS = as.data.frame(temp16_NMFS)

temp17_NMFS = ATFcpueGAMse %>% 
  right_join(mycenter_NMFS)
temp17_NMFS = as.data.frame(temp17_NMFS)

temp18_NMFS = StdBioPH %>% 
  right_join(mycenter_NMFS)
temp18_NMFS = as.data.frame(temp18_NMFS)

temp19_NMFS = StdBioATF %>% 
  right_join(mycenter_NMFS)
temp19_NMFS = as.data.frame(temp19_NMFS)

temp21_NMFS = StdBioPHse %>% 
  right_join(mycenter_NMFS)
temp21_NMFS = as.data.frame(temp21_NMFS)

temp22_NMFS = StdBioATFse %>% 
  right_join(mycenter_NMFS)
temp22_NMFS = as.data.frame(temp22_NMFS)

# Convert the shapefile polygons to ggplot polygons
goa.df = fortify(clip3_NMFS, region='id')
#-------------------------------------------------------------
#  Now for plotting - Diet Data, INPFC
#-------------------------------------------------------------
# Join the data that contains the grid cell ids with the shapefile.
colnames(goa.df)[colnames(goa.df) == "id"] = "EEZgrid"
newdat6_NMFS = goa.df %>% left_join(temp6_NMFS)
newdat6_NMFS = subset(newdat6_NMFS, meanD >=0)

newdat7_NMFS = goa.df %>% left_join(temp7_NMFS)
newdat7_NMFS = subset(newdat7_NMFS, stdD >=0)

newdat8_NMFS = goa.df %>% left_join(temp8_NMFS)
newdat8_NMFS = subset(newdat8_NMFS, meanD >=0)

newdat9_NMFS = goa.df %>% left_join(temp9_NMFS)
newdat9_NMFS = subset(newdat9_NMFS, stdD >=0)

newdat10_NMFS = goa.df %>% left_join(temp10_NMFS)
newdat10_NMFS = subset(newdat10_NMFS, meanPA >=0)

newdat11_NMFS = goa.df %>% left_join(temp11_NMFS)
newdat11_NMFS = subset(newdat11_NMFS, sePA >=0)

newdat12_NMFS = goa.df %>% left_join(temp12_NMFS)
newdat12_NMFS = subset(newdat12_NMFS, meanPA >=0)

newdat13_NMFS = goa.df %>% left_join(temp13_NMFS)
newdat13_NMFS = subset(newdat13_NMFS, sePA >=0)

newdat14_NMFS = goa.df %>% left_join(temp14_NMFS)
newdat14_NMFS = subset(newdat14_NMFS, meanCPUE >=0)

newdat15_NMFS = goa.df %>% left_join(temp15_NMFS)
newdat15_NMFS = subset(newdat15_NMFS, seCPUE >=0)

newdat16_NMFS = goa.df %>% left_join(temp16_NMFS)
newdat16_NMFS = subset(newdat16_NMFS, meanCPUE >=0)

newdat17_NMFS = goa.df %>% left_join(temp17_NMFS)
newdat17_NMFS = subset(newdat17_NMFS, seCPUE >=0)

newdat18_NMFS = goa.df %>% left_join(temp18_NMFS)
newdat18_NMFS = subset(newdat18_NMFS, meanBio >=0)

newdat19_NMFS = goa.df %>% left_join(temp19_NMFS)
newdat19_NMFS = subset(newdat19_NMFS, meanBio >=0)

newdat21_NMFS = goa.df %>% left_join(temp21_NMFS)
newdat21_NMFS = subset(newdat21_NMFS, seBio >=0)

newdat22_NMFS = goa.df %>% left_join(temp22_NMFS)
newdat22_NMFS = subset(newdat22_NMFS, seBio >=0)

IPHC_shape = readOGR(".", "GOA_Den")
IPHC_shape = spTransform(IPHC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
IPHC_2C = subset(IPHC_shape, REG_AREA=="2C")
IPHC_3A = subset(IPHC_shape, REG_AREA=="3A")
IPHC_3B = subset(IPHC_shape, REG_AREA=="3B")
IPHC_4A = subset(IPHC_shape, REG_AREA=="4A")
#################################################################

# ANCOVAs
#################################################################
# Classify cells by INPFC stat areas:
ANCOVAdata = PH_ATF_spatialGAMred
coordinates(ANCOVAdata) = ~ Long + Lat

INPFC_shape = readOGR(".", "GOA_Shapes")
INPFC_shape = spTransform(INPFC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
proj4string(ANCOVAdata) = proj4string(INPFC_shape)
ANCOVAdata$INPFC = over(ANCOVAdata, INPFC_shape)

IPHC_shape = readOGR(".", "GOA_Den")
IPHC_shape = spTransform(IPHC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ANCOVAdata$IPHC = over(ANCOVAdata, IPHC_shape)

ANCOVAdf = as.data.frame(ANCOVAdata)

# Log-transform SpatOver to meet assumption of heterogeneity:
ANCOVAdf$logSpatOver = log(ANCOVAdata$SpatOver + 1)
ANCOVAdf$year = as.numeric(ANCOVAdf$year)

INPFCdf = subset(ANCOVAdf, INPFC.REP_AREA!="659")
INPFCdf = subset(INPFCdf, INPFC.REP_AREA!="NA")

IPHCdf = subset(ANCOVAdf, IPHC.REG_AREA!="2B")
IPHCdf = subset(IPHCdf, IPHC.REG_AREA!="NA")
IPHCdf$IPHC.REG_AREA = ordered(IPHCdf$IPHC.REG_AREA, levels=c("4A", "3B", "3A", "2C"))

# Run ANCOVA:
INPFCmodel1 = aov(SpatOver ~ year*INPFC.REP_AREA, data=INPFCdf)
summary(INPFCmodel1)
# plot(INPFCmodel1)

# Drop non-significant interaction term:
INPFCmodel2 = aov(SpatOver ~ year + INPFC.REP_AREA, data=INPFCdf)
summary(INPFCmodel2)
coefficients(INPFCmodel2)
# plot(INPFCmodel2)
termplot(INPFCmodel2, se=T, partial.resid=T)
# TukeyHSD(INPFCmodel2, "year")
TukeyHSD(INPFCmodel2, "INPFC.REP_AREA")

# Combine areas without significant differences in S:
INPFCdf$statAREA = INPFCdf$INPFC.REP_AREA
levels(INPFCdf$statAREA)
levels(INPFCdf$statAREA) = c("610", "620", "630", "eastGOA", "649", "eastGOA", "659")

INPFCdf %>%
  group_by(year) %>%
  summarize(meanS = mean(SpatOver)) 

INPFCdf %>%
  group_by(year) %>%
  summarize(sdS = sd(SpatOver)) 

INPFCdf %>%
  group_by(INPFCdf$statAREA) %>%
  summarize(meanS = mean(SpatOver)) 

INPFCdf %>%
  group_by(INPFCdf$statAREA) %>%
  summarize(sdS = sd(SpatOver)) 


# Run ANCOVA:
IPHCmodel1 = aov(SpatOver ~ year*IPHC.REG_AREA, data=IPHCdf)
summary(IPHCmodel1)
coefficients(IPHCmodel1)
# plot(IPHCmodel1)

# Drop non-significant interaction term:
IPHCmodel2 = aov(SpatOver ~ year + IPHC.REG_AREA, data=IPHCdf)
summary(IPHCmodel2)
# plot(IPHCmodel2)
termplot(IPHCmodel2, se=T, partial.resid=T)
# TukeyHSD(IPHCmodel2, "year")
TukeyHSD(IPHCmodel2, "IPHC.REG_AREA")

# Combine areas without significant differences in S:
IPHCdf$statAREA = IPHCdf$IPHC.REG_AREA
levels(IPHCdf$statAREA)
levels(IPHCdf$statAREA) = c("westGOA", "westGOA", "3A", "2C")

IPHCdf %>%
  group_by(IPHCdf$statAREA) %>%
  summarize(meanS = mean(SpatOver)) 

IPHCdf %>%
  group_by(IPHCdf$statAREA) %>%
  summarize(sdS = sd(SpatOver)) 
################################################################

################################################################
# Plot Empty Grid:
EmptyGrid = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat10_NMFS, aes(x=long, y=lat, group=group), fill=NA, col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.5) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

EmptyGrid
ggsave(filename="EmptyGrid_100KM.png", plot=EmptyGrid, dpi=500, width=12, height=8, units="in")

EmptyGrid_AFSC = ggplot() +
  geom_polygon(data=newdat10_NMFS, aes(x=long, y=lat, group=group), fill=NA, col="black", lwd = 0.25) + 
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.5) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray90", col="black", lwd=0.25) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

EmptyGrid_AFSC
ggsave(filename="EmptyGrid_100KM-AFSC.png", plot=EmptyGrid_AFSC, dpi=500, width=12, height=8, units="in")

EmptyGrid_AFSC_noGrid = ggplot() +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.5) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray90", col="black", lwd=0.25) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

EmptyGrid_AFSC_noGrid
ggsave(filename="EmptyGrid_100KM-AFSC_noGrid.png", plot=EmptyGrid_AFSC_noGrid, dpi=500, width=12, height=8, units="in")

# Plot mean annual dietary overlap by cell:
df = data.frame(x = -134.4, y = 61.7, text = c("Mean Interannual Spatial Overlap"))

# Plot mean PA by grid cell, PH:
textPA = data.frame(x = -135.3, y = 61.65, text = c("Mean P/A Predictions"))

# PACIFIC HALIBUT P/A:
PAplot_PH = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat10_NMFS, aes(x=long, y=lat, group=group, fill=meanPA), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textPA, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

PAplot_PH
ggsave(filename="PAplot_PH_30_69.png", plot=PAplot_PH, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat10_NMFS$year = as.factor(newdat10_NMFS$year)

for(var in unique(newdat10_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat10_NMFS[newdat10_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanPA), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textPA, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("PAplot_PH_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

textPAse = data.frame(x = -135.3, y = 61.65, text = c("StdDev P/A Predictions"))

PAplot_PHse = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat11_NMFS, aes(x=long, y=lat, group=group, fill=sePA), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,0.1), breaks=c(0,0.03, 0.06, 0.09)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textPAse, aes(label=text, x=-133.4, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent")) 

PAplot_PHse
ggsave(filename="PAplot_PHse_30_69.png", plot=PAplot_PHse, dpi=500, width=12, height=8, units="in")

for(var in unique(newdat11_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat11_NMFS[newdat11_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=sePA), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="") +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textPAse, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("PAplotSE_PH_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

### ARROWTOOTH FLOUNDER P/A:
PAplot_ATF = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat12_NMFS, aes(x=long, y=lat, group=group, fill=meanPA), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textPA, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

PAplot_ATF
ggsave(filename="PAplot_ATF_30_69.png", plot=PAplot_ATF, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat12_NMFS$year = as.factor(newdat12_NMFS$year)

for(var in unique(newdat12_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat12_NMFS[newdat12_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanPA), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textPA, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("PAplot_ATF_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

PAplot_ATFse = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat13_NMFS, aes(x=long, y=lat, group=group, fill=sePA), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,0.1), breaks=c(0,0.03, 0.06, 0.09)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textPAse, aes(label=text, x=-133.4, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent")) 

PAplot_ATFse
ggsave(filename="PAplot_ATFse_30_69.png", plot=PAplot_ATFse, dpi=500, width=12, height=8, units="in")

for(var in unique(newdat13_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat13_NMFS[newdat13_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=sePA), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="") +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textPAse, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("PAplotSE_ATF_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

# Plot mean CPUE by grid cell, PH:
textCPUE = data.frame(x = -139.3, y = 61.65, text = c("Mean CPUE Predictions"))

CPUEplot_PH = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat14_NMFS, aes(x=long, y=lat, group=group, fill=meanCPUE), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(3, 6.5), breaks=c(3,4,5,6)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textCPUE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

CPUEplot_PH
ggsave(filename="CPUEplot_PH_30_69.png", plot=CPUEplot_PH, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat14_NMFS$year = as.factor(newdat14_NMFS$year)

for(var in unique(newdat14_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat14_NMFS[newdat14_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanCPUE), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(3.6, 7.2), breaks=c(4.0, 5.0, 6.0, 7.0)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textCPUE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("CPUEplot_PH_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

textSE = data.frame(x = -135.3, y = 61.65, text = c("Std Dev CPUE Predictions"))

CPUEplot_PHse = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat15_NMFS, aes(x=long, y=lat, group=group, fill=seCPUE), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0, 2.1), breaks=c(0,0.5,1,1.5,2)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textSE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

CPUEplot_PHse
ggsave(filename="CPUEplot_PHse_30_69.png", plot=CPUEplot_PHse, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat15_NMFS$year = as.factor(newdat15_NMFS$year)

for(var in unique(newdat15_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat15_NMFS[newdat15_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=seCPUE), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="") +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textSE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("CPUEplot_PHse_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

### ARROWTOOTH FLOUNDER CPUE:
CPUEplot_ATF = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat16_NMFS, aes(x=long, y=lat, group=group, fill=meanCPUE), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(2, 10), breaks=c(2, 4, 6, 8, 10)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textCPUE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

CPUEplot_ATF
ggsave(filename="CPUEplot_ATF_30_69.png", plot=CPUEplot_ATF, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat16_NMFS$year = as.factor(newdat16_NMFS$year)

for(var in unique(newdat16_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat16_NMFS[newdat16_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanCPUE), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(2, 10), breaks=c(2, 4, 6, 8, 10)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textCPUE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("CPUEplot_ATF_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

CPUEplot_ATFse = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat17_NMFS, aes(x=long, y=lat, group=group, fill=seCPUE), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1.5), breaks=c(0,0.5,1,1.5)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textSE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

CPUEplot_ATFse
ggsave(filename="CPUEplot_ATFse_30_69.png", plot=CPUEplot_ATFse, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat17_NMFS$year = as.factor(newdat17_NMFS$year)

for(var in unique(newdat17_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat17_NMFS[newdat17_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=seCPUE), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="") +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textSE, aes(label=text, x=-135.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("CPUEplot_ATFse_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

### Std Biomass ###
text1 = data.frame(x = -139, y = 61.65, text = c("Std. Abundance (no. per ha)"))

### PACIFIC HALIBUT:
StdBioPHplot = ggplot() +
  geom_polygon(data=newdat18_NMFS, aes(x=long, y=lat, group=group, fill=meanBio), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00, 0.25, 0.50, 0.75, 1.0), guide="colourbar") +
  theme_bw() +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
  # geom_text(data=text1, aes(label=text, x=-134.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.957, 0.830)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=15), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(5.5, "mm")) +
  theme(legend.key.height = unit(8.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=16)) +
  theme(axis.text.x = element_text(family="Arial", size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

StdBioPHplot
ggsave(filename="StdBioPHplot_30_69.png", plot=StdBioPHplot, dpi=500, width=12, height=8, units="in")

StdBioPHplot_AFSC = ggplot() +
  geom_polygon(data=newdat18_NMFS, aes(x=long, y=lat, group=group, fill=meanBio), col="black", lwd = 0.25) + 
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00, 0.25, 0.50, 0.75, 1.0), guide="colourbar") +
  theme_bw() +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
  # geom_text(data=text1, aes(label=text, x=-134.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.957, 0.830)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=15), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(5.5, "mm")) +
  theme(legend.key.height = unit(8.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=16)) +
  theme(axis.text.x = element_text(family="Arial", size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))


StdBioPHplot_AFSC
ggsave(filename="StdBioPHplot_30_69_AFSC.png", plot=StdBioPHplot_AFSC, dpi=500, width=12, height=8, units="in")


# Probability of Occurrence by Year:
newdat18_NMFS$year = as.factor(newdat18_NMFS$year)

for(var in unique(newdat18_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat18_NMFS[newdat18_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanBio), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00, 0.25, 0.50, 0.75, 1.0)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=text1, aes(label=text, x=-134.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("StdBioPHplot_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

### ARROWTOOTH FLOUNDER:
StdBioATFplot = ggplot() +
  geom_polygon(data=newdat19_NMFS, aes(x=long, y=lat, group=group, fill=meanBio), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00, 0.25, 0.50, 0.75, 1.0)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray90", col="black", lwd=0.25) +
  # geom_text(data=text1, aes(label=text, x=-134.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.957, 0.830)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=15), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(5.5, "mm")) +
  theme(legend.key.height = unit(8.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=16)) +
  theme(axis.text.x = element_text(family="Arial", size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

StdBioATFplot
ggsave(filename="StdBioATFplot_30_69.png", plot=StdBioATFplot, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat19_NMFS$year = as.factor(newdat19_NMFS$year)

for(var in unique(newdat19_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat19_NMFS[newdat19_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanBio), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00, 0.25, 0.50, 0.75, 1.0)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=text1, aes(label=text, x=-134.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("StdBioATFplot_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

text2 = data.frame(x = -139, y = 61.65, text = c("StdDev Std Biomass"))

### PACIFIC HALIBUT:
StdBioPHplotse = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat21_NMFS, aes(x=long, y=lat, group=group, fill=seBio), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0.0,0.21), breaks=c(0.0,0.05, 0.10, 0.15, 0.2)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=text2, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

StdBioPHplotse
ggsave(filename="StdBioPHplotse_30_69.png", plot=StdBioPHplotse, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat21_NMFS$year = as.factor(newdat21_NMFS$year)

for(var in unique(newdat21_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat21_NMFS[newdat21_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=seBio), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0.0,2.21), breaks=c(0.00, 0.45, 0.90, 1.35, 1.80)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=text2, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("StdBioPHplot_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

### ARROWTOOTH FLOUNDER:
StdBioATFplotse = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat22_NMFS, aes(x=long, y=lat, group=group, fill=seBio), col="black", lwd = 0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0.0,0.21), breaks=c(0.0,0.05, 0.10, 0.15, 0.2)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=text2, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

StdBioATFplotse
ggsave(filename="StdBioATFplotse_30_69.png", plot=StdBioATFplotse, dpi=500, width=12, height=8, units="in")

# Probability of Occurrence by Year:
newdat22_NMFS$year = as.factor(newdat22_NMFS$year)

for(var in unique(newdat22_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat22_NMFS[newdat22_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=seBio), col="black", lwd=0.25) +   geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0.0,2.21), breaks=c(0.00, 0.45, 0.90, 1.35, 1.80)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=text2, aes(label=text, x=-133.15, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))   
  
  ggsave(p, filename=paste("StdBioATFplot_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

### Spatial Overlap ###
textZ = data.frame(x = -139, y = 61.65, text = c("Spatial Overlap"))

Spatial_Overlap_meanYrGAMM = ggplot() +
  geom_polygon(data=newdat6_NMFS, aes(x=long, y=lat, group=group, fill=meanD), col="black", lwd=0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray90", col="black", lwd=0.25) +
  # geom_text(data=textZ, aes(label=text, x=-132.35, y=61.65, size=12), show.legend=FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(panel.background = element_rect(fill="lightblue2", colour="black", size=1, line="solid")) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.957, 0.830)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=15), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(5.5, "mm")) +
  theme(legend.key.height = unit(8.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=16)) +
  theme(axis.text.x = element_text(family="Arial", size=16)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

Spatial_Overlap_meanYrGAMM
ggsave(filename="Spatial_Overlap_meanYrGAMM-Mod_30_69.png", plot=Spatial_Overlap_meanYrGAMM, dpi=500, width=12, height=8, units="in")

# Plot spatial overlap by grid cell and year (facet):
Spatial_Overlap_mean_YrGAMM = ggplot() +
  geom_polygon(data=newdat8_NMFS, aes(x=long, y=lat, group=group, fill=meanD), col="black") + 
  scale_fill_distiller(palette="Spectral", name=" S", limits=c(0.0, 01.0), breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), guide=FALSE) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.5) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.5) +
  theme_bw() +
  ggtitle("") +
  facet_wrap(~ year, ncol=3, dir="v") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.655, 0.225)) +
  theme(legend.title = element_text(size=14, vjust=2)) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial"), legend.text.align = 0) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0)) +
  theme(legend.background = element_rect(fill="transparent"))

Spatial_Overlap_mean_YrGAMM
ggsave(filename="Spatial_Overlap_mean_YrGAMM_30_69.png", plot=Spatial_Overlap_mean_YrGAMM, dpi=500, width=16, height=10.5, units="in")


# Plot spatial overlap by grid cell and year (for loop):
newdat6_NMFS$year = as.factor(newdat6_NMFS$year)

for(var in unique(newdat6_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat6_NMFS[newdat6_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=meanD), col="black", lwd=0.25) + 
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textZ, aes(label=text, x=-132.35, y=61.65, size=12), show.legend=FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(panel.background = element_rect(fill="lightblue2", colour="black", size=1, line="solid")) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))
  
  ggsave(p, filename=paste("Spatial_Overlap_Plot_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}

Spatial_Overlap_sdYrGAMM = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat7_NMFS, aes(x=long, y=lat, group=group, fill=stdD), col="black", lwd=0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,0.2), breaks=c(0.00, 0.05, 0.10, 0.15, 0.20)) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textZ, aes(label=text, x=-132.35, y=61.65, size=12), show.legend=FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(panel.background = element_rect(fill="lightblue2", colour="black", size=1, line="solid")) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.968, 0.863)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

Spatial_Overlap_sdYrGAMM
ggsave(filename="Spatial_Overlap_sdYrGAMM-Mod_30_69.png", plot=Spatial_Overlap_sdYrGAMM, dpi=500, width=12, height=8, units="in")

# Plot spatial overlap by grid cell and year (for loop):
newdat7_NMFS$year = as.factor(newdat7_NMFS$year)

for(var in unique(newdat7_NMFS$year)) {
  p = ggplot() +
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
    geom_polygon(data=newdat7_NMFS[newdat7_NMFS$year==var,], aes(x=long, y=lat, group=group, fill=stdD), col="black", lwd=0.25) + 
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
    geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_text(data=textZ, aes(label=text, x=-132.35, y=61.65, size=12), show.legend=FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(panel.background = element_rect(fill="lightblue2", colour="black", size=1, line="solid")) +
    geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
    theme_bw() +
    ggtitle("") +
    theme(legend.position = c(0.968, 0.863)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.title = element_blank()) +
    theme(legend.title.align = 0) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(1.9, "mm")) +
    theme(legend.key.height = unit(6.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=14)) +
    theme(axis.text.x = element_text(family="Arial", size=14)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
    theme(axis.title.y = element_text(vjust=1.1, size=16)) +
    labs(x="Longitude (dd)", y="Latitude (dd)") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))
  
  ggsave(p, filename=paste("Spatial_Overlap_Plot_sd_30_69_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}
#################################################################
# Re-group overlap values to show "significance":
### < 0.40, 0.40-0.59, > 0.60
newdat6_NMFS$OverBin = cut(newdat6_NMFS$meanD, breaks = c(-0.001, 0.199999, 0.399999, 0.599999, 1.0))
levels(newdat6_NMFS$OverBin) = c(" < 0.20", "   0.20 - 0.39", "   0.40 - 0.59", " > 0.60")
# Check bins:
table(newdat6_NMFS$meanD, newdat6_NMFS$OverBin)
#################################################################
Spatial_Overlap_Sig = ggplot() +
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="cadetblue1", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="deepskyblue2", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="blue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill="midnightblue", show.legend=FALSE, col=NA, linetype="solid", lwd=1.2) +  
  geom_polygon(data=newdat6_NMFS, aes(x=long, y=lat, group=group, fill=OverBin), col="black", lwd=0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_manual(values=c("green1", NA, "yellow1", "red1")) +
  geom_polygon(data=NMFS_shape, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_text(data=textZ, aes(label=text, x=-132.35, y=61.65, size=12), show.legend=FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(panel.background = element_rect(fill="lightblue2", colour="black", size=1, line="solid")) +
  geom_text(data=NMFS_cent, aes(group=MgmtArea, label=MgmtArea, x=lon, y=lat, size=12), show.legend = FALSE) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.94, 0.915)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(4.5, "mm")) +
  theme(legend.key.height = unit(4.5, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=16)) +
  theme(axis.title.y = element_text(vjust=1.1, size=16)) +
  labs(x="Longitude (dd)", y="Latitude (dd)") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

Spatial_Overlap_Sig 
ggsave(filename="Spatial_Overlap_Sig_30_69.png", plot=Spatial_Overlap_Sig , dpi=500, width=12, height=8, units="in")