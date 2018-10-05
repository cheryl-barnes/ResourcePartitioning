# This script file includes the code necessary to calculate species-specific proportions of prey by weight in each survey year and grid cell. We then used proportions of prey by weight (Chipps and Garvey 2007) to calculate Schoener's index of dietary overlap (Schoener 1968).

# We used standardized survey data procured from the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]). Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible here: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php. All the data necessary to complete the following analyses can be found in the 'Data' folder. See von Szalay and Raring (2016) and Livingston et al. (2017) for data collection methods.

# References:
# Chipps SR, Garvey JE. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. Guy CS, Brown ML, editors. Bethesda: American Fisheries Society; 2007. pp. 473–514.
# Livingston PA, Aydin K, Buckley TW, Lang GM, Yang MS, Miller BS. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes. 2017;100(4):443–470.
# Schoener TW. The anolis lizards of Bimini: resource partitioning in a complex fauna. Ecol. 1968;49(4):704–726.
# von Szalay PG, Raring NW. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle (WA): National Oceanic and Atmospheric Administration; 2016. Technical Memorandum: NMFS-AFSC-325. Sponsored by the US Department of Commerce.

rm(list=ls())
graphics.off()

setwd("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/ResourcePartitioning/ResourcePartitioning")

# DATA PREPARATION (START)
################################################################
# Format trawl data:
trawl = read.csv("barnes120717.csv")
trawl$STATIONID = as.factor(trawl$STATIONID)
trawl$SPECIES_CODE = as.factor(trawl$SPECIES_CODE)
levels(trawl$SPECIES_CODE)
levels(trawl$SPECIES_CODE) = c("ATF", "PH", "SBL", "PC", "WEP")

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

# Sum NUMBER_FISH by species, management area, and year:
require(plyr)

# Exclude incomplete entries not used in modeling (i.e., rows with missing depth, BT data): NO NEED TO DO THIS!
# trawl_comp = trawl[complete.cases(trawl), ]
# Exclude species not of interest:
trawl_comp = subset(trawl, SPECIES_CODE!="PC")
trawl_comp = subset(trawl_comp, SPECIES_CODE!="SBL")
trawl_comp = subset(trawl_comp, SPECIES_CODE!="WEP")

n_area_yr = trawl_comp %>% 
  group_by(SPECIES_CODE, MgmtArea, YEAR) %>% 
  summarise(SumFish = sum(NUMBER_FISH))
print(n_area_yr, n=400)

n_area = trawl_comp %>% 
  group_by(SPECIES_CODE, MgmtArea) %>% 
  summarise(SumFish = sum(NUMBER_FISH))
print(n_area, n=400)

n_yr = trawl_comp %>% 
  group_by(SPECIES_CODE, YEAR) %>% 
  summarise(SumFish = sum(NUMBER_FISH))
print(n_yr, n=400)
################################################################
##### Format food habits data:
### SEAK, 96-115 GW only:



# rm(list=ls())
setwd("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/Diet/")

# DATA PREPARATION (START)
##############################################################
##### Format food habits data:
preyWT = get(load("preyWTGOA_2013.Rdata"))
# write.csv(preyWT_wide, "preyWT_wide_1990-2013.csv")

# Rename preyWT_wide column to match preyWT_long (proceeding code previously written):
colnames(preyWT)[6] = "Predator"
preyWT = preyWT[,c(1,4,6,9:13,15:23,26:148)]

# Create own "HaulJoin" column by concatenating Vessel, Cruise, and Haul:
preyWT$Haul_Join = paste(preyWT$Vessel, preyWT$Cruise, preyWT$Haul, sep="")

# Assign Management Areas (values) based upon Survey Strata (index):
index = unique(preyWT$Strata)

values = c("610", "610", "610", "610", "630", "630", "630", "630", "630", "630", "640", "640", "640", "620", "620", "620", 
           "610", "610", "610", "610", "610", "620", "620", "620", "620", "630", "640", "640", "640", "640", "640", "640", 
           "650", "650", "650", "650", "650", "620", "630", "630", "630", "640", "630", "630", "630", "630", "620", "630", 
           "620", "630", "650", "650", "610", "640", "620") 

preyWT$MgmtArea = values[match(preyWT$Strata, index)]

# Make sure that all assignments were made correctly.
table(preyWT$Strata, preyWT$MgmtArea)

# Exclude 1981, 1984, and 1987 data:
preyWT = subset(preyWT, Yr!=1981)
preyWT = subset(preyWT, Yr!=1984)
preyWT = subset(preyWT, Yr!=1987)

# Calculate gape width based on fork length (empirical data collected in by Barnes and Gile, 2015):
PH = subset(preyWT, Predator=="PACIFIC HALIBUT")
ATF = subset(preyWT, Predator=="ARROWTOOTH FLOUNDR")

PH$GW = with(PH, PredL * 1.1966)
ATF$GW = with(ATF, PredL * 2.0639)

# Round gape width to nearest mm:
PH$GW = round(PH$GW, 0)
ATF$GW = round(ATF$GW, 0)

# Recombine for plotting:
PredPreyData = rbind(PH, ATF)

PredPreyData = subset(PredPreyData, GW >95)
PredPreyData = subset(PredPreyData, GW <116)
SEAK = subset(PredPreyData, MgmtArea == "650")

# Total number of stomachs sampled (w/ contents and w/o):
table(SEAK$Predator)
table(SEAK$Predator, SEAK$Yr)

# Remove all empty stomachs:
SEAK = SEAK[ , c(1:55, 57:143)]

# Reshape data, wide to long (selecting only columns of interest):
### Need Haul_Join as unique identifier - others rows will be duplicated in future...
require(reshape2)

colnames(SEAK)[31]
ncol(SEAK)
colnames(SEAK)[139]

preyWT_long = melt(SEAK, id.vars = c("ID", "Yr", "MgmtArea", "Station", "Haul_Join", "Predator", "PredL", "W_use", "GW"), measure.vars = 31:139, variable.name = "PreySpecies", value.name = "WT")

# Rename predators:
levels(as.factor(preyWT_long$Predator))
preyWT_long$Predator[preyWT_long$Predator == "ARROWTOOTH FLOUNDR"] = "ATF"
preyWT_long$Predator[preyWT_long$Predator == "PACIFIC HALIBUT"] = "PH"
levels(factor(preyWT_long$Predator))

# Remove all empty stomachs:
preyWT_long = subset(preyWT_long, PreySpecies!="Empty")

# Eliminate PreySpecies WT = 0 (removes empty stomachs from data frame):
preyWTcontents = subset(preyWT_long, WT > 0, select = c(ID, Predator, PredL, MgmtArea, Yr, Haul_Join, Station, PreySpecies, WT, GW))

# Convert back to wide format to calculate number of stomachs with contents:
preyWTcontents$WT = as.numeric(preyWTcontents$WT)
preyWTwide = dcast(preyWTcontents, Haul_Join + Predator + PredL + Yr + MgmtArea + Station ~ PreySpecies, value.var = "WT", fun.aggregate = sum)

# Calculate sample sizes by predator, management area, and year:
table(preyWTwide$MgmtArea, preyWTwide$Yr, preyWTwide$Predator)

# Calculate sample sizes by predator and management area:
table(preyWTwide$Predator, preyWTwide$MgmtArea)

# Calculate sample sizes by predator and year:
table(preyWTwide$Predator, preyWTwide$Yr)

# Calculate sample sizes by predator:
table(preyWTwide$Predator)

preyWTwide$Predator = as.factor(preyWTwide$Predator)
preyWTwide$Predator = factor(preyWTwide$Predator)

levels(preyWTcontents$PreySpecies)
# Rename prey species constituting less than 1% total diet:
levels(preyWTcontents$PreySpecies) = list(Other="Algae", Pleuronectiformes="AK.Plaice", Chondrichthyes="AK.Skate", Chondrichthyes="Aleutian.Skate", Cnidaria="Anemones", Pleuronectiformes="Arrow.or.Kam", Atheresthes.stomias="Arrowtooth", Pleurogrammus.monopterygius="Atka", Chionoecetes.sp="Bairdi", Argentiniformes="Bathylagidae", Crustacea="Benth.Amph", Other="Benthic..Hydroid", Other="Benth..Urochordata", Chondrichthyes="Big.Skate", Other="Birds", Echinodermata="Brittle.Star", Scorpaeniformes="Canary.Rock", Mallotus.villosus="Capelin", Other="Chaeteg.etc.", Mollusca="Clam", Crustacea="Copepod", Pleuronectiformes="Dover.Sole", Crustacea="Dungeness", Scorpaeniformes="Dusky.Rock", Perciformes="Eelpout", Osmeriformes="Eulachon", Crustacea="Euphausiid", Pleuronectiformes="FH.Sole", Teleostei="Fish.Eggs", Pleuronectiformes="Gen.Rock.Sole", Scorpaeniformes="Gen.Thorny", Mollusca="Gen..Cephalopod", Clupeiformes="Gen..Clupeids", Crustacea="Gen..Crab", Crustacea="Gen..Crustacea", Echinodermata="Gen..Echinoderm", Teleostei="Gen..Fish", Pleuronectiformes="Gen..Flatfish", Gadiformes="Gen..Gadid", Scorpaeniformes="Gen..Hexagrammidae", Cnidaria="Gen..Hydrozoa", Mollusca="Gen..Mollusc", Other="Gen..Particulate", Scorpaeniformes="Gen..Rockfish", Inorganic.Material="Gen..Rocks.et.al", Scorpaeniformes="Gen..Sebastes", Osmeriformes="Gen..Smelt", Chondrichthyes="Gen..Shark.Skate", Gadiformes="Giant.Grenadier", Other="Glopp", Pleuronectiformes="Gr..Turbot", Scorpaeniformes="Greenlings", Gadiformes="Hake", Paguroidea="Hermit.Crab", Clupea.pallasii="Herring", Pleuronectiformes="Kamchat.fl", Crustacea="King.Crab", Scorpaeniformes="Lg.Sculpin", Gadiformes="Macrouridae", Teleostei="Managed.Forage", Crustacea="Misc..Crab", Crustacea="Misc..Crustacean", Pleuronectiformes="Misc..Flatfish", Other="Misc..Worm..Etc.", Myctophiformes="Myctophidae", Crustacea="Mysid", Pleuronectiformes="N.Rock.Sole", Pandalidae="NP.Shrimp", Octopoda="Octopus", Other="Offal", Other="Opilo", Osmeriformes="Other.pel..Smelt", Gadus.macrocephalus="P..Cod", Pleuronectiformes="P..Halibut", Pandalidae="Pandalidae", Crustacea="Pel.Amph", Cnidaria="Pel..Gel..Filter.Feeder", Other="Polychaete", Scorpaeniformes="POP", Other="Prickle.squish.deep", Other="Prickle.squish.round", Other="Protozoan", Mollusca="Pteropod", Chondrichthyes="Ratfish", Pleuronectiformes="Rex.sole", Pleuronectiformes="S.Rock.Sole", Scorpaeniformes="Sablefish", Salmoniformes="Salmon", Ammodytes.hexapterus="Sandlance", Scorpaeniformes="Sculpin", Cnidaria="Scypho.Jellies", Cnidaria="Sea.Pens", Echinodermata="Sea.Star", Scorpaeniformes="Sebastes", Scorpaeniformes="Sharpchin.Rock", Scorpaeniformes="Shortsp.Thorny", Mollusca="Snail", Porifera="Sponge", Teuthida="Squid", Chondrichthyes="Unid.Bathyraja", Chionoecetes.sp="Unid.Chion", Chondrichthyes="Unid.Rajidae", Echinodermata="Urchins.dollars.cucumbers", Gadus.chalcogrammus="W..Pollock", Chondrichthyes="WhtBlotch.Skate", Pleuronectiformes="YF.Sole")

# Order predators:
preyWTcontents$Predator = ordered(preyWTcontents$Predator, levels = c("PH", "ATF"))

# Order prey items by phylogeny:
preyWTcontents$PreySpecies = ordered(preyWTcontents$PreySpecies, levels = c("Porifera", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Crustacea", "Pandalidae", "Chionoecetes.sp", "Paguroidea", "Echinodermata", "Chondrichthyes", "Teleostei", "Clupeiformes", "Clupea.pallasii", "Salmoniformes", "Osmeriformes", "Mallotus.villosus", "Argentiniformes", "Myctophiformes", "Gadiformes", "Gadus.chalcogrammus", "Gadus.macrocephalus", "Perciformes", "Ammodytes.hexapterus", "Scorpaeniformes", "Pleurogrammus.monopterygius", "Pleuronectiformes", "Atheresthes.stomias", "Offal", "Other", "Inorganic.Material"))

# DATA PREPARATION (END)
################################################################
################################################################
require(dplyr)
require(stats)
options(max.print = 1000)

##### Need to see about getting proportional scale rather than percent:
require(scales)
require(ggplot2)
require(gtable)
require(grid)

### Proportion of prey by weight; all areas, years, and size classes combined:
preyWTall = ggplot(preyWTcontents, aes(x = factor(Predator), y = WT, fill = PreySpecies, order = -as.numeric(PreySpecies))) +
  geom_bar(position = "fill", stat = "identity") +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(legend.title = element_text(size=14, face="bold", vjust=1.5)) +
  theme(legend.text = element_text(family="Arial", size=10.5), legend.text.align = 0) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_blank()) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  scale_fill_manual(values = c("Teuthida"="firebrick3", "Octopoda"="red", "Crustacea"="orange", "Pandalidae"="khaki1", "Teleostei"="cadetblue1", "Clupea.pallasii"="blue", "Salmoniformes"="blue4", "Gadus.chalcogrammus"="mediumorchid1", "Gadus.macrocephalus"="darkviolet", "Scorpaeniformes"="purple4", "Pleuronectiformes"="azure4", "Other"="black"), name="Prey Taxa", labels = c("Teuthida", "Octopodidae", "Crustacea", "Pandalidae", "Teleostei", expression(paste(italic("Clupea pallasii"))), "Salmoniformes", expression(paste(italic("  G. chalcogrammus"))), expression(paste(italic(" G. macrocephalus"))), "Misc. Rockfishes", "Misc. Flatfishes", "Other"), drop=TRUE) +
  labs(x="", y="") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0.001, 0.001)) +
  scale_x_discrete(limits = c("PH", "ATF"), expand = c(0.01, 0.01)) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  theme(legend.background = element_rect(fill="transparent"))

preyWTall

ggsave(filename="preyWT_SEAK_96_115.png", plot=preyWTall, width=6, height=8.5, units="in")
###########################################################
preyWT = get(load("preyWTGOA_2013.Rdata"))
# write.csv(preyWT_wide, "preyWT_wide_1990-2013.csv")

# Rename preyWT_wide column to match preyWT_long (proceeding code previously written):
colnames(preyWT)[6] = "Predator"

# Create own "HaulJoin" column by concatenating Vessel, Cruise, and Haul:
trawl$Haul_Join = paste(trawl$VESSEL, trawl$CRUISE, trawl$HAUL, sep="")

preyWT$Haul_Join = paste(preyWT$Vessel, preyWT$Cruise, preyWT$Haul, sep="")

# Assign Management Areas (values) based upon Survey Strata (index):
index = unique(preyWT$Strata)

values = c("610", "610", "610", "610", "630", "630", "630", "630", "630", "630", "640", "640", "640", "620", "620", "620", 
           "610", "610", "610", "610", "610", "620", "620", "620", "620", "630", "640", "640", "640", "640", "640", "640", 
           "650", "650", "650", "650", "650", "620", "630", "630", "630", "640", "630", "630", "630", "630", "620", "630", 
           "620", "630", "650", "650", "610", "640", "620") 

preyWT$MgmtArea = values[match(preyWT$Strata, index)]

# Make sure that all assignments were made correctly.
table(preyWT$Strata, preyWT$MgmtArea)

# Exclude 1981, 1984, and 1987 data:
preyWT = subset(preyWT, Yr!=1981)
preyWT = subset(preyWT, Yr!=1984)
preyWT = subset(preyWT, Yr!=1987)

# Total number of all stomachs sampled (w/ contents and empty):
table(preyWT$Predator)
table(preyWT$Predator, preyWT$MgmtArea)
table(preyWT$MgmtArea, preyWT$Yr, preyWT$Predator)

# Select only fish 30 - 89 cm:
preyWT_30_89 = subset(preyWT, PredL>29)
preyWT_30_89 = subset(preyWT_30_89, PredL<70)

# Total number of all stomachs sampled (w/ contents and empty):
table(preyWT_30_89$Predator)

# Reshape data, wide to long (selecting only columns of interest):
### Need Haul_Join as unique identifier - others rows will be duplicated in future...
require(reshape2)

colnames(preyWT)[39]
ncol(preyWT)
colnames(preyWT)[150]
colnames(preyWT)[149]
colnames(preyWT)[148]

preyWT_long = melt(preyWT, id.vars = c("PPID", "Yr", "MgmtArea", "Station", "Haul_Join", "RLONG", "RLAT", "Predator", "PredL", "W_use", "PredFull"), measure.vars = 39:148, variable.name = "PreySpecies", value.name = "WT")

# Rename predators:
levels(as.factor(preyWT_long$Predator))
preyWT_long$Predator[preyWT_long$Predator == "ARROWTOOTH FLOUNDR"] = "ATF"
preyWT_long$Predator[preyWT_long$Predator == "PACIFIC COD"] = "PC"
preyWT_long$Predator[preyWT_long$Predator == "PACIFIC HALIBUT"] = "PH"
preyWT_long$Predator[preyWT_long$Predator == "WALLEYE POLLOCK"] = "WEP"

# Remove all empty stomachs:
preyWT_long = subset(preyWT_long, PreySpecies!="Empty")

# Add depth data from trawl data frame:
colnames(trawl)[2] = "Yr"
colnames(trawl)[17] = "Station"
colnames(trawl)[21] = "BT"
trawlRed = trawl[ , c("Yr", "Haul_Join", "GEAR_DEPTH", "BT")]
trawlRed$GEAR_DEPTH = as.numeric(trawlRed$GEAR_DEPTH)

trawlRedagg1 = aggregate(GEAR_DEPTH ~ Yr + Haul_Join, trawlRed, mean)
trawlRedagg2 = aggregate(BT ~ Yr + Haul_Join, trawlRed, mean)
trawlRedagg = merge(trawlRedagg1, trawlRedagg2, by=c("Yr", "Haul_Join"), all=TRUE)

preyWTmerge = merge(preyWT_long, trawlRedagg, by=c("Yr", "Haul_Join"), all.x=TRUE)

### Set fork length bins (20 cm incremements: 10-29, 30-49, 50-69, 70-89):
preyWTmerge$FLBin = cut(preyWTmerge$PredL, breaks = c(0, 9, 29, 49, 69, 89, 171))
levels(preyWTmerge$FLBin) = c("<10", "10-29", "30-49", "50-69", "70-89", ">=90")

# Check bins:
options(max.print=100000000)
table(preyWTmerge$Predator, preyWTmerge$PredL, preyWTmerge$FLBin)

### Set depth bins (100 m incremements: < 100, 100-199, 200-299, 300-399, 400-499, >=500):
preyWTmerge$DepthBin = cut(preyWTmerge$GEAR_DEPTH, breaks = c(0, 99, 199, 299, 399, 499, 600))
levels(preyWTmerge$DepthBin) = c("<100", "100-199", "200-299", "300-399", "400-499", ">=500")
# Check bins:
table(preyWTmerge$GEAR_DEPTH, preyWTmerge$DepthBin)

# Rename management areas:
preyWTmerge$MgmtArea = as.factor(preyWTmerge$MgmtArea)
levels(preyWTmerge$MgmtArea) = c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern")

# Eliminate PreySpecies WT = 0 (removes empty stomachs from data frame):
preyWTcontents = subset(preyWTmerge, WT > 0, select = c(PPID, Predator, PredL,  FLBin, MgmtArea, Yr, Haul_Join, Station, RLONG, RLAT, PreySpecies, WT, PredFull, BT, GEAR_DEPTH, DepthBin))

# Convert back to wide format to calculate number of stomachs with contents:
preyWTcontents$WT = as.numeric(preyWTcontents$WT)
# write.csv(preyWTcontents, "preyWTcontents.csv") 
preyWTwide = dcast(preyWTcontents, Haul_Join + Predator + PredL + FLBin + Yr + MgmtArea + Station + GEAR_DEPTH + DepthBin ~ PreySpecies, value.var = "WT", fun.aggregate = sum)

# Calculate sample sizes by predator, management area, and year:
table(preyWTwide$MgmtArea, preyWTwide$Yr, preyWTwide$Predator)

# Calculate sample sizes by predator and management area:
table(preyWTwide$Predator, preyWTwide$MgmtArea)

# Calculate sample sizes by predator and year:
table(preyWTwide$Predator, preyWTwide$Yr)

# Calculate sample sizes by predator:
table(preyWTwide$Predator)

# Select only fish 30 - 89 cm:
preyWTwide_30_89 = subset(preyWTwide, FLBin==c("30-49", "50-69", "70-89"))

# Calculate sample sizes by predator, management area, and year:
table(preyWTwide_30_89$MgmtArea, preyWTwide_30_89$Yr, preyWTwide_30_89$Predator)

# Calculate sample sizes by predator and management area:
table(preyWTwide_30_89$Predator, preyWTwide_30_89$MgmtArea)

# Calculate sample sizes by predator and year:
table(preyWTwide_30_89$Predator, preyWTwide_30_89$Yr)

# Calculate sample sizes by predator:
table(preyWTwide_30_89$Predator)

# DATA PREPARATION (END)
################################################################

################################################################
# Note on significance of dietary overlap (Link and Auster 2013); General rule: > 40% merits consideration, > 60% significant (see Ross 1986, Garrison pers comm).

# Set up matrix to hold various overlap values:
MatrixCols = c("Group", "SampleSize", "NoPH", "NoATF", "OverlapMeasure")
OverlapMatrix = matrix(NA, nrow = 1000, ncol = length(MatrixCols))
OverlapMatrix = as.data.frame(OverlapMatrix)
colnames(OverlapMatrix) = c("Group", "SampleSize", "NoPH", "NoATF", "OverlapMeas")

# OverlapMatrix$SampleSize = the number of prey taxa used to calculate overlap.
################################################################
### DIETARY OVERLAP, BY PREDATOR ONLY:    

# Calculate total number of stomachs sampled by predator:
table(preyWT$Predator)

# Calculate number of stomachs w/ contents by predator:
table(preyWTwide$Predator)

# Individually calculate percentage of prey by weight for each predator:
PH = subset(preyWTcontents, Predator=="PH")
ATF = subset(preyWTcontents, Predator=="ATF")

require(dplyr)    
# Calculate PH and ATF proportional data:
PHprop = PH %>%
  group_by(PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHprop) = c("PreySpecies", "sum_PH", "prop_PH")
PHprop$prop_PH = round(PHprop$prop_PH, digits=5)

ATFprop = ATF %>%
  group_by(PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFprop) = c("PreySpecies", "sum_ATF", "prop_ATF")
ATFprop$prop_ATF = round(ATFprop$prop_ATF, digits=5)

# Rename prey species constituting less than 1% total diet:
levels(preyWTcontents$PreySpecies)
levels(preyWTcontents$PreySpecies) = list(Other="Algae", Pleuronectiformes="AK.Plaice", Chondrichthyes="AK.Skate", Chondrichthyes="Aleutian.Skate", Cnidaria="Anemones", Pleuronectiformes="Arrow.or.Kam", Atheresthes.stomias="Arrowtooth", Pleurogrammus.monopterygius="Atka", Chionoecetes.sp="Bairdi", Argentiniformes="Bathylagidae", Crustacea="Benth.Amph", Other="Benthic..Hydroid", Other="Benth..Urochordata", Chondrichthyes="Big.Skate", Other="Birds", Echinodermata="Brittle.Star", Scorpaeniformes ="Canary.Rock", Mallotus.villosus="Capelin", Other="Chaeteg.etc.", Mollusca="Clam", Crustacea="Copepod", Pleuronectiformes="Dover.Sole", Crustacea="Dungeness", Scorpaeniformes="Dusky.Rock", Perciformes="Eelpout", Empty="Empty", Osmeriformes="Eulachon", Euphausiacea="Euphausiid", Pleuronectiformes="FH.Sole", Teleostei="Fish.Eggs", Pleuronectiformes="Gen.Rock.Sole", Scorpaeniformes="Gen.Thorny", Mollusca="Gen..Cephalopod", Clupeiformes="Gen..Clupeids", Crustacea="Gen..Crab", Crustacea="Gen..Crustacea", Echinodermata="Gen..Echinoderm", Teleostei="Gen..Fish", Pleuronectiformes="Gen..Flatfish", Gadiformes="Gen..Gadid", Scorpaeniformes="Gen..Hexagrammidae", Cnidaria="Gen..Hydrozoa", Mollusca="Gen..Mollusc", Other="Gen..Particulate", Scorpaeniformes="Gen..Rockfish", Inorganic.Material="Gen..Rocks.et.al", Scorpaeniformes="Gen..Sebastes", Osmeriformes="Gen..Smelt", Chondrichthyes="Gen..Shark.Skate", Gadiformes="Giant.Grenadier", Other="Glopp", Pleuronectiformes="Gr..Turbot", Scorpaeniformes="Greenlings", Gadiformes="Hake", Paguroidea="Hermit.Crab", Clupea.pallasii="Herring", Pleuronectiformes="Kamchat.fl", Crustacea="King.Crab", Scorpaeniformes="Lg.Sculpin", Gadiformes="Macrouridae", Teleostei="Managed.Forage", Crustacea="Misc..Crab", Crustacea="Misc..Crustacean", Pleuronectiformes="Misc..Flatfish", Other="Misc..Worm..Etc.", Myctophiformes="Myctophidae", Crustacea="Mysid", Pleuronectiformes="N.Rock.Sole", Pandalidae="NP.Shrimp", Octopoda="Octopus", Offal="Offal", Other="Opilo", Osmeriformes="Other.pel..Smelt", Gadus.macrocephalus="P..Cod", Pleuronectiformes="P..Halibut", Pandalidae="Pandalidae", Crustacea="Pel.Amph", Cnidaria="Pel..Gel..Filter.Feeder", Other="Polychaete", Sebastes.alutus="POP", Other="Prickle.squish.deep", Other="Prickle.squish.round", Other="Protozoan", Mollusca="Pteropod", Chondrichthyes="Ratfish", Pleuronectiformes="Rex.sole", Pleuronectiformes="S.Rock.Sole", Scorpaeniformes="Sablefish", Salmoniformes="Salmon", Ammodytes.hexapterus="Sandlance", Scorpaeniformes="Sculpin", Cnidaria="Scypho.Jellies", Cnidaria="Sea.Pens", Echinodermata="Sea.Star", Scorpaeniformes="Sebastes", Scorpaeniformes="Sharpchin.Rock", Scorpaeniformes="Shortsp.Thorny", Mollusca="Snail", Porifera="Sponge", Teuthida="Squid", Chondrichthyes="Unid.Bathyraja", Chionoecetes.sp="Unid.Chion", Chondrichthyes="Unid.Rajidae", Echinodermata="Urchins.dollars.cucumbers", Gadus.chalcogrammus="W..Pollock", Chondrichthyes="WhtBlotch.Skate", Pleuronectiformes="YF.Sole")

# Individually calculate percentage of prey by weight for each predator (again):
PH = subset(preyWTcontents, Predator=="PH")
ATF = subset(preyWTcontents, Predator=="ATF")

# Calculate PH and ATF proportional data:
PHprop = PH %>%
  group_by(PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHprop) = c("PreySpecies", "sum_PH", "prop_PH")
PHprop$prop_PH = round(PHprop$prop_PH, digits=5)

ATFprop = ATF %>%
  group_by(PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFprop) = c("PreySpecies", "sum_ATF", "prop_ATF")
ATFprop$prop_ATF = round(ATFprop$prop_ATF, digits=5)

# Merge PH and ATF proportional data:
WTprop = merge(PHprop, ATFprop, by = "PreySpecies", all = TRUE)
# Set NA values for prey species to 0 (absent):  
WTprop[is.na(WTprop)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTprop$overlap = abs(WTprop$prop_PH - WTprop$prop_ATF)

# Calculate dietary overlap (Schoener 1970, second step):
PredOverlap = 1-0.5*(sum(WTprop$overlap))

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[1, 1] = "Predator"
OverlapMatrix[1, 2] = nrow(WTprop)
OverlapMatrix[1, 3] = nrow(subset(preyWTwide, Predator == "PH"))
OverlapMatrix[1, 4] = nrow(subset(preyWTwide, Predator == "ATF"))
OverlapMatrix[1, 5] = PredOverlap

################################################################
### DIETARY OVERLAP, BY FORK LENGTH:  

# Calculate number of stomachs w/ contents by predator and FL:
table(preyWTwide$Predator, preyWTwide$FLBin)
preyWTwide = subset(preyWTwide, FLBin!="<10")
preyWTwide = subset(preyWTwide, FLBin!=">=90")
preyWTwide = droplevels(preyWTwide)
table(preyWTwide$Predator, preyWTwide$FLBin)

# Calculate total number of stomachs sampled by predator and FL:
preyWTred = subset(preyWTcontents, FLBin!="<10")
preyWTred = subset(preyWTred, FLBin!=">=90")
preyWTred = droplevels(preyWTred)
table(preyWTred$Predator, preyWTred$FLBin)

PH = subset(preyWTred, Predator=="PH")
ATF = subset(preyWTred, Predator=="ATF")

# Individually calculate percentage of prey by weight for each predator:
PHpropFL = PH %>%
  group_by(FLBin, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropFL) = c("FL", "PreySpecies", "sum_PH", "prop_PH")

ATFpropFL = ATF %>%
  group_by(FLBin, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropFL) = c("FL", "PreySpecies", "sum_ATF", "prop_ATF")

# Merge PH and ATF proportional data:
WTpropFL = merge(PHpropFL, ATFpropFL, by = c("FL", "PreySpecies"), all = TRUE)

# Set NA values for prey species to 0 (absent):  
WTpropFL[is.na(WTpropFL)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTpropFL$overlap_pre = abs(WTpropFL$prop_PH - WTpropFL$prop_ATF)

# Calculate dietary overlap (Schoener 1970, second step):
aggPredFL = aggregate(overlap_pre ~ FL, WTpropFL, sum)
aggPredFL$overlap = 1-0.5*aggPredFL$overlap_pre

# Count number of rows for inclusion in OverlapMatrix:
nrow(aggPredFL)

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[2:5, 1] = paste(aggPredFL$FL, sep = ".")
WTpropFL_n = aggregate(overlap_pre ~ FL, WTpropFL, length)
OverlapMatrix[2:5, 2] = WTpropFL_n$overlap_pre

FLbyPred = aggregate(PredL ~ Predator + FLBin, preyWTwide, length, drop = FALSE)
FLbyPred = dcast(FLbyPred, FLBin ~ Predator, value.var = "PredL", fun.aggregate = sum)
OverlapMatrix[2:5, 3] = FLbyPred$PH
OverlapMatrix[2:5, 4] = FLbyPred$ATF
OverlapMatrix[2:5, 5] = aggPredFL$overlap

################################################################
### DIETARY OVERLAP, BY MGMT AREA:    
table(preyWTwide$Predator, preyWTwide$MgmtArea)

PH = subset(preyWTcontents, Predator=="PH")
ATF = subset(preyWTcontents, Predator=="ATF")

# Individually calculate percentage of prey by weight for each predator:
PHpropMgmt = PH %>%
  group_by(MgmtArea, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropMgmt) = c("MgmtArea", "PreySpecies", "sum_PH", "prop_PH")

ATFpropMgmt = ATF %>%
  group_by(MgmtArea, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropMgmt) = c("MgmtArea", "PreySpecies", "sum_ATF", "prop_ATF")

# Merge PH and ATF proportional data:
WTpropMgmt = merge(PHpropMgmt, ATFpropMgmt, by = c("MgmtArea", "PreySpecies"), all = TRUE)

# Set NA values for prey species to 0 (absent):  
WTpropMgmt[is.na(WTpropMgmt)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTpropMgmt$overlap_pre = abs(WTpropMgmt$prop_PH - WTpropMgmt$prop_ATF)

# Calculate dietary overlap (Schoener 1970, second step):
aggPredMgmt = aggregate(overlap_pre ~ MgmtArea, WTpropMgmt, sum)
aggPredMgmt$overlap = 1-0.5*aggPredMgmt$overlap_pre

# Count number of rows for inclusion in OverlapMatrix:
nrow(aggPredMgmt)

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[6:10, 1] = paste(aggPredMgmt$MgmtArea, sep = ".")
WTpropMgmt_n = aggregate(overlap_pre ~ MgmtArea, WTpropMgmt, length)
OverlapMatrix[6:10, 2] = WTpropMgmt_n$overlap_pre
MgmtPred = aggregate(PredL ~ Predator + MgmtArea, preyWTwide, length, drop = FALSE)
MgmtPred = dcast(MgmtPred, MgmtArea ~ Predator, value.var = "PredL", fun.aggregate = sum)
OverlapMatrix[6:10, 3] = MgmtPred$PH
OverlapMatrix[6:10, 4] = MgmtPred$ATF
OverlapMatrix[6:10, 5] = aggPredMgmt$overlap

################################################################
### DIETARY OVERLAP, BY YEAR:    
table(preyWTwide$Predator, preyWTwide$Yr)

PH = subset(preyWTcontents, Predator=="PH")
ATF = subset(preyWTcontents, Predator=="ATF")

# Individually calculate percentage of prey by weight for each predator:
PHpropYr = PH %>%
  group_by(Yr, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropYr) = c("Yr", "PreySpecies", "sum_PH", "prop_PH")

ATFpropYr = ATF %>%
  group_by(Yr, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropYr) = c("Yr", "PreySpecies", "sum_ATF", "prop_ATF")

# Merge PH and ATF proportional data:
WTpropYr = merge(PHpropYr, ATFpropYr, by = c("Yr", "PreySpecies"), all = TRUE)

# Set NA values for prey species to 0 (absent):  
WTpropYr[is.na(WTpropYr)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTpropYr$overlap_pre = abs(WTpropYr$prop_PH - WTpropYr$prop_ATF)

# Calculate dietary overlap (Schoener 1970, second step):
aggPredYr = aggregate(overlap_pre ~ Yr, WTpropYr, sum)
aggPredYr$overlap = 1-0.5*aggPredYr$overlap_pre

# Count number of rows for inclusion in OverlapMatrix:
nrow(aggPredYr)

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[11:21, 1] = paste(aggPredYr$Yr, sep = ".")
WTpropYr_n = aggregate(overlap_pre ~ Yr, WTpropYr, length)
OverlapMatrix[11:21, 2] = WTpropYr_n$overlap_pre
YrPred = aggregate(PredL ~ Predator + Yr, preyWTwide, length, drop = FALSE)
YrPred = dcast(YrPred, Yr ~ Predator, value.var = "PredL", fun.aggregate = sum)
OverlapMatrix[11:21, 3] = YrPred$PH
OverlapMatrix[11:21, 4] = YrPred$ATF
OverlapMatrix[11:21, 5] = aggPredYr$overlap

################################################################
### DIETARY OVERLAP, BY FORK LENGTH AND MGMT AREA:    
table(preyWTwide$Predator, preyWTwide$FLBin, preyWTwide$MgmtArea)

FL_MG = preyWTred
FL_MG$UniqueID = paste(FL_MG$FLBin, FL_MG$MgmtArea, sep=".")

# Remove 10-29 cm FL Southeastern (no PH sampled):
FL_MG = subset(FL_MG, UniqueID!="10-29.Southeastern")

PH = subset(FL_MG, Predator=="PH")
ATF = subset(FL_MG, Predator=="ATF")

# Individually calculate percentage of prey by weight for each predator:
PHpropFL_Mgmt = PH %>%
  group_by(FLBin, MgmtArea, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropFL_Mgmt) = c("FL", "MgmtArea", "PreySpecies", "sum_PH", "prop_PH")

ATFpropFL_Mgmt = ATF %>%
  group_by(FLBin, MgmtArea, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropFL_Mgmt) = c("FL", "MgmtArea", "PreySpecies", "sum_ATF", "prop_ATF")

# Merge PH and ATF proportional data:
WTpropFL_Mgmt = merge(PHpropFL_Mgmt, ATFpropFL_Mgmt, by = c("FL", "MgmtArea", "PreySpecies"), all = TRUE)

# Set NA values for prey species to 0 (absent):  
WTpropFL_Mgmt[is.na(WTpropFL_Mgmt)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTpropFL_Mgmt$overlap_pre = abs(WTpropFL_Mgmt$prop_PH - WTpropFL_Mgmt$prop_ATF)

# Calculate dietary overlap (Schoener 1970, second step):
aggPredFL_Mgmt = aggregate(overlap_pre ~ FL + MgmtArea, WTpropFL_Mgmt, sum)
aggPredFL_Mgmt$overlap = 1-0.5*aggPredFL_Mgmt$overlap_pre

# Count number of rows for inclusion in OverlapMatrix:
nrow(aggPredFL_Mgmt)

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[22:40, 1] = paste(aggPredFL_Mgmt$FL, aggPredFL_Mgmt$MgmtArea, sep = ".")
WTpropFL_Mgmt_n = aggregate(overlap_pre ~ FL + MgmtArea, WTpropFL_Mgmt, length)

OverlapMatrix[22:40, 2] = WTpropFL_Mgmt_n$overlap_pre
FLMgmtPred = aggregate(PredL ~ Predator + FLBin + MgmtArea, preyWTwide, length)
FLMgmtPred = dcast(FLMgmtPred, FLBin + MgmtArea ~ Predator, value.var = "PredL", fun.aggregate = sum)

FLMgmtPred = subset(FLMgmtPred, PH!=0)
FLMgmtPred = subset(FLMgmtPred, ATF!=0)

OverlapMatrix[22:40, 3] = FLMgmtPred$PH
OverlapMatrix[22:40, 4] = FLMgmtPred$ATF
OverlapMatrix[22:40, 5] = aggPredFL_Mgmt$overlap

################################################################
### DIETARY OVERLAP, BY FORK LENGTH AND YEAR: 
table(preyWTwide$Predator, preyWTwide$FLBin, preyWTwide$Yr)

PH = subset(preyWTred, Predator=="PH")
ATF = subset(preyWTred, Predator=="ATF")

# Individually calculate percentage of prey by weight for each predator:
PHpropFL_Yr = PH %>%
  group_by(FLBin, Yr, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropFL_Yr) = c("FL", "Yr", "PreySpecies", "sum_PH", "prop_PH")

ATFpropFL_Yr = ATF %>%
  group_by(FLBin, Yr, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropFL_Yr) = c("FL", "Yr", "PreySpecies", "sum_ATF", "prop_ATF")

# Merge PH and ATF proportional data:
WTpropFL_Yr = merge(PHpropFL_Yr, ATFpropFL_Yr, by = c("FL", "Yr", "PreySpecies"), all = TRUE)

# Set NA values for prey species to 0 (absent):  
WTpropFL_Yr[is.na(WTpropFL_Yr)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTpropFL_Yr$overlap_pre = abs(WTpropFL_Yr$prop_PH - WTpropFL_Yr$prop_ATF)

# Calculate dietary overlap (Schoener 1970, second step):
aggPredFL_Yr = aggregate(overlap_pre ~ FL + Yr, WTpropFL_Yr, sum)
aggPredFL_Yr$overlap = 1-0.5*aggPredFL_Yr$overlap_pre

# Count number of rows for inclusion in OverlapMatrix:
nrow(aggPredFL_Yr)

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[41:84, 1] = paste(aggPredFL_Yr$FL, aggPredFL_Yr$Yr, sep = ".")
WTpropFL_Yr_n = aggregate(overlap_pre ~ FL + Yr, WTpropFL_Yr, length)
OverlapMatrix[41:84, 2] = WTpropFL_Yr_n$overlap_pre
FLYrPred = aggregate(PredL ~ Predator + FLBin + Yr, preyWTwide, length, drop = FALSE)
FLYrPred = dcast(FLYrPred, FLBin + Yr ~ Predator, value.var = "PredL", fun.aggregate = sum)

# Remove rows with no PH or ATF sampled:
FLYrPred = FLYrPred[!(FLYrPred$PH == 0 | FLYrPred$ATF == 0), ] 
OverlapMatrix[41:84, 3] = FLYrPred$PH
OverlapMatrix[41:84, 4] = FLYrPred$ATF
OverlapMatrix[41:84, 5] = aggPredFL_Yr$overlap

################################################################
### DIETARY OVERLAP, BY FORK LENGTH, MGMT AREA, AND YEAR:    
table(preyWTwide$Predator, preyWTwide$FLBin, preyWTwide$MgmtArea, preyWTwide$Yr)

FMY = FL_MG
FMY$UniqueID = as.factor(paste(FMY$FLBin, FMY$MgmtArea, FMY$Yr, sep="."))

table1 = as.data.frame(table(FMY$Predator, FMY$UniqueID))
colnames(table1) = c("Predator", "UniqueID", "Freq")

FMY = merge(FMY, table1, by=c("UniqueID", "Predator"), all=TRUE)

FMY = subset(FMY, Freq!=0)

PH = subset(FMY, Predator=="PH")
ATF = subset(FMY, Predator=="ATF")

# Individually calculate percentage of prey by weight for each predator:
PHpropFL_Mgmt_Yr = PH %>%
  group_by(FLBin, MgmtArea, Yr, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropFL_Mgmt_Yr) = c("FL", "MgmtArea", "Yr", "PreySpecies", "sum_PH", "prop_PH")

ATFpropFL_Mgmt_Yr = ATF %>%
  group_by(FLBin, MgmtArea, Yr, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropFL_Mgmt_Yr) = c("FL", "MgmtArea", "Yr", "PreySpecies", "sum_ATF", "prop_ATF")

# Merge PH and ATF proportional data:
WTpropFL_Mgmt_Yr = merge(PHpropFL_Mgmt_Yr, ATFpropFL_Mgmt_Yr, by = c("FL", "MgmtArea", "Yr", "PreySpecies"), all = TRUE)

# Set NA values for prey species to 0 (absent):  
WTpropFL_Mgmt_Yr[is.na(WTpropFL_Mgmt_Yr)] = 0

# Calculate dietary overlap (Schoener 1970, first step):
WTpropFL_Mgmt_Yr$overlap_pre = abs(WTpropFL_Mgmt_Yr$prop_PH - WTpropFL_Mgmt_Yr$prop_ATF)

# Remove combinations with PH or ATF not sampled:
WTpropFL_Mgmt_Yr$UniqueID = paste(WTpropFL_Mgmt_Yr$FL, WTpropFL_Mgmt_Yr$MgmtArea, WTpropFL_Mgmt_Yr$Yr, sep=".")

# Calculate dietary overlap (Schoener 1970, second step):
aggPredFL_Mgmt_Yr = aggregate(overlap_pre ~ FL + MgmtArea + Yr, WTpropFL_Mgmt_Yr, sum)
aggPredFL_Mgmt_Yr$overlap = 1-0.5*aggPredFL_Mgmt_Yr$overlap_pre
aggPredFL_Mgmt_Yr$UniqueID = paste(aggPredFL_Mgmt_Yr$FL, aggPredFL_Mgmt_Yr$MgmtArea, aggPredFL_Mgmt_Yr$Yr, sep=".")

# Count number of rows for inclusion in OverlapMatrix:
nrow(aggPredFL_Mgmt_Yr)

FLMgmtYrPred = dcast(FMY, UniqueID ~ Predator, value.var = "PredL", fun.aggregate = sum)
colnames(OverlapMatrix)[1] = "UniqueID"
colnames(FLMgmtYrPred)[2] = "NoATF"
colnames(FLMgmtYrPred)[3] = "NoPH"
FLMgmtYrPred$OverlapMeas = NA
FLMgmtYrPred$SampleSize = NA
FLMgmtYrPred = FLMgmtYrPred[c(1, 5, 3, 2, 4)]

WTpropFL_Mgmt_Yr_n = aggregate(overlap_pre ~ FL + MgmtArea + Yr, WTpropFL_Mgmt_Yr, length)
WTpropFL_Mgmt_Yr_n$UniqueID = as.factor(paste(WTpropFL_Mgmt_Yr_n$FL, WTpropFL_Mgmt_Yr_n$MgmtArea, WTpropFL_Mgmt_Yr_n$Yr, sep="."))

# Add overlap group(s), value(s) and sample size(s) to matrix:
OverlapMatrix[85:266, 1] = paste(aggPredFL_Mgmt_Yr$FL, aggPredFL_Mgmt_Yr$MgmtArea, aggPredFL_Mgmt_Yr$Yr, sep = ".")
OverlapMatrix[85:266, 2] = WTpropFL_Mgmt_Yr_n$overlap_pre
OverlapMatrix[85:266, 5] = aggPredFL_Mgmt_Yr$overlap

# Match values from different data frames:
OverlapMatrix$NoPH = FLMgmtYrPred[match(OverlapMatrix$UniqueID, FLMgmtYrPred$UniqueID), 3]
OverlapMatrix$NoATF = FLMgmtYrPred[match(OverlapMatrix$UniqueID, FLMgmtYrPred$UniqueID), 4]
################################################################
# Add back in values that were deleted because I'm too tired to figure it out:
# Predator Only:
OverlapMatrix[1, 1] = "Predator"
OverlapMatrix[1, 2] = nrow(WTprop)
OverlapMatrix[1, 3] = nrow(subset(preyWTwide, Predator == "PH"))
OverlapMatrix[1, 4] = nrow(subset(preyWTwide, Predator == "ATF"))
OverlapMatrix[1, 5] = PredOverlap
#Fork Length Only:
OverlapMatrix[2:5, 3] = FLbyPred$PH
OverlapMatrix[2:5, 4] = FLbyPred$ATF
OverlapMatrix[2:5, 5] = aggPredFL$overlap
# Mgmt Area Only:
OverlapMatrix[6:10, 1] = paste(aggPredMgmt$MgmtArea, sep = ".")
OverlapMatrix[6:10, 2] = WTpropMgmt_n$overlap_pre
OverlapMatrix[6:10, 3] = MgmtPred$PH
OverlapMatrix[6:10, 4] = MgmtPred$ATF
OverlapMatrix[6:10, 5] = aggPredMgmt$overlap
# Year Only:
OverlapMatrix[11:21, 1] = paste(aggPredYr$Yr, sep = ".")
OverlapMatrix[11:21, 2] = WTpropYr_n$overlap_pre
OverlapMatrix[11:21, 3] = YrPred$PH
OverlapMatrix[11:21, 4] = YrPred$ATF
OverlapMatrix[11:21, 5] = aggPredYr$overlap
# Fork Length and Mgmt Area:
OverlapMatrix[22:40, 1] = paste(aggPredFL_Mgmt$FL, aggPredFL_Mgmt$MgmtArea, sep = ".")
OverlapMatrix[22:40, 2] = WTpropFL_Mgmt_n$overlap_pre
OverlapMatrix[22:40, 3] = FLMgmtPred$PH
OverlapMatrix[22:40, 4] = FLMgmtPred$ATF
OverlapMatrix[22:40, 5] = aggPredFL_Mgmt$overlap
# Fork Length and Yr:
OverlapMatrix[41:84, 1] = paste(aggPredFL_Yr$FL, aggPredFL_Yr$Yr, sep = ".")
OverlapMatrix[41:84, 2] = WTpropFL_Yr_n$overlap_pre
OverlapMatrix[41:84, 3] = FLYrPred$PH
OverlapMatrix[41:84, 4] = FLYrPred$ATF
OverlapMatrix[41:84, 5] = aggPredFL_Yr$overlap

OverlapMatrix = na.omit(OverlapMatrix)
write.csv(OverlapMatrix, file = "OverlapMatrix.csv")
################################################################
# Extract FL, MgmtArea, and Yr data to graph:
FL_Mgmt_Yr_Overlap = OverlapMatrix[85:266, ]

require(splitstackshape)
FL_Mgmt_Yr_Overlap = cSplit(FL_Mgmt_Yr_Overlap, "UniqueID", ".")
colnames(FL_Mgmt_Yr_Overlap) = c("SampleSize", "NoPH", "NoATF", "OverlapMeas", "FLBin", "MgmtArea", "Yr")
FL_Mgmt_Yr_Overlap$FLBin = ordered(FL_Mgmt_Yr_Overlap$FLBin, levels = c("<10", "10-29", "30-49", "50-69", "70-89", ">=90"))
levels(FL_Mgmt_Yr_Overlap$FLBin) = c("FL: < 10 cm", "FL: 10-29 cm", "FL: 30-49 cm", "FL: 50-69 cm", "FL: 70-89 cm", "FL: >= 90 cm")
FL_Mgmt_Yr_Overlap$MgmtArea = ordered(FL_Mgmt_Yr_Overlap$MgmtArea, levels = c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern"))

# Plot dietary overlap through time, but mgmt area and fork length:

# Exclude groups with less than 10 PH or ATF:
FL_Mgmt_Yr_Overlap = subset(FL_Mgmt_Yr_Overlap, NoPH >= 10 & NoATF >= 10)

require(ggplot2)
OverlapPlot = ggplot(FL_Mgmt_Yr_Overlap, aes(x = Yr, y = OverlapMeas, group = MgmtArea, color = MgmtArea)) +
  geom_point(aes(color = MgmtArea)) + 
  geom_line(aes(color = MgmtArea, linetype = MgmtArea)) +
  facet_wrap(~ FLBin, ncol=1) +
  theme_bw() +
  ggtitle("") +
  theme(strip.text.y = element_text(hjust = 0)) +
  theme(strip.text.x = element_text(hjust = 0)) +
  theme(plot.title = element_text(hjust = 0)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(legend.title = element_text(size=14, face="bold", vjust=1.5)) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=14, hjust=1)) +
  labs(x="", y="Dietary Overlap Index (Schoener 1970)") +
  guides(color=guide_legend(title="INPFC\n  Stat Area")) +
  scale_color_manual(values=c("firebrick3", "orange", "green3", "blue", "purple")) +
  scale_linetype(guide = "none") +
  scale_y_continuous(breaks = c(0, 0.4, 0.8)) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) +
  theme(legend.background = element_rect(fill="transparent"))

OverlapPlot
ggsave(filename="OverlapPlot_INPFC.png", plot=OverlapPlot, dpi=500, width=9, height=8, units="in")

OverlapPlot = ggplot(FL_Mgmt_Yr_Overlap, aes(x = Yr, y = OverlapMeas, group = FLBin, color = FLBin)) +
  geom_point(aes(color = FLBin)) + 
  geom_line(aes(color = FLBin, linetype = FLBin)) +
  facet_wrap(~ MgmtArea, ncol=1) +
  theme_bw() +
  ggtitle("") +
  theme(strip.text.y = element_text(hjust = 0)) +
  theme(strip.text.x = element_text(hjust = 0)) +
  theme(plot.title = element_text(hjust = 0)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(legend.title = element_text(size=14, face="bold", vjust=1.5)) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=14, hjust=1)) +
  labs(x="", y="Dietary Overlap Index (Schoener 1970)") +
  guides(color=guide_legend(title="Size Bin")) +
  scale_color_manual(values=c("firebrick3", "orange", "green3", "blue", "purple")) +
  scale_linetype(guide = "none") +
  scale_y_continuous(breaks = c(0, 0.4, 0.8)) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) +
  theme(legend.background = element_rect(fill="transparent"))

OverlapPlot
ggsave(filename="OverlapPlot_FL.png", plot=OverlapPlot, dpi=500, width=9, height=8, units="in")
################################################################
### LEFT OFF HERE:

# Plot proportion of prey by weight through time:
WTpropFL_Mgmt_Yr_long = melt(WTpropFL_Mgmt_Yr, id.vars = c("UniqueID", "Yr", "MgmtArea", "FL", "PreySpecies"), measure.vars = c(6, 8), variable.name = "Predator", value.name = "PropWT")
levels(WTpropFL_Mgmt_Yr_long$Predator)
levels(WTpropFL_Mgmt_Yr_long$Predator) = c("PH", "ATF")

propWT50_69 = subset(WTpropFL_Mgmt_Yr_long, FL=="50-69")
aggpropWT = aggregate(PropWT ~ PreySpecies + Predator, propWT50_69, mean)
aggpropWT$PropWT = round(aggpropWT$PropWT, digits=5)

aggpropWT = subset(aggpropWT, PreySpecies == c("Paguroidea", "Crustacea", "Gadus.chalcogrammus", "Clupea.pallasii", "Chionoecetes.sp", "Sebastes.alutus", "Pleurogrammus.monopterygius", "Ammodytes.hexapterus", "Offal", "Pandalidae", "Mallotus.villosus", "Euphausiacea", "Teleostei", "Osmeriformes", "Pleuronectiformes", "Gadus.macrocephalus"))

PropWTyr = ggplot(aggpropWT, aes(x = Yr, y = PropWT, group = PreySpecies, color = PreySpecies)) + 
  geom_line(aes(x = Yr, y = PropWT, group = PreySpecies)) +
  facet_wrap(~ MgmtArea, ncol=1) +
  theme_bw() +
  ggtitle("") +
  theme(strip.text.y = element_text(hjust = 0)) +
  theme(plot.title = element_text(hjust = 0)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(legend.title = element_text(size=14, face="bold", vjust=1.5)) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_color_manual(values = c("Porifera"="lightpink3", "Cnidaria"="lavenderblush2", "Mollusca"="brown4", "Teuthida"="firebrick3", "Octopoda"="red", "Crustacea"="chocolate", "Euphausiacea"="chocolate1", "Pandalidae"="orange", "Chionoecetes.sp"="gold", "Paguroidea"="khaki1", "Echinodermata"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupeiformes"="forestgreen", "Clupea.pallasii"="darkgreen", "Salmoniformes"="lightcyan1", "Osmeriformes"="cadetblue1", "Mallotus.villosus"="turquoise3", "Argentiniformes"="deepskyblue3", "Myctophiformes"="deepskyblue4", "Gadiformes"="blue", "Gadus.chalcogrammus"="mediumblue", "Gadus.macrocephalus"="blue4", "Perciformes"="mediumpurple1", "Ammodytes.hexapterus"="purple4", "Scorpaeniformes"="darkviolet", "Sebastes.alutus"="black", "Pleurogrammus.monopterygius"="hotpink1", "Pleuronectiformes"="mediumvioletred", "Atheresthes.stomias"="maroon4", "Offal"="gray91", "Other"="azure4", "Inorganic.Material"="black"), name="Prey Taxa", labels = c("Porifera", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Crustacea", "Euphausiacea", "Pandalidae", expression(paste(italic("Chionoecetes spp."))), "Paguroidea", "Echinodermata", "Chondrichthyes", "Teleostei", "Clupeiformes", expression(paste(italic("   Clupea pallasii"))), "Salmoniformes", "Osmeriformes", expression(paste(italic("   Mallotus villosus"))), "Argentiniformes", "Myctophiformes", "Gadiformes", expression(paste(italic("   Gadus chalcogrammus"))), expression(paste(italic("   Gadus macrocephalus"))), "Perciformes", expression(paste(italic("   Ammodytes hexapterus"))), "Scorpaeniformes", expression(paste(italic("   Sebastes alutus"))), expression(paste(italic("   Pleurogrammus monopterygius"))), "Pleuronectiformes", expression(paste(italic("   Atheresthes stomias"))), "Offal", "Other", "Inorganic Material"), drop=FALSE) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=14, hjust=1)) +
  labs(x="", y="Dietary Overlap Index (Schoener 1970)") +
  guides(color=guide_legend(title="Mgmt Area")) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8)) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011)) +
  theme(legend.background = element_rect(fill="transparent"))

PropWTyr
ggsave(filename="PropWTyr.png", plot=PropWTyr, dpi=500, width=9, height=8, units="in")

PropWTyr_smooth = ggplot(propWT50_69, aes(x = Yr, y = PropWT, group = PreySpecies, color = PreySpecies)) + 
  geom_smooth(aes(x = Yr, y = PropWT, group = PreySpecies, color = PreySpecies), method = lm, se = FALSE) +
  facet_wrap(~ MgmtArea, ncol=1) +
  theme_bw() +
  ggtitle("") +
  theme(strip.text.y = element_text(hjust = 0)) +
  theme(plot.title = element_text(hjust = 0)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(legend.title = element_text(size=14, face="bold", vjust=1.5)) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_color_manual(values = c("Porifera"="lightpink3", "Cnidaria"="lavenderblush2", "Mollusca"="brown4", "Teuthida"="firebrick3", "Octopoda"="red", "Crustacea"="chocolate", "Euphausiacea"="chocolate1", "Pandalidae"="orange", "Chionoecetes.sp"="gold", "Paguroidea"="khaki1", "Echinodermata"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupeiformes"="forestgreen", "Clupea.pallasii"="darkgreen", "Salmoniformes"="lightcyan1", "Osmeriformes"="cadetblue1", "Mallotus.villosus"="turquoise3", "Argentiniformes"="deepskyblue3", "Myctophiformes"="deepskyblue4", "Gadiformes"="blue", "Gadus.chalcogrammus"="mediumblue", "Gadus.macrocephalus"="blue4", "Perciformes"="mediumpurple1", "Scorpaeniformes"="darkviolet", "Ammodytes.hexapterus"="purple4", "Pleurogrammus.monopterygius"="hotpink1", "Pleuronectiformes"="mediumvioletred", "Atheresthes.stomias"="maroon4", "Offal"="gray91", "Other"="azure4", "Inorganic.Material"="black"), name="Prey Taxa", labels = c("Porifera", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Crustacea", "Euphausiacea", "Pandalidae", expression(paste(italic("Chionoecetes spp."))), "Paguroidea", "Echinodermata", "Chondrichthyes", "Teleostei", "Clupeiformes", expression(paste(italic("   Clupea pallasii"))), "Salmoniformes", "Osmeriformes", expression(paste(italic("   Mallotus villosus"))), "Argentiniformes", "Myctophiformes", "Gadiformes", expression(paste(italic("   Gadus chalcogrammus"))), expression(paste(italic("   Gadus macrocephalus"))), "Perciformes", expression(paste(italic("   Ammodytes hexapterus"))), "Scorpaeniformes", expression(paste(italic("   Pleurogrammus monopterygius"))), "Pleuronectiformes", expression(paste(italic("   Atheresthes stomias"))), "Offal", "Other", "Inorganic Material"), drop=FALSE) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=14, hjust=1)) +
  labs(x="", y="Dietary Overlap Index (Schoener 1970)") +
  guides(color=guide_legend(title="Mgmt Area")) +
  guides(fill=guide_legend(ncol=1, keyheight = 0.9)) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8)) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011)) +
  theme(legend.background = element_rect(fill="transparent"))

PropWTyr_smooth
ggsave(filename="PropWTyr_smooth.png", plot=PropWTyr_smooth, dpi=500, width=9, height=8, units="in")

### EXPLORING OVERLAP INDICES
# Eliminate rows with less than 10 PH or ATF:
OverlapMatrix = subset(OverlapMatrix, NoPH >= 10 & NoATF >= 10)
which.max(OverlapMatrix$OverlapMeas)
OverlapMatrix[128, ]
which.min(OverlapMatrix$OverlapMeas)
OverlapMatrix[74, ]
