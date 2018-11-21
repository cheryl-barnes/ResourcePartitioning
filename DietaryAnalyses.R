# This script file relies on results (e.g., spatial data frame manipulations, construction of uniform grid) from the 'SpatialAnalyses' file and includes code necessary to calculate species-specific proportions of prey by weight in each survey year and grid cell. We then used proportions of prey by weight (Chipps and Garvey 2007) to calculate Schoener's index of dietary overlap (Schoener 1968).

# We analyzed standardized survey data procured from the Alaska Fisheries Science Center (National Marine Fisheries Service, National Oceanic and Atmospheric Association [NOAA]). Food habits data (1990 to 2013) were provided by the Resource Ecology and Ecosystem Modeling (REEM) Program and are publicly accessible here: https://access.afsc.noaa.gov/REEM/WebDietData/DietDataIntro.php. All the data necessary to complete the following analyses can be found in the 'Data' folder. See von Szalay and Raring (2016) and Livingston et al. (2017) for data collection methods.

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) have been excluded.

# References:
# Chipps SR, Garvey JE. Assessment of diets and feeding patterns. In: Analysis and interpretation of freshwater fisheries data. Guy CS, Brown ML, editors. Bethesda: American Fisheries Society; 2007. pp. 473–514.
# Livingston PA, Aydin K, Buckley TW, Lang GM, Yang MS, Miller BS. Quantifying food web interactions in the North Pacific – a data-based approach. Environ Biol Fishes. 2017;100(4):443–470.
# Schoener TW. The anolis lizards of Bimini: resource partitioning in a complex fauna. Ecol. 1968;49(4):704–726.
# von Szalay PG, Raring NW. Data report: 2015 Gulf of Alaska bottom trawl survey. Seattle (WA): National Oceanic and Atmospheric Administration; 2016. Technical Memorandum: NMFS-AFSC-325. Sponsored by the US Department of Commerce.

setwd("~/Documents/UAF/Dissertation/GitHub/ResourcePartitioning/")
  # source(SpatialAnalyses.R)
################################################################
### INITIAL DATA PREPARATION ###
################################################################
# Read in and format food habits data:
preyWT = get(load("Data/preyWTGOA_2013.Rdata"))

# Create unique haul identifier by concatenating VESSEL, CRUISE, and HAUL (allows for joining to trawl survey data below):
preyWT$Haul_Join = paste(preyWT$Vessel, preyWT$Cruise, preyWT$Haul, sep="")

# Manually assign International North Pacific Fisheries Commissionstatistical areas (i.e., values) based on survey strata (i.e., index). Note: The second number from the right corresponds with individual statistical areas (e.g., STRATUM 120 = StatArea 620, STRATUM 251 = StatArea 650).
index = unique(preyWT$Strata)
values = c("610", "610", "610", "610", "630", "630", "630", "630", "630", "630", "640", "640", "640", "620", "620", "620", 
           "610", "610", "610", "610", "610", "620", "620", "620", "620", "630", "640", "640", "640", "640", "640", "640", 
           "650", "650", "650", "650", "650", "620", "630", "630", "630", "640", "630", "630", "630", "630", "620", "630", 
           "620", "630", "650", "650", "610", "640", "620") 
preyWT$StatArea = values[match(preyWT$Strata, index)]
table(preyWT$Strata, preyWT$StatArea) # Check

# Rename INPFC statistical areas:
preyWT$StatArea = as.factor(preyWT$StatArea)
levels(preyWT$StatArea) = c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "Southeastern")

# Exclude data from 1981, 1984 and 1987 (survey methods were standardized in 1990):
preyWT = subset(preyWT, Yr!=1981)
preyWT = subset(preyWT, Yr!=1984)
preyWT = subset(preyWT, Yr!=1987)

# Select only predators of interest:
preyWT = subset(preyWT, Species_name==c("PACIFIC HALIBUT", "ARROWTOOTH FLOUNDR"))

# Limit data to size classes of interest (30 to 69 cm):
preyWT_red = subset(preyWT, PredL > 29)
preyWT_red = subset(preyWT_red, PredL < 70)
table(preyWT_red$Species_name)

# Set 10 cm fork length bins (30-39, 40-49, 50-59, 60-69):
preyWT_red$FLBin = cut(preyWT_red$PredL, breaks = c(0, 29, 39, 49, 59, 69, 171))
levels(preyWT_red$FLBin) = c("<30", "30-39", "40-49", "50-59", "60-69", ">=70")
table(preyWT_red$Species_name, preyWT_red$PredL, preyWT_red$FLBin) # Check

# Remove all empty stomachs:
preyWT_red$sumWT = rowSums(preyWT_red[,39:148])
preyWT_red = preyWT_red[ , c(1:63, 65:152)]
preyWTcontents = subset(preyWT_red, preyWT_red$sumWT > 0)

# Compute sample sizes for proportion of prey by weight plot:
table(preyWTcontents$StatArea, preyWTcontents$FLBin, preyWTcontents$Species_name)
################################################################
# Join spatial and diet data:
require(dplyr)
require(stats)
require(scales)
require(ggplot2)
require(gtable)
require(grid)
options(max.print = 1000)

# Find where data overlap with the clipped grid and convert back to data frame:
coordinates(preyWTcontents) = c("RLONG", "RLAT")
proj4string(preyWTcontents) = proj4string(clipTrawl)

# Find where data overlap with the clipped grid and convert back to data frame:
tempdat = data.frame(myrows=names(over(preyWTcontents, clipTrawl)), mygrid=over(preyWTcontents, clipTrawl))
preyWTcontents$id2 = over(preyWTcontents, clipTrawl)
preyWTcontents = as.data.frame(preyWTcontents)

# Create data frame with id values and grid cell coordinates:
mycenter_INPFC = as.data.frame(gCentroid(clip_INPFC, byid=TRUE)) %>% 
  mutate(id2=1:n(),
         EEZgrid=rownames(.))
colnames(mycenter_INPFC)[3] = "id"
xy=mycenter_INPFC[,c(1,2)]
mycenter_INPFC_sp = SpatialPointsDataFrame(coords=xy, data=mycenter_INPFC, proj4string=CRS("+proj=utm +zone=5, +datum=WGS84 +units=km +no_defs"))
mycenter_INPFC_Pr = spTransform(mycenter_INPFC_sp, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "))
mycenter_INPFC_df = as.data.frame(mycenter_INPFC_Pr) 
names(mycenter_INPFC_df) = c("x.UTM", "y.UTM", "id2", "EEZgrid", "x.lon", "y.lat")

# Join trawl and grid cell data (and clean up):
mycenter_INPFC_D = mycenter_INPFC_df %>% right_join(preyWTcontents)

# Relabel predator names:
levels(mycenter_INPFC_D$Species_name) = c("ATF", "PH")

# Reshape data, wide to long (selecting only columns of interest):
require(reshape2)
preyWT_long = melt(mycenter_INPFC_D, id.vars = c("id2", "EEZgrid", "Yr", "StatArea", "Haul_Join", "RLONG", "RLAT", "Species_name", "PredL", "FLBin"), measure.vars = 45:153, variable.name = "PreySpecies", value.name = "WT")
################################################################
### CALCULATE PREY DIVERSITY AND EVENNESS ###
################################################################
# Calculate prey diversity and evenness, by predator:
require(vegan)
Hsum = preyWT_long %>%
  group_by(Species_name, PreySpecies) %>%
  summarize(sumWT = sum(WT))
Hcalc = spread(Hsum, key = PreySpecies, value = sumWT)
Hcalc[is.na(Hcalc)] = 0

# PACIFIC HALIBUT #
H_PH = subset(Hcalc, Species_name=="PH")
H_PH = H_PH[,2:110]
Hprime_PH = diversity(H_PH)
Hprime_PH

H_PH$S = rowSums(H_PH > 0)
J_PH = Hprime_PH/log(ncol(H_PH))
J_PH

# ARROWTOOTH FLOUNDER #
H_ATF = subset(Hcalc, Species_name=="ATF")
H_ATF = H_ATF[,2:110]
Hprime_ATF = diversity(H_ATF)
Hprime_ATF

H_ATF$S = rowSums(H_ATF > 0)
J_ATF = Hprime_ATF/log(ncol(H_ATF))
J_ATF
################################################################
### CALCULATE AND PLOT PROPORTIONS OF PREY BY WEIGHT ###
################################################################
# Subset diet by predator for calculating proportions of prey by weight:
PH_all = subset(preyWT_long, Species_name=="PH")
ATF_all = subset(preyWT_long, Species_name=="ATF")

PHprop_all = PH_all %>%
  group_by(PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHprop_all) = c("PreySpecies", "sum_PH", "prop_PH")
PHprop_all$prop_PH = round(PHprop_all$prop_PH, digits=5)

ATFprop_all = ATF_all %>%
  group_by(PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFprop_all) = c("PreySpecies", "sum_ATF", "prop_ATF")
ATFprop_all$prop_ATF = round(ATFprop_all$prop_ATF, digits=5)

# Regroup prey taxa constituting less than 1% total diet by weight (for plotting only):
preyWTgrouped = preyWT_long
levels(preyWTgrouped$PreySpecies) = list(Other="Algae", Pleuronectiformes="AK.Plaice", Chondrichthyes="AK.Skate", Chondrichthyes="Aleutian.Skate", Cnidaria="Anemones", Pleuronectiformes="Arrow.or.Kam", Atheresthes.stomias="Arrowtooth", Pleurogrammus.monopterygius="Atka", Chionoecetes.sp="Bairdi", Argentiniformes="Bathylagidae", Crustacea="Benth.Amph", Other="Benthic..Hydroid", Other="Benth..Urochordata", Chondrichthyes="Big.Skate", Other="Birds", Echinodermata="Brittle.Star", Scorpaeniformes="Canary.Rock", Mallotus.villosus="Capelin", Other="Chaeteg.etc.", Mollusca="Clam", Crustacea="Copepod", Pleuronectiformes="Dover.Sole", Crustacea="Dungeness", Scorpaeniformes="Dusky.Rock", Perciformes="Eelpout", Osmeriformes="Eulachon", Euphausiacea="Euphausiid", Pleuronectiformes="FH.Sole", Teleostei="Fish.Eggs", Pleuronectiformes="Gen.Rock.Sole", Scorpaeniformes="Gen.Thorny", Mollusca="Gen..Cephalopod", Clupeiformes="Gen..Clupeids", Crustacea="Gen..Crab", Crustacea="Gen..Crustacea", Echinodermata="Gen..Echinoderm", Teleostei="Gen..Fish", Pleuronectiformes="Gen..Flatfish", Gadiformes="Gen..Gadid", Scorpaeniformes="Gen..Hexagrammidae", Cnidaria="Gen..Hydrozoa", Mollusca="Gen..Mollusc", Other="Gen..Particulate", Scorpaeniformes="Gen..Rockfish", Inorganic.Material="Gen..Rocks.et.al", Scorpaeniformes="Gen..Sebastes", Osmeriformes="Gen..Smelt", Chondrichthyes="Gen..Shark.Skate", Gadiformes="Giant.Grenadier", Other="Glopp", Pleuronectiformes="Gr..Turbot", Scorpaeniformes="Greenlings", Gadiformes="Hake", Paguroidea="Hermit.Crab", Clupea.pallasii="Herring", Pleuronectiformes="Kamchat.fl", Crustacea="King.Crab", Scorpaeniformes="Lg.Sculpin", Gadiformes="Macrouridae", Teleostei="Managed.Forage", Crustacea="Misc..Crab", Crustacea="Misc..Crustacean", Pleuronectiformes="Misc..Flatfish", Other="Misc..Worm..Etc.", Myctophiformes="Myctophidae", Crustacea="Mysid", Pleuronectiformes="N.Rock.Sole", Pandalidae="NP.Shrimp", Octopoda="Octopus", Offal="Offal", Other="Opilo", Osmeriformes="Other.pel..Smelt", Gadus.macrocephalus="P..Cod", Pleuronectiformes="P..Halibut", Pandalidae="Pandalidae", Crustacea="Pel.Amph", Cnidaria="Pel..Gel..Filter.Feeder", Other="Polychaete", Sebastes.alutus="POP", Other="Prickle.squish.deep", Other="Prickle.squish.round", Other="Protozoan", Mollusca="Pteropod", Chondrichthyes="Ratfish", Pleuronectiformes="Rex.sole", Pleuronectiformes="S.Rock.Sole", Scorpaeniformes="Sablefish", Salmoniformes="Salmon", Ammodytes.hexapterus="Sandlance", Scorpaeniformes="Sculpin", Cnidaria="Scypho.Jellies", Cnidaria="Sea.Pens", Echinodermata="Sea.Star", Scorpaeniformes="Sebastes", Scorpaeniformes="Sharpchin.Rock", Scorpaeniformes="Shortsp.Thorny", Mollusca="Snail", Porifera="Sponge", Teuthida="Squid", Chondrichthyes="Unid.Bathyraja", Chionoecetes.sp="Unid.Chion", Chondrichthyes="Unid.Rajidae", Echinodermata="Urchins.dollars.cucumbers", Gadus.chalcogrammus="W..Pollock", Chondrichthyes="WhtBlotch.Skate", Pleuronectiformes="YF.Sole")

# Order prey items by phylogeny:
preyWTgrouped$PreySpecies = ordered(preyWTgrouped$PreySpecies, levels = c("Porifera", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Crustacea", "Euphausiacea", "Pandalidae", "Chionoecetes.sp", "Paguroidea", "Echinodermata", "Chondrichthyes", "Teleostei", "Clupeiformes", "Clupea.pallasii", "Salmoniformes", "Osmeriformes", "Mallotus.villosus", "Argentiniformes", "Myctophiformes", "Gadiformes", "Gadus.chalcogrammus", "Gadus.macrocephalus", "Perciformes", "Ammodytes.hexapterus", "Scorpaeniformes", "Sebastes.alutus", "Pleurogrammus.monopterygius", "Pleuronectiformes", "Atheresthes.stomias", "Offal", "Other", "Inorganic.Material"))
################################################################
### CALCULATE AND PLOT PROPORTIONS OF PREY BY WEIGHT ###
################################################################# 
# Calculate proportions for each survey year-grid cell combination (using new prey taxa groupings - for plot only):
PH_grouped = subset(preyWTgrouped, Species_name=="PH")
ATF_grouped = subset(preyWTgrouped, Species_name=="ATF")

detach(package:plyr)
PHprop_grouped = PH_grouped %>%
  group_by(FLBin, StatArea, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHprop_grouped) = c("FLBin", "StatArea", "PreySpecies", "n", "prop")
PHprop_grouped$prop = round(PHprop_grouped$prop, digits=5)
PHprop_grouped$Species_name = "PH"

ATFprop_grouped = ATF_grouped %>%
  group_by(FLBin, StatArea, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFprop_grouped) = c("FLBin", "StatArea", "PreySpecies", "n", "prop")
ATFprop_grouped$prop = round(ATFprop_grouped$prop, digits=5)
ATFprop_grouped$Species_name = "ATF"

# Merge PH and ATF proportional data:
WTprop_grouped = rbind(PHprop_grouped, ATFprop_grouped)

# Order and relabel predators:
WTprop_grouped$Species_name = ordered(WTprop_grouped$Species_name, levels = c("PH", "ATF"))
levels(WTprop_grouped$Species_name) = c("Pacific Halibut", "Arrowtooth Flounder")

# Plot proportions of prey by predator, size class, and stat area:
quartz()
preyWTplot = ggplot(WTprop_grouped, aes(x = FLBin, y = prop, fill = PreySpecies, order = -as.numeric(PreySpecies))) + 
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(Species_name ~ StatArea) + 
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(colour="black", size=1, line="solid")) +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  theme(panel.spacing.y = unit(0.8, "lines")) +
  theme(legend.title = element_text(size=11, face="bold", vjust=1.5)) + 
  theme(legend.text = element_text(family="Arial", size=9.5), legend.text.align = 0) +
  theme(axis.text = element_text(family="Arial", size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(colour="white",fill="white")) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  scale_fill_manual(values = c("Porifera"="lightpink3", "Cnidaria"="lavenderblush2", "Mollusca"="brown4", "Teuthida"="firebrick3", "Octopoda"="red", "Crustacea"="chocolate", "Euphausiacea"="chocolate1", "Pandalidae"="orange", "Chionoecetes.sp"="gold", "Paguroidea"="khaki1", "Echinodermata"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupeiformes"="forestgreen", "Clupea.pallasii"="darkgreen", "Salmoniformes"="lightcyan1", "Osmeriformes"="cadetblue1", "Mallotus.villosus"="turquoise3", "Argentiniformes"="deepskyblue3", "Myctophiformes"="deepskyblue4", "Gadiformes"="blue", "Gadus.chalcogrammus"="mediumblue", "Gadus.macrocephalus"="blue4", "Perciformes"="mediumpurple1", "Ammodytes.hexapterus"="darkviolet", "Scorpaeniformes"="purple4", "Sebastes.alutus"="black", "Pleurogrammus.monopterygius"="hotpink1", "Pleuronectiformes"="mediumvioletred", "Atheresthes.stomias"="maroon4", "Offal"="gray91", "Other"="azure4", "Inorganic.Material"="black"), name="", labels = c("Porifera", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Crustacea", "Euphausiacea", "Pandalidae", expression(paste(italic("Chionoecetes spp."))), "Paguroidea", "Echinodermata", "Chondrichthyes", "Teleostei", "Clupeiformes", expression(paste(italic("   Clupea pallasii"))), "Salmoniformes", "Osmeriformes", expression(paste(italic("   Mallotus villosus"))), "Argentiniformes", "Myctophiformes", "Gadiformes", expression(paste(italic("   Gadus chalcogrammus"))), expression(paste(italic("   Gadus macrocephalus"))), "Perciformes", expression(paste(italic("   Ammodytes hexapterus"))), "Scorpaeniformes", expression(paste(italic("   Sebastes alutus"))), expression(paste(italic("   Pleurogrammus monopterygius"))), "Pleuronectiformes", expression(paste(italic("   Atheresthes stomias"))), "Offal", "Other", "Inorganic Material"), drop=FALSE) +
scale_fill_manual(values = c("Porifera"="lightpink3", "Cnidaria"="lavenderblush2", "Mollusca"="brown4", "Teuthida"="firebrick3", "Octopoda"="red", "Crustacea"="chocolate", "Euphausiacea"="chocolate1", "Pandalidae"="orange", "Chionoecetes.sp"="gold", "Paguroidea"="khaki1", "Echinodermata"="lemonchiffon", "Chondrichthyes"="darkolivegreen1", "Teleostei"="chartreuse2", "Clupeiformes"="forestgreen", "Clupea.pallasii"="darkgreen", "Salmoniformes"="lightcyan1", "Osmeriformes"="cadetblue1", "Mallotus.villosus"="turquoise3", "Argentiniformes"="deepskyblue3", "Myctophiformes"="deepskyblue4", "Gadiformes"="blue", "Gadus.chalcogrammus"="mediumblue", "Gadus.macrocephalus"="blue4", "Perciformes"="mediumpurple1", "Ammodytes.hexapterus"="darkviolet", "Scorpaeniformes"="purple4", "Sebastes.alutus"="black", "Pleurogrammus.monopterygius"="hotpink1", "Pleuronectiformes"="mediumvioletred", "Atheresthes.stomias"="maroon4", "Offal"="gray91", "Other"="azure4", "Inorganic.Material"="black"), name="Prey Taxa", labels = c("Porifera", "Cnidaria", "Mollusca", "Teuthida", "Octopoda", "Crustacea", "Euphausiacea", "Pandalidae", expression(paste(italic("Chionoecetes spp."))), "Paguroidea", "Echinodermata", "Chondrichthyes", "Teleostei", "Clupeiformes", expression(paste(italic("   Clupea pallasii"))), "Salmoniformes", "Osmeriformes", expression(paste(italic("   Mallotus villosus"))), "Argentiniformes", "Myctophiformes", "Gadiformes", expression(paste(italic("   Gadus chalcogrammus"))), expression(paste(italic("   Gadus macrocephalus"))), "Perciformes", expression(paste(italic("   Ammodytes hexapterus"))), "Scorpaeniformes", expression(paste(italic("   Sebastes alutus"))), expression(paste(italic("   Pleurogrammus monopterygius"))), "Pleuronectiformes", expression(paste(italic("   Atheresthes stomias"))), "Offal", "Other", "Inorganic Material"), drop=FALSE) +
  labs(x="Fork Length (cm)", y="Proportion by Weight") +
  scale_y_continuous(breaks = c(0, 0.5, 1.0), expand = c(0, 0.01)) +
  scale_x_discrete(expand = c(0.03, 0.03)) +
  guides(fill=guide_legend(nrow=5)) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box.margin = margin(0,0,0,0), legend.background = element_rect(fill="transparent")) +
  theme(legend.direction = "horizontal", legend.position = "bottom", legend.box.margin = margin(-20,-10,0,-20), legend.background = element_rect(fill="transparent"))
preyWTplot

ggsave(filename="Plots/preyWTplot.png", plot=preyWTplot, dpi=500, width=11.5, height=7, units="in")
################################################################
### CALCULATE AND PLOT MEAN DIETARY OVERLAP ###
#################################################################
# Eliminate survey year-grid cell combinations with less than three fish sampled for either species:
Nfish = mycenter_INPFC_D %>% 
  group_by(Yr, id2, Species_name) %>% 
  summarise(NoFish=length(unique(PPID)))
Nfish_wide = spread(Nfish, key = Species_name, value = NoFish)
Nfish_wide[is.na(Nfish_wide)] = 0
Nfish_sub = subset(Nfish_wide, ATF >= 3)
Nfish_sub = subset(Nfish_sub, PH >= 3)
preyWT_longGrid = Nfish_sub %>% left_join(preyWT_long)
################################################################
# Summarize diet data (ungrouped prey taxa for overlap calculations):

# PACIFIC HALIBUT #
PH = subset(preyWT_longGrid, Species_name =="PH")
PHpropGrid = PH %>%
  group_by(Yr, id2, EEZgrid, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n)) 
colnames(PHpropGrid) = c("Yr", "id2", "EEZgrid", "PreySpecies", "n_PH", "prop_PH")
PHpropGrid$prop_PH = round(PHpropGrid$prop_PH, digits=5)

# ARROWTOOTH FLOUNDER #
ATF = subset(preyWT_longGrid, Species_name =="ATF")
ATFpropGrid = ATF %>%
  group_by(Yr, id2, EEZgrid, PreySpecies) %>% 
  summarize(n = sum(WT)) %>%
  mutate(proportion = n / sum(n))
colnames(ATFpropGrid) = c("Yr", "id2", "EEZgrid", "PreySpecies", "n_ATF", "prop_ATF")
ATFpropGrid$prop_ATF = round(ATFpropGrid$prop_ATF, digits=5)

# Merge PH and ATF proportional data:
WTprop = merge(PHpropGrid, ATFpropGrid, by = c("id2", "EEZgrid", "Yr", "PreySpecies"), all = TRUE)
WTprop[is.na(WTprop)] = 0 # Set NA values to 0 (absent) 

# Remove prey taxa not consumed by either predator:
WTprop$sumWT = rowSums(WTprop[,c("prop_PH", "prop_ATF")])
WTprop = subset(WTprop, sumWT > 0)
################################################################
# Calculate dietary overlap (Schoener 1968):
WTprop$overlap = abs(WTprop$prop_PH - WTprop$prop_ATF)

D_Grid = WTprop %>%
  group_by(Yr, id2, EEZgrid) %>%
  summarize(n = sum(overlap)) %>%
  mutate(D = 1-0.5*(n))
D_Grid$D = round(D_Grid$D, digits=3)
D_Grid = D_Grid[ , c(1:3, 5)]

# Write CSV for analyses of resource partitioning:
write.csv(D_Grid, "Data/PH_ATF_D.csv")
  # D_overlap = read.csv("Data/PH_ATF_D.csv")

# Summarize by grid cell (all years combined):
D_Grid_mean = D_Grid %>%
  group_by(id2) %>%
  mutate(meanD = mean(D))
D_Grid_mean = D_Grid_mean[ , c(2:3,5)]
D_overlap = unique(D_Grid_mean)

# Join summary information and spatial data: 
goa.df = fortify(clip2_INPFC, region='id')
colnames(goa.df)[colnames(goa.df) == "id"] = "EEZgrid"
plot_PH_ATF_D = goa.df %>% right_join(D_overlap)
plot_PH_ATF_Dyr = goa.df %>% right_join(D_Grid)
################################################################
# Plot mean dietary overlap (all years combined):

### Dietary Overlap ###
textD = data.frame(text = c("Dietary Overlap"))

PH_ATF_Dplot = ggplot() +
  geom_polygon(data=plot_PH_ATF_D, aes(x=long, y=lat, group=group, fill=meanD), col="black", lwd=0.25) + 
  geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
  geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
  geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
  scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
  geom_polygon(data=INPFC_plot, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
  geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
  geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
  geom_text(data=textD, aes(label=text, x=-132.48, y=61.65, size=12), show.legend=FALSE) +
  geom_text(data=INPFC_cent, aes(group=StatArea, label=StatArea, x=START_LONGITUDE, y=START_LATITUDE, size=12), show.legend = FALSE) +
  geom_text(data=IPHC_2C, aes(label=text, x=-134.84, y=56.06, size=12), col="cadetblue1", fontface="bold", show.legend = FALSE) +
  geom_text(data=IPHC_3A, aes(label=text, x=-145.40, y=59.9, size=12), col="deepskyblue3", fontface="bold", show.legend = FALSE) +
  geom_text(data=IPHC_3B, aes(label=text, x=-157.39, y=55.65, size=12), col="blue", fontface="bold", show.legend = FALSE) +
  geom_text(data=IPHC_4A, aes(label=text, x=-167.75, y=53.05, size=12), col="midnightblue", fontface="bold", show.legend = FALSE) +
  theme_bw() +
  ggtitle("Dietary Overlap between Pacific Halibut and Arrowtooth Flounder") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(legend.position = c(0.957, 0.830)) +
  theme(plot.title = element_text(hjust = 0.5, size=14)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 1) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(6.0, "mm")) +
  theme(legend.key.height = unit(8.0, "mm")) +
  theme(legend.spacing.x = unit(2.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=12)) +
  theme(axis.text.x = element_text(family="Arial", size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  theme(axis.title.x = element_text(hjust=0.53, size=12)) +
  theme(axis.title.y = element_text(hjust=0.46, size=12)) +
  labs(x="Longitude", y="Latitude") +
  scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
  scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
  theme(legend.background = element_rect(fill="transparent"))

PH_ATF_Dplot
ggsave(filename="Plots/D_PH_ATF.png", plot=PH_ATF_Dplot, dpi=500, width=12, height=8, units="in")

# Plot dietary overlap by grid cell and year (for loop):
for(var in unique(plot_PH_ATF_Dyr$YEAR)) {
  p = ggplot() +
    geom_polygon(data=plot_PH_ATF_Dyr[plot_PH_ATF_Dyr$YEAR==var,], aes(x=long, y=lat, group=group, fill=D), col="black", lwd=0.25) + 
    geom_polygon(data=Area2C_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="cadetblue1", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="deepskyblue2", linetype="solid", lwd=1) +  
    geom_polygon(data=Area3B_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="blue", linetype="solid", lwd=1) +  
    geom_polygon(data=Area4A_df, alpha=0.3, aes(x=X, y=Y, group=factor(PID)), fill=NA, show.legend=FALSE, col="midnightblue", linetype="solid", lwd=1) +  
    scale_fill_distiller(palette="Spectral", name="", limits=c(0,1), breaks=c(0.00,0.25,0.50,0.75,1.00)) +
    geom_polygon(data=INPFC_plot, aes(x=X, y=Y, group=factor(PID)), fill=NA, col="black", linetype="solid", lwd=0.25) +
    geom_polygon(data=world2, aes(x=X, y=Y, group=factor(PID)), fill="lightgrey", col="black", lwd=0.25) +
    geom_polygon(data=Canada, aes(x=X, y=Y, group=factor(PID)), fill="gray91", col="black", lwd=0.25) +
    geom_text(data=textD, aes(label=text, x=-132.48, y=61.65, size=12), show.legend=FALSE) +
    geom_text(data=INPFC_cent, aes(group=StatArea, label=StatArea, x=START_LONGITUDE, y=START_LATITUDE, size=12), show.legend = FALSE) +
    geom_text(data=IPHC_2C, aes(label=text, x=-134.84, y=56.06, size=12), col="cadetblue1", fontface="bold", show.legend = FALSE) +
    geom_text(data=IPHC_3A, aes(label=text, x=-145.40, y=59.9, size=12), col="deepskyblue3", fontface="bold", show.legend = FALSE) +
    geom_text(data=IPHC_3B, aes(label=text, x=-157.39, y=55.65, size=12), col="blue", fontface="bold", show.legend = FALSE) +
    geom_text(data=IPHC_4A, aes(label=text, x=-167.75, y=53.05, size=12), col="midnightblue", fontface="bold", show.legend = FALSE) +
    theme_bw() +
    ggtitle("Dietary Overlap between Pacific Halibut and Arrowtooth Flounder") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
    theme(legend.position = c(0.957, 0.830)) +
    theme(plot.title = element_text(hjust = 0.5, size=14)) +
    theme(legend.title = element_blank()) +
    theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 1) +
    theme(legend.background = element_rect(fill="transparent")) + 
    theme(legend.key = element_blank()) +
    theme(legend.key.width = unit(6.0, "mm")) +
    theme(legend.key.height = unit(8.0, "mm")) +
    theme(legend.spacing.x = unit(2.0, "mm")) +
    theme(axis.text.y = element_text(family="Arial", size=12)) +
    theme(axis.text.x = element_text(family="Arial", size=12)) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(axis.title.x = element_text(hjust=0.53, size=12)) +
    theme(axis.title.y = element_text(hjust=0.46, size=12)) +
    labs(x="Longitude", y="Latitude") +
    scale_x_continuous(limits = c(lonmin, lonmax), expand = c(0,0)) +
    scale_y_continuous(limits = c(latmin, latmax), expand = c(0,0), labels = function(x) round(as.numeric(x), digits=0)) +
    theme(legend.background = element_rect(fill="transparent"))
  
  ggsave(p, filename=paste("Plots/D_PH_ATF_", var, ".png", sep=""), dpi=500, width=12, height=8, units="in")
}
#################################################################
### TEST FOR SPATIOTEMPORAL CHANGES IN DIETARY OVERLAP ###
#################################################################
# Run ANCOVAs to test for relationships among dietary overlap, survey year, and INPFC statistical area or IPHC regulatory area:

# Prepare overlap data:
dietary_Grid = unique(PH_ATF_spatialGAMred[,c("YEAR", "id2", "EEZgrid", "START_LONGITUDE", "START_LATITUDE")])
colnames(D_Grid) = c("YEAR", "id2", "EEZgrid", "D")
D_Grid$YEAR = as.factor(D_Grid$YEAR)
DietOverData = merge(D_Grid, dietary_Grid, all.x=TRUE)
DietOverData = na.omit(DietOverData)
coordinates(DietOverData) = ~ START_LONGITUDE + START_LATITUDE

# Join overlap and spatial data to identify INPFC/IPHC areas:
proj4string(DietOverData) = proj4string(INPFC_shape)
DietOverData$INPFC = over(DietOverData, INPFC_shape)
DietOverData$IPHC = over(DietOverData, IPHC_shape)
DietOverDf = as.data.frame(DietOverData)
DietOverDf$YEAR = as.factor(DietOverDf$YEAR)
DietOverDf = na.omit(DietOverDf)
#################################################################
# INPFC statistical areas:
DietINPFC = DietOverDf[,c(1:2,4,8)]
DietINPFCdf = unique(DietINPFC)
DietINPFCdf = subset(DietINPFCdf, INPFC.REP_AREA!="649")

# With interaction:
INPFCmodel_A = aov(D ~ YEAR * INPFC.REP_AREA, data=DietINPFCdf)
summary(INPFCmodel_A)
# Without (non-significant) interaction:
INPFCmodel_B = aov(D ~ YEAR + INPFC.REP_AREA, data=DietINPFCdf)
summary(INPFCmodel_B)
TukeyHSD(INPFCmodel_B, "YEAR")
plot(DietINPFCdf$D ~ DietINPFCdf$YEAR)
#################################################################
# IPHC regulatory areas:
DietIPHC = DietOverDf[,c(1:2,4,17)]
DietIPHCdf = unique(DietIPHC)
DietIPHCdf = subset(DietIPHCdf, IPHC.REG_AREA!="2B")

# With interaction:
IPHCmodel_A = aov(D ~ YEAR * IPHC.REG_AREA, data=DietIPHCdf)
summary(IPHCmodel_A)
# Without (non-significant) interaction:
IPHCmodel_B = aov(D ~ YEAR + IPHC.REG_AREA, data=DietIPHCdf)
summary(IPHCmodel_B)
TukeyHSD(IPHCmodel_B, "YEAR")
plot(DietIPHCdf$D ~ DietIPHCdf$Yr)
#################################################################
ggsave(filename="Plots/D_PH_ATF.png", plot=PH_ATF_Dplot, dpi=500, width=12, height=8, units="in")