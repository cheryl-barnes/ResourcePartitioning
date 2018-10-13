# This script file relies on results (e.g., spatial data frame manipulations, construction of uniform grid, overlap measure calculations) from the 'SpatialAnalyses' and 'DietaryAnalyses' files. Here, we combine spatial overlap and dietary overlap to quantify the correlation between the two measures and thus the degree of resource partitioning between Pacific Halibut and Arrowtooth Flounder in the Gulf of Alaska. Descriptions of resource partitioning can be found in Schoener (1974) and Ross (1986).

# Rcode developed by: Cheryl Barnes
# cheryl.barnes@alaska.edu
# Preliminary code (e.g., exploratory data analyses, sample size calculations, summary statistics, diagnostics) were excluded.

# References:
# Ross ST. Resource partitioning in fish assemblage: review of field studies. Copeia. 1986;1986(2):352–388.
# Schoener TW. Resource partitioning in ecological communities. Science. 1974;185(4145):27–39.

setwd("~/Documents/UAF/Dissertation/GitHub/ResourcePartitioning/")
# source(SpatialAnalyses.R)
# source(DietaryAnalyses.R)
################################################################
### DATA PREPARATION AND ANALYSES ###
################################################################
# Read in and prepare overlap data:
spatial = read.csv("Data/PH_ATF_S.csv")
dietary = read.csv("Data/PH_ATF_D.csv")
colnames(dietary) = c("X", "YEAR", "id2", "EEZgrid", "D")

# Join spatial and dietary overlap measures:
require(dplyr)
overlap = spatial %>% right_join(dietary, by=c("YEAR", "id2"))
overlap = overlap[complete.cases(overlap), ]
overlap = unique(overlap)
overlap = overlap[,c(2:7,10)]
colnames(overlap) = c("YEAR", "id2", "EEZgrid", "START_LONGITUDE", "START_LATITUDE", "S", "D")

# Test for basin-wide correlation between spatial overlap and dietary overlap (i.e., degree of resource partioning - proxy for competition):
cor.test(overlap$D, overlap$S, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)

# Test for region-specific correlations:
setwd("~/Documents/UAF/Dissertation/GitHub/ResourcePartitioning/Data/")

require(rgdal)
require(sp)
coordinates(overlap) = ~ START_LONGITUDE + START_LATITUDE
INPFC_shape = readOGR(".", "GOA_Shapes")
INPFC_shape = spTransform(INPFC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
proj4string(overlap) = proj4string(INPFC_shape)
overlap$INPFC = over(overlap, INPFC_shape)

IPHC_shape = readOGR(".", "GOA_Den")
IPHC_shape = spTransform(IPHC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
overlap$IPHC = over(overlap, IPHC_shape)

setwd("~/Documents/UAF/Dissertation/GitHub/ResourcePartitioning/")
overlap = as.data.frame(overlap)

require(ggplot2)
overlap = fortify(overlap)
overlap = overlap[complete.cases(overlap), ]
overlap = subset(overlap, INPFC.REP_AREA!="649")
overlap = subset(overlap, INPFC.REP_AREA!="659")
overlap = subset(overlap, IPHC.REG_AREA!="2A")
overlap = subset(overlap, IPHC.REG_AREA!="2B")
overlap = subset(overlap, IPHC.REG_AREA!="4B")

library(plyr)
levels(overlap$INPFC.REP_AREA) = c("Shumagin", "Chirikof", "Kodiak", "Yakutat", "NA", "Southeastern", "NA")
overlap$IPHC.REG_AREA = ordered(overlap$IPHC.REG_AREA, levels=c("4A", "3B", "3A", "2C"))

# INPFC statistical area coefficients:
Shum = subset(overlap, INPFC.REP_AREA == "Shumagin")
Chir = subset(overlap, INPFC.REP_AREA == "Chirikof")
Kod = subset(overlap, INPFC.REP_AREA == "Kodiak")
Yak = subset(overlap, INPFC.REP_AREA == "Yakutat")
SE = subset(overlap, INPFC.REP_AREA == "Southeastern")

cor.test(Shum$S, Shum$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(Chir$S, Chir$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(Kod$S, Kod$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(Yak$S, Yak$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(SE$S, SE$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)

# IPHC regulatory area coefficients:
IPHC4A = subset(overlap, IPHC.REG_AREA == "4A")
IPHC3B = subset(overlap, IPHC.REG_AREA == "3B")
IPHC3A = subset(overlap, IPHC.REG_AREA == "3A")
IPHC2C = subset(overlap, IPHC.REG_AREA == "2C")

cor.test(IPHC4A$S, IPHC4A$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(IPHC3B$S, IPHC3B$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(IPHC3A$S, IPHC3A$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
cor.test(IPHC2C$S, IPHC2C$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.90)
################################################################
# Plot relationship between spatial overlap and dietary overlap:
# Convert from wide to long format for plotting:
overlap_long = melt(overlap, id.vars = c("YEAR", "id2", "EEZgrid", "START_LONGITUDE", "START_LATITUDE"), measure.vars = 6:7, variable.name = "OverlapIndex", value.name = "measure")
library(ggplot2)
overlap_summ = ddply(overlap_long, c("OverlapIndex", "YEAR"), summarise,
                     N    = length(measure),
                     mean = mean(measure),
                     sd   = sd(measure),
                     se   = sd / sqrt(N))
overlap_summ
overlap_summ$YEAR = as.numeric(as.character(overlap_summ$YEAR))
pd = position_dodge(0.1) 
overlap_summ$CI = 1.96 * overlap_summ$sd

NicheOverlap = ggplot(data=overlap_summ, aes(x=YEAR, y=mean, col=OverlapIndex)) +
  geom_line(aes(linetype=OverlapIndex)) +
  geom_point(aes(shape=OverlapIndex), size=2.5) +
  geom_errorbar(aes(ymin=mean-CI, ymax=mean+CI), width=0.5) +
  scale_color_manual(values=c("blue", "red")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.05, 0.93)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(legend.spacing.x = unit(2.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=12)) +
  theme(axis.text.x = element_text(family="Arial", size=11)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Survey Year", y="Overlap Measure") +
  theme(axis.title.x = element_text(vjust=-0.13, size=12)) +
  theme(axis.title.y = element_text(vjust=1.1, size=12)) +
  theme(strip.text.x = element_text(family="Arial", size=12)) +
  theme(strip.text = element_text(hjust=0.5)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) +
  scale_y_continuous(limits=c(-0.5,1.0))

NicheOverlap
ggsave(filename="Plots/NicheOverlap.png", plot=NicheOverlap, dpi=500, width=6, height=5, units="in")

plot(overlap$D ~ overlap$S)