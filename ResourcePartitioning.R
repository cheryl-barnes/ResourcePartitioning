# This script file combines spatial overlap and dietary overlap to quantify the correlation between the two measures and thus the degree of resource partitioning between Pacific Halibut and Arrowtooth Flounder in the Gulf of Alaska. Descriptions of resource partitioning can be found in Schoener (1974) and Ross (1986).

# References:
# Ross ST. Resource partitioning in fish assemblage: review of field studies. Copeia. 1986;1986(2):352–388.
# Schoener TW. Resource partitioning in ecological communities. Science. 1974;185(4145):27–39.

rm(list=ls())
graphics.off()

# test

setwd("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/ResourcePartitioning/ResourcePartitioning")

INPFCgrid = read.csv("clip_NMFS_100KM.csv")

setwd("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/")

spatial = read.csv("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/Spatial/PH_ATF_spatialGAMred_yr_max_100KM.csv")

colnames(spatial) = c("X", "UniqueID", "year", "id2", "EEZ_grid", "depth", "BT", "Long", "Lat", "PHpa.fit", "PHpa.se", "PHpa.low", "PHpa.hi", "PHcpueFit", "PHcpueSE", "PHcpueLow", "PHcpueHi", "PHpredBio", "PHstdBio", "ATFpa.fit", "ATFpa.se", "ATFpa.low", "ATFpa.hi", "ATFcpueFit", "ATFcpueSE", "ATFcpueLow", "ATFcpueHi", "ATFpredBio", "ATFstdBio", "S")

spatial = spatial[ , c(2:30)]
spatial$UniqueID = as.character(spatial$UniqueID)
spatial$EEZ_grid = as.character(spatial$EEZ_grid)

dietary = read.csv("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/Diet/DietaryOverlapGrid_30_69_100KM.csv")
dietary$UniqueID = as.character(paste(dietary$id2, dietary$Yr, sep="_"))
dietary = dietary[ , c(9:11, 14)]
colnames(dietary) = c("year", "id2", "D", "UniqueID")

require(dplyr)
overlap = spatial %>% right_join(dietary, by="UniqueID")

O1 = overlap %>%
  group_by(year) %>%
  summarize(meanO1 = mean(S)) 
O2 = overlap %>%
  group_by(year) %>%
  summarize(meanO2 = mean(D)) 
tracking = cbind(O1, O2)

# Remove unnecessary columns and delete enviro data to get one entry per grid cell:
overlap = overlap[ , c(1:6, 7:8, 29, 32)]
colnames(overlap) = c("UniqueID", "year", "id2", "EEZgrid", "depth", "BT", "lon", "lat", "S", "D")
overlap = overlap[complete.cases(overlap), ]
overlap = unique(overlap)




require(sp)
coordinates(overlap) = ~ lon + lat

setwd("~/Documents/UAF/Dissertation/Analyses/MappingInfo")

require(rgdal)
INPFC_shape = readOGR(".", "GOA_Shapes")
INPFC_shape = spTransform(INPFC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
proj4string(overlap) = proj4string(INPFC_shape)
overlap$INPFC = over(overlap, INPFC_shape)

IPHC_shape = readOGR(".", "GOA_Den")
IPHC_shape = spTransform(IPHC_shape, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
overlap$IPHC = over(overlap, IPHC_shape)

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

require(lattice)
xyplot(overlap$D ~ overlap$S)
xyplot(overlap$D ~ overlap$S | overlap$year)
xyplot(overlap$D ~ overlap$S | overlap$INPFC.REP_AREA)


setwd("~/Documents/UAF/Dissertation/Analyses/Ch_2_Manuscript/")
overlap$year = as.factor(overlap$year)

pdf("S_yr.pdf")
plot(x = overlap$year, y = overlap$S)
dev.off()

pdf("D_yr.pdf")
plot(x = overlap$year, y = overlap$D)
dev.off()

# write.csv(overlap, "overlap.csv")

require(reshape2)
overlap_red = overlap[ , c(1:10, 19)]

cor.test(overlap_red$D, overlap_red$S, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)

require(quantreg)
quantReg25D = rq(D ~ S, data=overlap, tau=0.25)
summary(quantReg25D, se="boot")

quantReg50D = rq(D ~ S, data=overlap, tau=0.5)
summary(quantReg50D, se="boot")

quantReg75D = rq(D ~ S, data=overlap, tau=0.75)
summary(quantReg75D, se="boot")

quantRegAll_D = rq(D ~ S, data=overlap, tau=c(0.25, 0.50, 0.75))
summary(quantRegAll_D, se="boot")

quantile.regressions_D = data.frame(t(coef(quantRegAll_D)))
colnames(quantile.regressions_D) = c("intercept", "slope")
quantile.regressions_D$quantile = rownames(quantile.regressions_D)
quantile.regressions_D

anova(quantReg25D, quantReg75D) # no difference - should not use quantile regression

quantRegPlot_D = rq(D ~ S, data=overlap, tau=seq(0.05, 0.95, by=0.05))
quantRegPlotSumm_D = summary(quantRegPlot_D)
plot(quantRegPlotSumm_D)

library(ggplot2)
scatterplot_D = qplot(x=S, y=D, data=overlap)
scatterplot_D + geom_abline(aes(intercept=intercept, slope=slope,
                                colour=quantile), data=quantile.regressions_D)

quantReg25S = rq(S ~ D, data=overlap, tau=0.25)
summary(quantReg25S, se="boot")

quantReg50S = rq(S ~ D, data=overlap, tau=0.5)
summary(quantReg50S, se="boot")

quantReg75S = rq(S ~ D, data=overlap, tau=0.75)
summary(quantReg75S, se="boot")

quantRegAll_S = rq(S ~ D, data=overlap, tau=c(0.25, 0.50, 0.75))
summary(quantRegAll_S, se="boot")

quantile.regressions_S = data.frame(t(coef(quantRegAll_S)))
colnames(quantile.regressions_S) = c("intercept", "slope")
quantile.regressions_S$quantile = rownames(quantile.regressions_S)
quantile.regressions_S

anova(quantReg25S, quantReg75S) # no difference - should not use quantile regression

quantRegPlot_S = rq(S ~ D, data=overlap, tau=seq(0.05, 0.95, by=0.05))
quantRegPlotSumm_S = summary(quantRegPlot_S)
plot(quantRegPlotSumm_S)

scatterplot_S = qplot(x=D, y=S, data=overlap)
scatterplot_S + geom_abline(aes(intercept=intercept, slope=slope,
                                colour=quantile), data=quantile.regressions_S)

overlap_long = melt(overlap_red, id.vars = c("UniqueID", "year", "id2", "EEZgrid", "lon", "lat", "INPFC.REP_AREA", "IPHC.REG_AREA"), measure.vars = 7:8, variable.name = "OverlapIndex", value.name = "measure")

summData = ddply(overlap_long, c("OverlapIndex", "year"), summarise,
                 N    = length(measure),
                 mean = mean(measure),
                 sd   = sd(measure),
                 se   = sd / sqrt(N))
summData
summData$year = as.numeric(as.character(summData$year))
pd = position_dodge(0.1) 
summData$CI = 1.96 * summData$sd

mean_over_line = ggplot(data=summData, aes(x=year, y=mean, col=OverlapIndex)) +
  geom_line() +
  geom_point() +
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
  labs(x="", y="Overlap Measure") +
  theme(axis.title.x = element_text(vjust=-0.13, size=12)) +
  theme(axis.title.y = element_text(vjust=1.1, size=12)) +
  theme(strip.text.x = element_text(family="Arial", size=12)) +
  theme(strip.text = element_text(hjust=0.5)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) 

mean_over_line
ggsave(filename="mean_over_line.png", plot=mean_over_line, dpi=500, width=6, height=5, units="in")

mean_over_line_CI = ggplot(data=summData, aes(x=year, y=mean, col=OverlapIndex)) +
  geom_line(aes(linetype=OverlapIndex)) +
  geom_point(aes(shape=OverlapIndex), size=2.5) +
  geom_errorbar(aes(ymin=mean-CI, ymax=mean+CI), width=0.5) +
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
  labs(x="", y="Overlap Measure") +
  theme(axis.title.x = element_text(vjust=-0.13, size=12)) +
  theme(axis.title.y = element_text(vjust=1.1, size=12)) +
  theme(strip.text.x = element_text(family="Arial", size=12)) +
  theme(strip.text = element_text(hjust=0.5)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) 

mean_over_line_CI
ggsave(filename="mean_over_line_CI.png", plot=mean_over_line_CI, dpi=500, width=6, height=5, units="in")

summData_INPFC = ddply(overlap_long, c("OverlapIndex", "year",  "INPFC.REP_AREA"), summarise,
                       N    = length(measure),
                       mean = mean(measure),
                       sd   = sd(measure),
                       se   = sd / sqrt(N))
summData_INPFC
summData_INPFC$year = as.numeric(as.character(summData_INPFC$year))
summData_INPFC$CI = 1.96 * summData_INPFC$sd
pd = position_dodge(0.1) 

mean_over_line_INPFCarea = ggplot(data=summData_INPFC, aes(x=year, y=mean, col=OverlapIndex)) +
  geom_errorbar(aes(ymin=mean-CI, ymax=mean+CI), width=0.5, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd) +
  scale_color_manual(values=c("blue", "red")) +
  facet_wrap(~ INPFC.REP_AREA, ncol = 1) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.0375, 0.965)) +
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
  labs(x="", y="Overlap Measure") +
  theme(axis.title.x = element_text(vjust=-0.13, size=12)) +
  theme(axis.title.y = element_text(vjust=1.1, size=12)) +
  theme(strip.text.x = element_text(family="Arial", size=12)) +
  theme(strip.text = element_text(hjust=0.5)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013)) 

mean_over_line_INPFCarea
ggsave(filename="mean_over_line_INPFCarea.png", plot=mean_over_line_INPFCarea, dpi=500, width=6, height=8, units="in")

summData_IPHC = ddply(overlap_long, c("OverlapIndex", "year",  "IPHC.REG_AREA"), summarise,
                      N    = length(measure),
                      mean = mean(measure),
                      sd   = sd(measure),
                      se   = sd / sqrt(N))
summData_IPHC
summData_IPHC$year = as.numeric(as.character(summData_IPHC$year))
summData_IPHC$CI = 1.96 * summData_IPHC$sd

mean_over_line_IPHCarea = ggplot(data=summData_IPHC, aes(x=year, y=mean, col=OverlapIndex)) +
  geom_errorbar(aes(ymin=mean-CI, ymax=mean+CI), width=0.5, position=pd, show.legend=FALSE) +
  geom_line(position=pd, show.legend=FALSE) +
  geom_point(position=pd, show.legend=FALSE) +
  scale_color_manual(values=c("blue", "red")) +
  facet_wrap(~ IPHC.REG_AREA, ncol = 1) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.0375, 0.965)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_blank()) +
  theme(legend.title.align = 0) +
  theme(legend.text = element_text(family="Arial", size=12), legend.text.align = 0) +
  theme(legend.background = element_rect(fill="transparent")) + 
  theme(legend.key = element_blank()) +
  theme(legend.key.width = unit(1.9, "mm")) +
  theme(legend.key.height = unit(6.0, "mm")) +
  theme(axis.text.y = element_text(family="Arial", size=12)) +
  theme(axis.text.x = element_text(family="Arial", size=12)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="", y="Overlap Measure") +
  theme(axis.title.x = element_text(vjust=-0.13, size=12)) +
  theme(axis.title.y = element_text(vjust=1.1, size=12)) +
  theme(strip.text.x = element_text(family="Arial", size=12)) +
  theme(strip.text = element_text(hjust=0.5)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(breaks = c(1990, 1993, 1996, 1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013))

mean_over_line_IPHCarea
ggsave(filename="mean_over_line_IPHCarea.png", plot=mean_over_line_IPHCarea, dpi=500, width=6, height=8, units="in")

summS = ddply(overlap, c("year", "INPFC.REP_AREA"), summarise,
              S_mean = mean(S))
summS

summD = ddply(overlap, c("year",  "INPFC.REP_AREA"), summarise,
              D_mean = mean(D))
summD

summOverlap = merge(summS, summD)

cor.test(overlap$S, overlap$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)

# Correlation coefficients:
Shum = subset(overlap, INPFC.REP_AREA == "Shumagin")
Chir = subset(overlap, INPFC.REP_AREA == "Chirikof")
Kod = subset(overlap, INPFC.REP_AREA == "Kodiak")
Yak = subset(overlap, INPFC.REP_AREA == "Yakutat")
SE = subset(overlap, INPFC.REP_AREA == "Southeastern")

cor.test(Shum$S, Shum$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(Chir$S, Chir$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(Kod$S, Kod$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(Yak$S, Yak$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(SE$S, SE$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)

# Correlation coefficients:
IPHC4A = subset(overlap, IPHC.REG_AREA == "4A")
IPHC3B = subset(overlap, IPHC.REG_AREA == "3B")
IPHC3A = subset(overlap, IPHC.REG_AREA == "3A")
IPHC2C = subset(overlap, IPHC.REG_AREA == "2C")

cor.test(IPHC4A$S, IPHC4A$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(IPHC3B$S, IPHC3B$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(IPHC3A$S, IPHC3A$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)
cor.test(IPHC2C$S, IPHC2C$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)

ResPart_overall_A = ggplot(data=overlap, aes(x=S, y=D, col=year)) +
  geom_point() +
  #geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  scale_color_manual(values=c("brown4", "red", "chocolate1", "gold","lemonchiffon", "chartreuse2", "darkgreen", "lightcyan1", "blue", "midnightblue", "mediumpurple1", "purple4", "black")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.950, 0.813)) +
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
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_overall_A
ggsave(filename="ResPart_overall_A_30_69.png", plot=ResPart_overall_A, dpi=500, width=8, height=8, units="in")

ResPart_overall_B = ggplot(data=overlap, aes(x=S, y=D, col=INPFC.REP_AREA)) +
  geom_point() +
  scale_color_manual(values=c("midnightblue", "blue", "deepskyblue2", "cadetblue1", "cadetblue1")) +
  #geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.905, 0.910)) +
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
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_overall_B
ggsave(filename="ResPart_overall_B.png_30_69.png", plot=ResPart_overall_B, dpi=500, width=8, height=8, units="in")

ResPart_overall_C = ggplot(data=overlap, aes(x=S, y=D, col=IPHC.REG_AREA, order=-as.numeric(IPHC.REG_AREA))) +
  geom_point() +
  #geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  scale_color_manual(values=c("midnightblue", "blue", "deepskyblue2", "cadetblue1")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.962, 0.927)) +
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
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_overall_C
ggsave(filename="ResPart_overall_C_30_69.png", plot=ResPart_overall_C, dpi=500, width=8, height=8, units="in")

ResPart_overall_D = ggplot(data=overlap, aes(x=S, y=D, col=DepthBin, order=-as.numeric(DepthBin))) +
  geom_point() +
  #geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  scale_color_manual(values=c("yellowgreen", "cadetblue1", "deepskyblue2", "blue", "midnightblue", "black")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.930, 0.887)) +
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
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))
ResPart_overall_D
ggsave(filename="ResPart_overall_D_30_69.png", plot=ResPart_overall_D, dpi=500, width=8, height=8, units="in")

ResPart_overall_E = ggplot(data=overlap, aes(x=S, y=D, col=TempBin, order=-as.numeric(TempBin))) +
  geom_point() +
  #geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  scale_color_manual(values=c("purple4", "blue2", "paleturquoise1", "green3", "yellow", "darkorange", "firebrick3")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.930, 0.887)) +
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
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_overall_E
ggsave(filename="ResPart_overall_E_30_69.png", plot=ResPart_overall_E, dpi=500, width=8, height=8, units="in")

ResPart_overall_F = ggplot(data=overlap, aes(x=S, y=D)) +
  geom_point() +
  geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="loess", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_overall_F
ggsave(filename="ResPart_overall_F_30_69.png", plot=ResPart_overall_F, dpi=500, width=8, height=, units="in")

ResPart_overall_G = ggplot(data=overlap, aes(x=S, y=D)) +
  geom_point() +
  #geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  scale_color_manual(values=c("black")) +
  theme_bw() +
  ggtitle("") +
  theme(legend.position = c(0.930, 0.887)) +
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
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_overall_G
ggsave(filename="ResPart_overall_G_30_69.png", plot=ResPart_overall_G, dpi=500, width=8, height=8, units="in")

ResPart_plot_yr = ggplot(data=overlap, aes(x=S, y=D, col=year)) +
  geom_point() +
  geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  facet_wrap(~ year, ncol=2, dir="h") +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.text.x = element_text(family="Arial", size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0,1), breaks=c(0.0,0.5,1.0)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0.0,0.5,1.0)) +
  theme(legend.background = element_rect(fill="transparent"))

ResPart_plot_yr
ggsave(filename="ResPart_plot_yr_30_69.png", plot=ResPart_plot_yr, dpi=500, width=12, height=8, units="in")

ResPart_plot_area = ggplot(data=overlap, aes(x=S, y=D, col=year), alpha=0.8) +
  geom_point() +
  geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  theme_bw() +
  ggtitle("") +
  facet_wrap(~ INPFC.REP_AREA, ncol=1, dir="v") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  scale_x_continuous(limits = c(0,1), breaks=c(0.0,0.5,1.0)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0.0,0.5,1.0)) +
  theme(legend.background = element_rect(fill="transparent"))

ResPart_plot_area
ggsave(filename="ResPart_plot_area_30_69.png", plot=ResPart_plot_area, dpi=500, width=4, height=8, units="in")

ResPart_plot_area_yr = ggplot(data=overlap, aes(x=S, y=D, col=year), alpha=0.8) +
  geom_point() +
  geom_smooth(data=overlap, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  theme_bw() +
  ggtitle("") +
  facet_grid(year ~ INPFC.REP_AREA) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  scale_x_continuous(limits = c(0,1), breaks=c(0.5)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0.5)) +
  theme(legend.background = element_rect(fill="transparent"))

ResPart_plot_area_yr
ggsave(filename="ResPart_plot_area_yr_30_69.png", plot=ResPart_plot_area_yr, dpi=500, width=8, height=8, units="in")

# Combine INPFC 620 and 630:
sub625 = subset(overlap, INPFC.REP_AREA!="610", drop=TRUE)
sub625 = subset(sub625, INPFC.REP_AREA!="640")
sub625 = subset(sub625, INPFC.REP_AREA!="650")

ResPart_plot_620_630 = ggplot(data=sub625, aes(x=S, y=D), alpha=0.8) +
  geom_point() +
  geom_smooth(data=sub625, aes(x=S, y=D), inherit.aes=FALSE, method="loess", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  scale_x_continuous(limits = c(0,1), breaks=c(0.5)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0.5)) +
  theme(legend.background = element_rect(fill="transparent"))

ResPart_plot_620_630

ResPart_plot_620_630_yr = ggplot(data=sub625, aes(x=S, y=D, col=year), alpha=0.8) +
  geom_point() +
  geom_smooth(data=sub625, aes(x=S, y=D), inherit.aes=FALSE, method="glm", se=TRUE, show.legend=FALSE, na.rm=TRUE) +
  theme_bw() +
  ggtitle("") +
  facet_wrap(~ year, ncol=2, dir="v") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.background = element_rect(fill="white", colour="black", size=1, line="solid")) +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_text(family="Arial", size=14)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  scale_x_continuous(limits = c(0,1), breaks=c(0.5)) +
  scale_y_continuous(limits = c(0,1), breaks=c(0.5)) +
  theme(legend.background = element_rect(fill="transparent"))

ResPart_plot_620_630_yr
ggsave(filename="ResPart_plot_620_630_yr_30_69.png", plot=ResPart_plot_620_630_yr, dpi=500, width=8, height=8, units="in")

##################################################################
# remove noisy eastern GOA
overlap_wGOA = subset(overlap, INPFC.REP_AREA!="Yakutat")
overlap_wGOA = subset(overlap_wGOA, INPFC.REP_AREA!="Southeastern")                  
View(overlap_wGOA)

ResPart_ByAreaYr = ggplot(data=overlap, aes(x=S, y=D, col=year)) +
  geom_point() +
  scale_color_manual(values=c("brown4", "red", "chocolate1", "gold","lemonchiffon", "chartreuse2", "darkgreen", "lightcyan1", "blue", "midnightblue", "mediumpurple1", "purple4", "black")) +
  theme_bw() +
  ggtitle("") +
  facet_wrap(~ INPFC.REP_AREA, ncol=1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_ByAreaYr
ggsave(filename="ResPart_ByAreaYr_30_69.png", plot=ResPart_ByAreaYr, dpi=500, width=5, height=5, units="in")

ResPart_ByArea = ggplot(data=overlap, aes(x=S, y=D)) +
  geom_point() +
  theme_bw() +
  ggtitle("") +
  facet_wrap(~ INPFC.REP_AREA, ncol=1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_ByArea
ggsave(filename="ResPart_ByArea_30_69.png", plot=ResPart_ByArea, dpi=500, width=4, height=8, units="in")

overlap_610_620 = subset(overlap_wGOA, INPFC.REP_AREA!="Shumagin")

ResPart_ByAreaPollYr_610_620 = ggplot(data=overlap_610_620, aes(x=S, y=D, col=year)) +
  geom_point() +
  scale_color_manual(values=c("green", "green", "gray", "black", "red", "red", "red", "red", "gray", "black", "green", "green", "black")) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.y = element_text(family="Arial", size=14)) +
  theme(axis.text.x = element_text(family="Arial", size=14)) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  labs(x="Spatial Overlap (S)", y="Dietary Overlap (D)") +
  theme(axis.title.x = element_text(vjust=-0.13, size=14)) +
  theme(axis.title.y = element_text(vjust=1.1, size=14)) +
  theme(strip.background = element_rect(fill="white", colour="white", size=1, line="solid")) +
  scale_x_continuous(limits = c(0, 1.03), expand = c(0,0)) +
  scale_y_continuous(limits = c(0, 1.03), expand = c(0,0))

ResPart_ByAreaPollYr_610_620 
ggsave(filename="ResPart_ByAreaPollYr_610_620_30_69.png", plot=ResPart_ByAreaPollYr_610_620 , dpi=500, width=6, height=5, units="in")

# Corr Test
cor.test(overlap_wGOA$S, overlap_wGOA$D, alternative="two.sided", method="pearson", exact=TRUE, conf.level = 0.95)

x = aggregate(D ~ year, data=overlap, FUN=mean)
aggregate(D ~ year, data=overlap, FUN=sd)

plot(x)

aggregate(S ~ IPHC.REG_AREA, data=overlap, FUN=mean)
aggregate(S ~ IPHC.REG_AREA, data=overlap, FUN=sd)
