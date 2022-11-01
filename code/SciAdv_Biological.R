######################################################################################################
######################################################################################################
library(plyr); library(dplyr);library(broom);library(ggplot2); library(lubridate);library(LoLinR); library(stringr)
library(lme4); library(lmerTest); library(multcomp); library(phytotools); library(googledrive); library(Rmisc)
library(tibble); library(ggpubr); library(wesanderson); library(tidyverse);library(vegan); library(ggfortify);
library(lsmeans); library(RLRsim); library(remef); library(nlme);
library(timeDate)
######################################################################################################
######################################################################################################
rm(list = ls())

# This is going to read in weight data: 
datUrchin = read.csv('data/Winter2021Wts.csv')
gp2 <- wes_palettes$Zissou1
datUrchin$Site <- factor(datUrchin$Site, levels = c("HA","PV","WP","NO","PA","VD"))
datUrchin$POPTrt <- paste(datUrchin$POP, datUrchin$Trt)

#Calculating growth metrics: relative growth rate
datUrchin$RGR2 <- log(datUrchin$WetWtF/datUrchin$WetWtI)*100 #growth
datUrchin$RGRB2 <- log(datUrchin$BuoyWtF/datUrchin$BuoyWtI)*100 #calcification Buoy = Buoyant weight

# Only select urchins that survived until the end of the experiment
datUrchinHealthy = datUrchin %>%
  filter(Disease == "N") 

#Level so that southern CA sites are grouped together
datUrchinHealthy$factor <- factor(datUrchinHealthy$Site, levels = c("HA","PV","WP","NO","PA","VD"))
datUrchinHealthy <- datUrchinHealthy[complete.cases(datUrchinHealthy[ ,14]),]

#Log transform covariate to linearize
datUrchinHealthy$logWI = log(datUrchinHealthy$WetWtI)  
datUrchinHealthy$logBI = log(datUrchinHealthy$BuoyWtI)  
datUrchinHealthy$logWF = log(datUrchinHealthy$WetWtF)  

#Subset just the "local adaptation" experiment data
datUrchinHealthy_current <- subset(datUrchinHealthy, datUrchinHealthy$CF == 'C')

#Mixed model testing differences in growth across treatments and populations
GrowthMdl <- lmer(RGR2 ~ POP * Trt + logWI + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_current)
summary(GrowthMdl)
anova(GrowthMdl)

#Calc sample sizes for trts
Sum_count <- datUrchinHealthy %>%
  count(POP, Trt, sort = TRUE) 

#Compare means at y-int
GMdl <- emmeans(GrowthMdl, ~ Trt*POP)

#Create custom contrasts
GNNC <- c(1,0,0,0)
GNSC <- c(0,1,0,0)
GSNC <- c(0,0,1,0)
GSSC <- c(0,0,0,1)

#Pairwise contrasts of interest for y-int
contrast(GMdl, method = list(GSSC - GNSC))
contrast(GMdl, method = list(GSNC - GNNC))
contrast(GMdl, method = list((GSNC - GNNC)-(GSSC - GNSC)))

#Extract growth at mean covariate
Growth <- data.frame(summary(GMdl)$emmean)
colnames(Growth) <- 'Growth'
GrowthSE <- data.frame(summary(GMdl)$SE)
colnames(GrowthSE) <- 'GrowthSE'
Growth <- cbind(Growth, GrowthSE)
Growth$Trt <- c('NC','SC','NC','SC')
Growth$POP <- c('N','N','S','S')

#Plot of mean growth
iG = ggplot(Growth, aes(x = Trt, y = Growth, fill = POP, group = POP)) +
  geom_line(color = "black", position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Growth-GrowthSE, ymax = Growth+GrowthSE), position=position_dodge(width=0.5), 
                color = 'black', width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  annotate("text", x = .875, y = 12.22+1.24+.4, label = "(52)") +
  annotate("text", x = 1.125, y = 16.43+1.47+.4, label = "(32)") +
  annotate("text", x = 1.875, y = 3.26+1.32+.4, label = "(44)") +
  annotate("text", x = 2.125, y = 14.89+1.63+.4, label = "(25)") +
  geom_point(shape = 21, size = 5, position=position_dodge(width=0.5)) +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
       legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x=element_blank()) +
  ylab(expression(atop("Growth",paste("(Log(",frac(~ww[F],~ww[I]),")*100)")))) +
  scale_x_discrete(limits = c('NC', 'SC'), labels = c('Strong Upwelling', 'Weak Upwelling')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
iG

#Subset data for just "climate change" experiment
datUrchinHealthy_CF <- subset(datUrchinHealthy, datUrchinHealthy$POPTrt == 'NORTH NC' | datUrchinHealthy$POPTrt == 'NORTH NF' |
                                datUrchinHealthy$POPTrt == 'SOUTH SC' | datUrchinHealthy$POPTrt == 'SOUTH SF')

#Mixed model testing differences in growth across treatments and populations
GrowthMdlCF <- lmer(RGR2 ~ POP * CF + logWI + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_CF)
summary(GrowthMdlCF)
anova(GrowthMdlCF)

GMdlF <- emmeans(GrowthMdlCF, ~ CF*POP)

GrowthCF <- data.frame(summary(GMdlF)$emmean)
colnames(GrowthCF) <- 'Growth'
GrowthCFSE <- data.frame(summary(GMdlF)$SE)
colnames(GrowthCFSE) <- 'GrowthSE'

GrowthCF <- cbind(GrowthCF, GrowthCFSE)
GrowthCF$Trt <- c('C','F','C','F')
GrowthCF$POP <- c('N','N','S','S')

#Plots of slope and intercepts
j = ggplot(GrowthCF, aes(x = Trt, y = Growth, fill = POP, group = POP)) +
  geom_errorbar(aes(ymin = Growth-GrowthSE, ymax = Growth+GrowthSE), position=position_dodge(width=0.5),
                color = 'black', width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5, position=position_dodge(width=0.5)) +
  annotate("text", x = .875, y = 13.91+2.74+.4, label = "(52)") +
  annotate("text", x = 1.125, y = 17.05+2.94+.4, label = "(25)") +
  annotate("text", x = 1.875, y = 6.88+2.78+.4, label = "(44)") +
  annotate("text", x = 2.125, y = 6.69+3.17+.4, label = "(15)") +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x=element_blank()) +
  ylab(expression(atop("Growth",paste("(Log(",frac(~ww[F],~ww[I]),")*100)")))) +
  scale_x_discrete(labels = c('Current', 'Future')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
j

############################################################################################################################
############################################################################################################################
#Mixed model testing differences in calcification across treatments and populations
CalcMdl <- lmer(RGRB2 ~ POP * Trt + logBI + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_current)
summary(CalcMdl)
anova(CalcMdl)

#Compare means at y-int
CMdl <- emmeans(CalcMdl, ~ Trt*POP)

#Pairwise contrasts of interest for y-int
contrast(CMdl, method = list(GSNC - GNNC))
contrast(CMdl, method = list(GSSC - GNSC))
contrast(CMdl, method = list((GSNC - GNNC)-(GSSC - GNSC)))

CalcP <- data.frame(summary(CMdl)$emmean)
colnames(CalcP) <- 'Calc'
CalcSE <- data.frame(summary(CMdl)$SE)
colnames(CalcSE) <- 'CalcSE'

CalcP <- cbind(CalcP, CalcSE)
CalcP$Trt <- c('NC','SC','NC','SC')
CalcP$POP <- c('N','N','S','S')

#Plots of slope and intercepts
k = ggplot(CalcP, aes(x = Trt, y = Calc,fill = POP, group = POP)) +
  geom_line(color = "black",position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Calc-CalcSE, ymax = Calc+CalcSE), position=position_dodge(width=0.5),
                color = 'black', width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5,position=position_dodge(width=0.5)) +
  annotate("text", x = .875, y = 5.8454+2.22+.5, label = "(52)") +
  annotate("text", x = 1.125, y = 16.3319+2.35+.5, label = "(32)") +
  annotate("text", x = 1.875, y = -.0898+2.26+.5, label = "(44)") +
  annotate("text", x = 2.125, y = 13.353+2.45+.5, label = "(25)") +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x="Rearing Environment") +
  ylab(expression(atop("Calcification",paste("(Log(",frac(~bw[F],~bw[I]),")*100)")))) +
  scale_x_discrete(limits = c('NC', 'SC'), labels = c('Strong Upwelling', 'Weak Upwelling')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
k

#Mixed model testing differences in growth across treatments and populations
CalcMdlCF <- lmer(RGRB2 ~ POP * CF + logBI + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_CF)
summary(CalcMdlCF)
anova(CalcMdlCF)

#Compare means at y-int
CMdlCF <- emmeans(CalcMdlCF, ~ CF*POP)

CalcCF <- data.frame(summary(CMdlCF)$emmean)
colnames(CalcCF) <- 'Growth'
CalcCFSE <- data.frame(summary(CMdlCF)$SE)
colnames(CalcCFSE) <- 'GrowthSE'

CalcCF <- cbind(CalcCF, CalcCFSE)
CalcCF$Trt <- c('C','F','C','F')
CalcCF$POP <- c('N','N','S','S')

#Plots of slope and intercepts
l = ggplot(CalcCF, aes(x = Trt, y = Growth, fill = POP, group = POP))+
  geom_errorbar(aes(ymin = Growth-GrowthSE, ymax = Growth+GrowthSE), position=position_dodge(width=0.5), 
                color = 'black', width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5,position=position_dodge(width=0.5)) +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  annotate("text", x = .875, y = 7.54+3.33+.5, label = "(52)") +
  annotate("text", x = 1.125, y = 15.21+3.51+.5, label = "(25)") +
  annotate("text", x = 1.875, y = 1.43+3.36+.5, label = "(44)") +
  annotate("text", x = 2.125, y = 5.75+3.73+.5, label = "(15)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x="Treatment") +
  ylab(expression(atop("Calcification",paste("(Log(",frac(~bw[F],~bw[I]),")*100)")))) +
  scale_x_discrete(labels = c('Current', 'Future')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
l

######################################################################################################
######################################################################################################
#Biomass correcting grazing data and changing all negative values to zero since -growth of kelp is not possible
datUrchinHealthy$GrazingCor <- datUrchinHealthy$Grazing/datUrchinHealthy$WetWtF
datUrchinHealthy$GrazingCor <- ifelse(datUrchinHealthy$GrazingCor < 0, 0, datUrchinHealthy$GrazingCor)
grazing <- datUrchinHealthy[complete.cases(datUrchinHealthy$GrazingCor), ] 

#calculate sample sizes for treatments
Sum_count <- grazing %>%
  count(POP, Trt, sort = TRUE) 

#Calculating the mean and adding 10% of mean to values in order to include zeros
GMean <- mean(grazing$GrazingCor)
Add <- GMean*.1
grazing$logGCor <- (log(grazing$GrazingCor+Add))

grazMod = lm(logGCor ~ WetWtF, data = grazing)
summary(grazMod)

ij = ggplot(grazing, aes(x = WetWtF, y = logGCor))+
  geom_point() + theme_bw() + theme(axis.text = element_text(size = 6)) +
  labs(x="Wet Weight") +
  ylab(bquote('Grazing Rate (g kelp ww' ~hr^-1*'g urch'~ww^-1*')'))
ij

grazing_current <- subset(grazing, grazing$CF == 'C')

#Mixed model testing differences in grazing across treatments and populations
GrazeMdl <- lmer(logGCor ~ Trt * POP + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = grazing_current)
summary(GrazeMdl)
anova(GrazeMdl)

Graze <- emmeans(GrazeMdl, ~ POP * Trt)

#Create custom contrasts
GNC <- c(1,1,0,0)
GSC <- c(0,0,1,1)

#Pairwise contrasts of interest
contrast(Graze, method = list(GNC - GSC))

GrazeC <- data.frame(summary(Graze)$emmean)
colnames(GrazeC) <- 'Graze'
GrazeCSE <- data.frame(summary(Graze)$SE)
colnames(GrazeCSE) <- 'GrazeSE'

GrazeC <- cbind(GrazeC, GrazeCSE)
GrazeC$Trt <- c('NC','SC','NC','SC')
GrazeC$POP <- c('N','N','S','S')

#Plots of slope and intercepts
m = ggplot(GrazeC, aes(x = Trt, y = Graze, fill = POP, group = POP))+
  geom_line(color = "black", position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Graze-GrazeSE, ymax = Graze+GrazeSE), position=position_dodge(width=0.5), color = 'black', width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5, position=position_dodge(width=0.5)) +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  annotate("text", x = .875, y = -3.49+.157+.03, label = "(52)") +
  annotate("text", x = 1.125, y = -3.26+.166+.03, label = "(32)") +
  annotate("text", x = 1.875, y = -3.32+.186+.035, label = "(44)") +
  annotate("text", x = 2.125, y = -2.63+.215+.03, label = "(22)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x= 'Rearing Environment') +
  ylab(bquote(atop('Grazing Rate','(g kelp ww' ~hr^-1*'g urch'~ww^-1*')'))) +
  scale_x_discrete(limits = c('NC', 'SC'), labels = c('Strong Upwelling', 'Weak Upwelling')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
m

grazing_CF <- subset(grazing, grazing$POPTrt == 'NORTH NC' | grazing$POPTrt == 'NORTH NF' |
                       grazing$POPTrt == 'SOUTH SC' | grazing$POPTrt == 'SOUTH SF')
#Mixed model testing differences in growth across treatments and populations
GrazeMdlF <- lmer(logGCor ~ POP * CF + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = grazing_CF)
summary(GrazeMdlF)
anova(GrazeMdlF)

GrazeF <- emmeans(GrazeMdlF, ~ POP*CF)

GrazeCF <- data.frame(summary(GrazeF)$emmean)
colnames(GrazeCF) <- 'Growth'
GrazeCFSE <- data.frame(summary(GrazeF)$SE)
colnames(GrazeCFSE) <- 'GrowthSE'

GrazeCF <- cbind(GrazeCF, GrazeCFSE)
GrazeCF$POP <- c('N','S','N','S')
GrazeCF$Trt <- c('C','C','F','F')

#Plots of slope and intercepts
n = ggplot(GrazeCF, aes(x = Trt, y = Growth, fill = POP))+
  geom_errorbar(aes(ymin = Growth-GrowthSE, ymax = Growth+GrowthSE), position=position_dodge(width=0.5), width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5,position=position_dodge(width=0.5)) +
  annotate("text", x = .875, y = -3.49+.14+.04, label = "(52)") +
  annotate("text", x = 1.125, y = -2.63+.205+.035, label = "(22)") +
  annotate("text", x = 1.875, y = -3.64+.152+.04, label = "(44)") +
  annotate("text", x = 2.125, y = -3.05+0.245+.04, label = "(15)") +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x="Treatment") +
  ylab(bquote(atop('Grazing Rate','(g kelp ww' ~hr^-1*'g urch'~ww^-1*')'))) +
  scale_x_discrete(labels = c('Current', 'Future'))
n

######################################################################################################
######################################################################################################
#Moving on to respirometry data
Resp <- datUrchin[!is.na(datUrchin$RUN),] #Remove any rows that doesn't have respirometry data
TIME8 = as.numeric((hm(Resp$Time)))
dat_long = data.frame(Vial = NA,Sample.ID=NA, Time = NA, Oxygen = NA, Temp = NA, RUN = NA, TRT = NA, CS = NA, POP = NA, QC = NA)  #Initialize empty dataframe

for (b in seq(from=21, to=41, by=3)){ #rbind the proper columns, so column 4 (ID) and the others. 
  # This takes the time from the column, and changes it to "seconds since capped"
  Resp[,b] = as.numeric((hm(Resp[,b])))-TIME8
  temp = Resp[,c(18,7,((b):(b+2)),19, 6, 20, 2, 42)] # This is how to do 4,5,6,7,8 -- 4,9,10,11,12
  # This adds a "vial" column - which you can ignore if you already have one. 
  colnames(temp) = colnames(dat_long) #column names for rbind
  dat_long = rbind(dat_long, temp)
}

dat_long = dat_long[complete.cases(dat_long$Oxygen),]
dat_long$Oxygen = as.numeric(dat_long$Oxygen)

dat_long_edit <- dat_long %>%
  filter(QC == "G") 

ggplot(dat_long, aes(x = Time, y = Oxygen, color = factor(RUN))) + geom_point() + 
  geom_line(aes(group = paste(Sample.ID,RUN)))

# plots all controls and all samples and separates by run, remember that controls should be 
# fairly flat lined, whereas samples should be going down
p <- ggplot(dat_long, aes(x = Time, y = Oxygen, color = POP)) + 
  geom_point(size = 1) +
  theme_classic() +
  geom_line(aes(group = Sample.ID), alpha = 1/2) +
  scale_color_manual(values = wes_palette("Zissou1", 2, "continuous")) +
  facet_grid(cols = vars(CS), rows = vars(TRT)) +
  theme(legend.position = "none") 
p

dat_long$combined <- paste(dat_long$Sample.ID,dat_long$RUN)
u <-data.frame(unique(dat_long$combined))
f = data.frame(Slope = NA,Sig=NA, Temp = NA, Sample.ID = NA, RUN = NA)  #Initialize empty dataframe

for (i in 1:nrow(u)) { 
  x <- dplyr::select(filter(dat_long, combined == u[i,]),c(Vial, Sample.ID, Time, Oxygen, Temp, RUN))
  mean_temp = mean(x$Temp, na.rm = TRUE)
  regression = rankLocReg(xall = x$Time, yall = x$Oxygen, alpha = .7)
  summary(regression)
  a = regression$allRegs
  f[i,] <- c(a[1,]$b1 * -1,a[1,]$b1LoCI > 0 | a[1,]$b1UpCI<0, mean_temp, x$Sample.ID[1], x$RUN[1])
}

dat_sum1 <- merge(Resp[,1:53],f, by=c("Sample.ID"))

#Create columns for calculating out rates
dat_sum1$resp_mg_sec <- NA
dat_sum1$resp_g_sec <- NA
dat_sum1$resp_mol_sec <- NA
dat_sum1$resp_umol_min <- NA

dat_sum1$Slope = as.numeric(dat_sum1$Slope)

#Filter for controls only
Controls1 = dat_sum1 %>%
  filter(CS == "C") 

#Calc mean of controls
Con_Sum1 <- summarySE(data=Controls1, measurevar = "Slope",
                      groupvars = c("RUN.x"), na.rm = TRUE)

#Subset only samples
Samp1 = dat_sum1 %>%
  filter(CS == "S")

#Adds run specific ave control to each sample
Samp1 <- merge(Con_Sum1, Samp1, by.x = c("RUN.x"), by.y = c("RUN.x"))
#Correct slope of samples by subracting the ave slope of controls in a run
Samp1$Slope.y = Samp1$Slope.y-Samp1$Slope.x

# Volume of vials used; S = small, L = large
Svial_vol_ml = 99.957
Lvial_vol_ml = 229.7059

# Slopes are in mg/L per second. convert to something else: ####
# firstly multiply by volume of vial in liters to convert to mg per second: 
for (i in 1:nrow(Samp1)) {
  if(Samp1$VTYPE[i]=='S') {
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Svial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }else{
    Samp1$resp_mg_sec[i] = Samp1$Slope.y[i] * (Lvial_vol_ml/1000)
    Samp1$resp_g_sec[i] = Samp1$resp_mg_sec[i] / 1000
    Samp1$resp_mol_sec[i] = Samp1$resp_g_sec[i] / 31.998
    Samp1$resp_umol_min[i] = Samp1$resp_mol_sec[i] * 1000000 * 60
  }
}   

dat_sum1 <- Samp1

dat_sum1$resp_corr_2 = dat_sum1$resp_umol_min/dat_sum1$WetWtF

dat_sum1$log_resp = log(dat_sum1$resp_umol_min)
dat_sum1$log_massWW = log(dat_sum1$WetWtF)

dat_Resp = dat_sum1 %>%
  filter(Problem == "G")

#calculate sample sizes for treatments
Sum_count <- dat_Resp %>%
  count(POP, Trt, sort = TRUE) 

###Single mass correction for all data#####################################################################
lmer_resp_corr_2 = lm(log_resp ~ log_massWW, data = dat_Resp)
summary(lmer_resp_corr_2)
coef(lmer_resp_corr_2)
coef = 0.4860836
# 
mean_mass = mean(dat_Resp$WetWtF)
# 
dat_Resp$resp_corr_2 = (dat_Resp$resp_umol_min/dat_Resp$WetWtF) *
  (dat_Resp$WetWtF/mean_mass)^(1-coef)
############################################################################################################

ggplot(dat_Resp, aes(x = WetWtF, y = resp_corr_2))+geom_point()+geom_smooth(method = 'lm') +
  theme_bw()+
  labs(x = 'Wet weight (g)', y = 'respiration rate (corr), umol/hr')

dat_Resp$Site <- factor(dat_Resp$Site, levels = c("HA", "PV", "WP","PA","VD","NO"))
q <- ggplot(dat_Resp, aes(x = Trt, y = resp_corr_2)) +
  geom_boxplot() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "bottom") +
  scale_fill_brewer(palette = "Dark2") +
  annotate(geom="text", x = 1, y = 0.01765, label = "a", color = "black") +
  annotate(geom="text", x = 2, y = 0.014, label = "a", color = "black") +
  annotate(geom="text", x = 3, y = 0.021, label = "b", color = "black") +
  annotate(geom="text", x = 4, y = 0.0206, label = "b", color = "black") +
  scale_x_discrete(labels = c('North Current', 'North Future', 'South Current', 'South Future')) +
  labs(y = 'mass specific respiration rate (umol/hr/g)', x = 'Treatment')
q

dat_RespC <- subset(dat_Resp, dat_Resp$Trt == 'NC' | dat_Resp$Trt == 'SC')

RMdl <- lmer(resp_corr_2 ~ Trt * POP + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = dat_RespC)
summary(RMdl)
anova(RMdl)

Rm <- emmeans(RMdl, ~ Trt*POP)

RespC <- data.frame(summary(Rm)$emmean)
colnames(RespC) <- 'Growth'
RespCSE <- data.frame(summary(Rm)$SE)
colnames(RespCSE) <- 'GrowthSE'

RespC <- cbind(RespC, RespCSE)
RespC$POP <- c('N','N','S','S')
RespC$Trt <- c('NC','SC','NC','SC')

#Plots of slope and intercepts
o = ggplot(RespC, aes(x = Trt, y = Growth, fill = POP, group = POP))+
  geom_line(color = "black",position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Growth-GrowthSE, ymax = Growth+GrowthSE), position=position_dodge(width=0.5),
                width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5,position=position_dodge(width=0.5)) +
  annotate("text", x = .875, y = 0.01189+.000841+0.00021, label = "(31)") +
  annotate("text", x = 1.125, y = .00899+.000885+0.00021, label = "(23)") +
  annotate("text", x = 1.875, y = .01456+.000859+0.00021, label = "(28)") +
  annotate("text", x = 2.125, y = .01335+0.000955+0.00021, label = "(16)") +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x= element_blank()) +
  ylab(bquote(atop('Respiration Rate', '('*mu~'mol' ~hr^-1 ~ind^-1*')')))+ 
  scale_x_discrete(limits = c('NC', 'SC'), labels = c('Strong Upwelling', 'Weak Upwelling')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
o

dat_RespCF <- subset(dat_Resp, dat_Resp$POPTrt == 'NORTH NC' | dat_Resp$POPTrt == 'NORTH NF' |
                       dat_Resp$POPTrt == 'SOUTH SC' | dat_Resp$POPTrt == 'SOUTH SF')
#Mixed model testing differences in growth across treatments and populations
RMdlCF <- lmer(resp_corr_2 ~ POP * CF + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = dat_RespCF)
summary(RMdlCF)
anova(RMdlCF)

######################################################################################################
######################################################################################################
Mortality <- read.csv('data/Mortality.csv')
Mortality$SITE <- factor(Mortality$SITE, levels = c("HA", "PV", "WP","PA","VD","NO"))

Sum_count <- datUrchinHealthy %>%
  count(POP, Trt, sort = TRUE) 

#Only current treatments
Mortality$C <- str_sub(Mortality$TRT, 2, 2)
Mortality_Current <- subset(Mortality, Mortality$C == 'C')

MortMdl <- lmer(Deaths ~ TIMESTEP * TRT * POP + AveInitialWt + (1|POP:SITE) , data = Mortality_Current)
summary(MortMdl)
anova(MortMdl)

MMdlC <- emtrends(MortMdl, ~TRT*POP, var = "TIMESTEP")

#Create custom contrasts
MNNC <- c(1,0,0,0)
MNSC <- c(0,1,0,0)
MSNC <- c(0,0,1,0)
MSSC <- c(0,0,0,1)

#Pairwise contrasts of interest
contrast(MMdlC, method = list(MNNC - MSNC))
contrast(MMdlC, method = list(MNSC - MSSC))

MortC <- data.frame(summary(MMdlC)$TIMESTEP.trend)
colnames(MortC) <- 'Mort'
MortCSE <- data.frame(summary(MMdlC)$SE)
colnames(MortCSE) <- 'MortSE'

MortC <- cbind(MortC, MortCSE)
MortC$Trt <- c('NC','SC','NC','SC')
MortC$POP <- c('N','N','S','S')

#Plots of slope and intercepts
p = ggplot(MortC, aes(x = Trt, y = Mort, fill = POP, group = POP)) +
  geom_line(color = "black",position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Mort-MortSE, ymax = Mort+MortSE), position=position_dodge(width=0.5), color = "black", width = 0.2)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5, position=position_dodge(width=0.5)) +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  annotate("text", x = .875, y = 0.000036+.00112+0.001, label = "(52)") +
  annotate("text", x = 1.125, y = .004491+.0012+0.004, label = "(35)") +
  annotate("text", x = 1.875, y = .02594+.00112+0.001, label = "(52)") +
  annotate("text", x = 2.125, y = .016816+0.00124+0.001, label = "(34)") +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  ylab(expression(atop("Mortality",paste("(",Delta,frac(italic("N"),Time),")")))) +
  labs(x=element_blank()) +
  scale_x_discrete(limits = c('NC', 'SC'), labels = c('Strong Upwelling', 'Weak Upwelling')) +
  guides(shape = guide_legend(title="Population", override.aes = list(size = 2)))
p

Mortality$POPTrt <- str_c(Mortality$POP,Mortality$TRT)

Mortality_CF <- subset(Mortality, Mortality$POPTrt == 'NNC' | Mortality$POPTrt == 'NNF' |
                         Mortality$POPTrt == 'SSC' | Mortality$POPTrt == 'SSF')

MortMdlC <- lmer(Deaths ~ TIMESTEP * C * POP + AveInitialWt + (1|POP:SITE) , data = Mortality_CF)
summary(MortMdlC)
anova(MortMdlC)

MMdlCF <- emtrends(MortMdlC, ~C*POP, var = "TIMESTEP")

#Create custom contrasts
MNC <- c(1,0,0,0)
MNF <- c(0,1,0,0)
MSC <- c(0,0,1,0)
MSF <- c(0,0,0,1)

#Pairwise contrasts of interest
contrast(MMdlCF, method = list(MNC - MNF))
contrast(MMdlCF, method = list(MSC - MSF))
contrast(MMdlCF, method = list((MNC-MNF)-(MSC - MSF)))

MortCF <- data.frame(summary(MMdlCF)$TIMESTEP.trend)
colnames(MortCF) <- 'Mort'
MortCFSE <- data.frame(summary(MMdlCF)$SE)
colnames(MortCFSE) <- 'MortSE'

MortCF <- cbind(MortCF, MortCFSE)
MortCF$Trt <- c('C','F','C','F')
MortCF$POP <- c('N','N','S','S')

#Plots of slope and intercepts
q = ggplot(MortCF, aes(x = Trt, y = Mort, fill = POP, group = POP)) +
  geom_errorbar(aes(ymin = Mort-MortSE, ymax = Mort+MortSE), position=position_dodge(width=0.5), width = 0.2)+
  #geom_line(color = "black") +
  annotate("text", x = .875, y = 0.0000924+.00119+0.002, label = "(52)") +
  annotate("text", x = 1.125, y = .0123+.00129+0.002, label = "(34)") +
  annotate("text", x = 1.875, y = .012+.00119+0.002, label = "(52)") +
  annotate("text", x = 2.125, y = .043+0.0015+0.002, label = "(35)") +
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(shape = 21, size = 5, position=position_dodge(width=0.5)) +
  scale_fill_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  guides(fill=guide_legend(title="Population")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = "top",
        legend.box.background = element_rect(colour = "black", size = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  ylab(expression(atop("Mortality",paste("(",Delta,frac(italic("N"),Time),")")))) +
  labs(x=element_blank()) +
  scale_x_discrete(limits = c('C', 'F'), labels = c('Current', 'Future'))
q

######################################################################################################
######################################################################################################
#Look first at initial gonad:SOM ratios
# Only select urchins that did not have outward signs of disease
datUrchinInitial = datUrchin %>%
  filter(Trt == "I") 
datUrchinInitial$Gonad <- datUrchinInitial$Ganimal - datUrchinInitial$G_AFDW
datUrchinInitial$SOM <- datUrchinInitial$Sanimal - datUrchinInitial$S_AFDW
datUrchinInitial$Gonad <- ifelse(datUrchinInitial$Gonad < 0, 0, datUrchinInitial$Gonad)
datUrchinInitial$GSOM <- datUrchinInitial$Gonad/datUrchinInitial$SOM

datUrchinInitial = datUrchinInitial %>%
  filter(GSOM >= 0) 

Sum_count <- datUrchinInitial %>%
  count(POP, sort = TRUE) 

datUrchinInitial$POP <- factor(datUrchinInitial$POP, levels = c("NORTH", "SOUTH"))

r <- ggplot(datUrchinInitial, aes(x = WetWtI, y = GSOM, color = POP)) +
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(size = 5, alpha = 1/5, aes(color = POP)) +
  guides(fill = guide_legend(title = 'Population'), colour=guide_legend(title="Population")) +
  scale_color_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  geom_smooth(method = "lm", aes(fill = POP, color = POP)) +
  scale_fill_brewer(palette = "Dark2", direction = 1, labels = c("Strong Upwelling", "Weak Upwelling")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = c(0.88,0.9),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x=element_blank(), y='G:SOM')
r 

GonadIMdl <- lmer(GSOM ~ POP * WetWtI +  (1|POP:Site) , data = datUrchinInitial)
summary(GonadIMdl)
anova(GonadIMdl)

#Now move on to experimental results
datUrchinHealthy$Gonad <- datUrchinHealthy$Ganimal - datUrchinHealthy$G_AFDW
datUrchinHealthy$SOM <- datUrchinHealthy$Sanimal - datUrchinHealthy$S_AFDW
datUrchinHealthy$GSOM <- datUrchinHealthy$Gonad/datUrchinHealthy$SOM
datUrchinHealthy$Gonad <- ifelse(datUrchinHealthy$Gonad < 0, 0, datUrchinHealthy$Gonad)

datUrchinHealthy = datUrchinHealthy %>%
  filter(GSOM >= 0) 

datUrchinHealthy_C <- subset(datUrchinHealthy, datUrchinHealthy$CF == 'C')

#Full model
GonadMdl <- lmer(GSOM ~ POP * Trt * WetWtF + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_C)
summary(GonadMdl)
anova(GonadMdl)

#Reduced model
GonadMdl <- lmer(GSOM ~ POP * Trt + POP * WetWtF + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_C)
summary(GonadMdl)
anova(GonadMdl)

GMdlF <- emmeans(GonadMdl, ~ WetWtF)

GonadC <- data.frame(summary(GMdlF)$emmean)
colnames(GonadC) <- 'Growth'
GonadCSE <- data.frame(summary(GMdlF)$SE)
colnames(GonadCSE) <- 'GrowthSE'

GonadC <- cbind(GonadC, GonadCSE)
GonadC$POP <- c('N','S')

Sum_count <- datUrchinHealthy_C %>%
  count(POP, sort = TRUE) 

#Plots of slope and intercepts
s = ggplot(datUrchinHealthy_C, aes(x = WetWtF, y = GSOM), color = POP, fill = POP)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(size = 5, alpha = 1/5, aes(color = POP)) +
  guides(fill = guide_legend(title = 'Population'), colour=guide_legend(title="Population")) +
  scale_color_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  geom_smooth(method = "lm", aes(fill = POP, color = POP)) +
  scale_fill_brewer(palette = "Dark2", direction = 1, labels = c("Strong Upwelling", "Weak Upwelling")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = c(0.88,0.9),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x=element_blank(), y='G:SOM')
s

datUrchinHealthy_CF <- subset(datUrchinHealthy, datUrchinHealthy$POPTrt == 'NORTH NC' | datUrchinHealthy$POPTrt == 'NORTH NF' |
                                datUrchinHealthy$POPTrt == 'SOUTH SC' | datUrchinHealthy$POPTrt == 'SOUTH SF')

#Full model
GonadMdlCF <- lmer(GSOM ~ POP * CF * WetWtF + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_CF)
summary(GonadMdlCF)
anova(GonadMdlCF)

#Reduced model
GonadMdlCF <- lmer(GSOM ~ CF * POP +  POP * WetWtF + (1|Trt:Header) + (1|POP:Site) + (1|Trt:Header:Bin), data = datUrchinHealthy_CF)
summary(GonadMdlCF)
anova(GonadMdlCF)

Sum_count <- datUrchinHealthy_CF %>%
  count(POP, sort = TRUE) 

#Plots of slope and intercepts
t = ggplot(datUrchinHealthy_CF, aes(x = WetWtF, y = GSOM), color = POP, fill = POP)+
  theme_classic() + theme(axis.text = element_text(size = 12)) +
  geom_point(size = 5, alpha = 1/5, aes(color = POP)) +
  guides(fill = guide_legend(title = 'Population'), colour=guide_legend(title="Population")) +
  scale_color_brewer(palette = "Dark2", direction = 1,labels = c("Strong Upwelling", "Weak Upwelling")) +
  geom_smooth(method = "lm", aes(fill = POP, color = POP)) +
  scale_fill_brewer(palette = "Dark2", direction = 1, labels = c("Strong Upwelling", "Weak Upwelling")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5),legend.position = c(0.88,0.9),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 9))+
  labs(x='Wet Weight (g)', y='G:SOM')
t

LA_Supp <- ggarrange(p,o,iG,m,k, labels = c("(A)","(B)","(C)","(D)","(E)"),  align = "hv", ncol = 2, nrow = 3, common.legend = TRUE,
                     legend = "top")
LA_Supp

ggsave(plot = LA_Supp, file = "Fig3_rev.png", 
       type = "cairo-png",  bg = "white",
       width = 35, height = 40, units = "cm", dpi = 600)

CF_Supp <- ggarrange(q,j,n,l, labels = c("(A)","(B)","(C)","(D)"),  align = "hv", ncol = 2, nrow = 2, common.legend = TRUE)
CF_Supp

ggsave(plot = CF_Supp, file = "Fig4_rev.png", 
       type = "cairo-png",  bg = "white",
       width = 35, height = 30, units = "cm", dpi = 600)

GSOM_Fig <- ggarrange(r,s,t, labels = c("(a)","(b)","(c)"),  align = "hv", ncol = 1, nrow = 3, common.legend = TRUE)
GSOM_Fig

ggsave(plot = GSOM_Fig, file = "Fig_S2.png", 
       type = "cairo-png",  bg = "white",
       width = 18, height = 40, units = "cm", dpi = 600)

MortFig <- ggarrange(p,q, labels = c("(A)","(B)"),  align = "hv", ncol = 2, common.legend = FALSE)
MortFig
ggsave(plot = MortFig, file = "Fig_Mort.png", 
       type = "cairo-png",  bg = "white",
       width = 25, height = 12, units = "cm", dpi = 300)



