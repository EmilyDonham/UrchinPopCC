#Written by Emily Donham 11/1/2020
#Plots processed pH data

######################################################################################################
######################################################################################################
library(plyr); library(dplyr);library(broom);library(ggplot2); library(lubridate);library(LoLinR); library(stringr)
library(lme4); library(lmerTest); library(multcomp); library(phytotools); library(googledrive); library(Rmisc)
library(tibble); library(ggpubr); library(wesanderson); library(tidyverse);library(vegan);
library(lsmeans); library(RLRsim); library(pracma); library(chron); library(caTools); library(TTR);
library(zoo); library(datetimeutils); library(metR); library(cowplot)
######################################################################################################
######################################################################################################

rm(list = ls())
######################################################################################################
pHDat = read.csv(paste('data/PA_all.csv'),
                 skip = 7)
pHDat$Lat = as.numeric(pHDat$Lat)
pHDat$Lon = as.numeric(pHDat$Lon)
pHDat$Depth = as.numeric(pHDat$Depth)
pHDat$pH = as.numeric(pHDat$pH)
pHDat$pH.Temp = as.numeric(pHDat$pH.Temp)
pHDat$QF_pH = as.numeric(pHDat$QF_pH )
pHDat$QF_TC = as.numeric(pHDat$QF_TC)
pHDat$DO..umol.kg. = as.numeric(pHDat$DO..umol.kg.)
pHDat$DO.mgL <- pHDat$DO..umol.kg./31.2512
pHDat$DO..Sat = as.numeric(pHDat$DO..Sat)
pHDat$DO.Temp = as.numeric(pHDat$DO.Temp)
pHDat$QF_DO = as.numeric(pHDat$QF_DO)
pHDat$year = as.numeric(str_sub(pHDat$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat$month = as.numeric(str_sub(pHDat$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat$day = as.numeric(str_sub(pHDat$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat$DT = as.POSIXct(paste(pHDat$day, pHDat$month, pHDat$year, 
                            pHDat$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat <- pHDat[complete.cases(pHDat$Lat),]
for (i in 1:nrow(pHDat)) {
  if(pHDat$DT[i] < '2018-05-03 6:00:00' & pHDat$DT[i] > '2017-11-02 22:00') {
    pHDat$DEP[i] = 1
  }else if(pHDat$DT[i] > '2018-05-03 22:00:00' & pHDat$DT[i] < '2018-7-10 6:00') {
    pHDat$DEP[i] = 2
  }else if(pHDat$DT[i] > '2018-07-10 22:00:00' & pHDat$DT[i] < '2018-9-26 6:00') {
    pHDat$DEP[i] = 3
  }else if(pHDat$DT[i] > '2018-09-27 22:00:00' & pHDat$DT[i] < '2018-12-20 6:00') {
    pHDat$DEP[i] = 4
  }else if(pHDat$DT[i] > '2019-04-01 22:00:00' & pHDat$DT[i] < '2019-11-07 6:00') {
    pHDat$DEP[i] = 5
  }else if(pHDat$DT[i] > '2019-11-07 22:00:00' & pHDat$DT[i] < '2020-5-5 6:00') {
    pHDat$DEP[i] = 6
  }else if(pHDat$DT[i] > '2020-07-16 22:00:00' & pHDat$DT[i] < '2020-10-07 6:00') {
    pHDat$DEP[i] = 7
  }else if(pHDat$DT[i] > '2020-10-07 7:30:00' & pHDat$DT[i] < '2021-07-27 8:00') {
    pHDat$DEP[i] = 8
  }else if(pHDat$DT[i] > '2021-07-27 22:00:00' & pHDat$DT[i] < '2021-11-18 6:00') {
    pHDat$DEP[i] = 9
  }else{pHDat$DEP[i] = NA
  }
}

#Average every hour
pH = pHDat %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] < '2018-05-03 6:00:00' & pH$DT[i] > '2017-11-02 22:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-05-03 22:00:00' & pH$DT[i] < '2018-7-10 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-07-10 22:00:00' & pH$DT[i] < '2018-9-26 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-09-27 22:00:00' & pH$DT[i] < '2018-12-20 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-04-01 22:00:00' & pH$DT[i] < '2019-11-07 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-11-07 22:00:00' & pH$DT[i] < '2020-5-5 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2020-07-16 22:00:00' & pH$DT[i] < '2020-10-07 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-10-07 7:30:00' & pH$DT[i] < '2021-07-27 8:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2021-07-27 22:00:00' & pH$DT[i] < '2021-11-18 6:00') {
    pH$DEP[i] = 9
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]


DO <- pHDat[complete.cases(pHDat$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

DO <- DO[complete.cases(DO$DEP),]

for (i in 1:nrow(pHDat)) {
  if(pHDat$QF_TC[i]==1) {
    pHDat$TC_Both[i] = pHDat$pH.Temp[i]
  }else{
    pHDat$TC_Both[i] = pHDat$DO.Temp[i]
  }
}

TC <- pHDat[complete.cases(pHDat$TC_Both),]

TC <- TC %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(TC_Both))
TC <- merge(x = TC, y = pHDat, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] < '2018-05-03 6:00:00') {
    TC$DEP[i] = 1
  }else{
    TC$DEP[i] = TC$DEP[i]
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllPA <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllPA <- merge(x = AllPA, y = TC, by = "DT", all = TRUE)
names(AllPA) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')


AllPA_Good <- AllPA
AllPA <- AllPA[complete.cases(AllPA$pH),]
AllPA <- AllPA[complete.cases(AllPA$DO),]
AllPA <- AllPA[complete.cases(AllPA$TC),]

#Turn all dates to same year in order to separate by month
AllPA <- AllPA %>%
  mutate(date=ymd_hm(format(AllPA$DT, "2017-%m-%d-%H:%M")))

pal = wes_palette("Zissou1",6,type = "continuous")


pHPA <- ggplot(AllPA_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
  pHPA 

TCPA <- ggplot(AllPA_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCPA 

DOPA <- ggplot(AllPA_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOPA 

PA <- ggarrange(pHPA, TCPA, DOPA, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
PA

#Regressions of time series 
pHDO <- lm(AllPA$DO ~ AllPA$pH)
summary(pHDO)
pHTC <- lm(AllPA$TC ~ AllPA$pH)
summary(pHTC)
DOTC <- lm(AllPA$TC ~ AllPA$DO)
summary(DOTC)

#Average values of time series
pHmu <- mean(AllPA$pH); pHmu;
pHstd <- std(AllPA$pH); pHstd;
TCmu <- mean(AllPA$TC); TCmu;
TCstd <- std(AllPA$TC); TCstd;
DOmu <- mean(AllPA$DO); DOmu;
DOstd <- std(AllPA$DO); DOstd;


PAscat <- ggplot(AllPA, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 1.5, position = "jitter") +
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y="Dissolved Oxygen (mg/L)", x="pH") +
  theme_classic() +
  xlim(7.3, 8.3) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
PAscat <- PAscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
PAscat


######################################################################################################
pHDat2 = read.csv(paste('data/VD_all.csv'),
                 skip = 7)
pHDat2$Lat = as.numeric(pHDat2$Lat)
pHDat2$Lon = as.numeric(pHDat2$Lon)
pHDat2$Depth = as.numeric(pHDat2$Depth)
pHDat2$pH = as.numeric(pHDat2$pH)
pHDat2$pH.Temp = as.numeric(pHDat2$pH.Temp)
pHDat2$QF_pH = as.numeric(pHDat2$QF_pH )
pHDat2$QF_TC = as.numeric(pHDat2$QF_TC)
pHDat2$DO..umol.kg. = as.numeric(pHDat2$DO..umol.kg.)
pHDat2$DO.mgL <- pHDat2$DO..umol.kg./31.2512
pHDat2$DO..Sat = as.numeric(pHDat2$DO..Sat)
pHDat2$DO.Temp = as.numeric(pHDat2$DO.Temp)
pHDat2$QF_DO = as.numeric(pHDat2$QF_DO)
pHDat2$year = as.numeric(str_sub(pHDat2$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat2$month = as.numeric(str_sub(pHDat2$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat2$day = as.numeric(str_sub(pHDat2$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat2$DT = as.POSIXct(paste(pHDat2$day, pHDat2$month, pHDat2$year, 
                            pHDat2$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat2 <- pHDat2[complete.cases(pHDat2$Lat),]

for (i in 1:nrow(pHDat2)) {
  if(pHDat2$DT[i] > '2017-11-03 22:00:00' & pHDat2$DT[i] < '2018-04-26 8:00') {
    pHDat2$DEP[i] = 1
  }else if(pHDat2$DT[i] > '2018-04-26 22:00:00' & pHDat2$DT[i] < '2018-07-08 8:00') {
    pHDat2$DEP[i] = 2
  }else if(pHDat2$DT[i] > '2018-07-08 22:00:00' & pHDat2$DT[i] < '2019-01-25 8:00') {
    pHDat2$DEP[i] = 3
  }else if(pHDat2$DT[i] > '2019-01-30 22:00:00' & pHDat2$DT[i] < '2019-04-04 8:00') {
    pHDat2$DEP[i] = 4
  }else if(pHDat2$DT[i] > '2019-04-04 22:00:00' & pHDat2$DT[i] < '2019-11-07 8:00') {
    pHDat2$DEP[i] = 5
  }else if(pHDat2$DT[i] > '2019-11-07 22:00:00' & pHDat2$DT[i] < '2020-05-05 8:00') {
    pHDat2$DEP[i] = 6
  }else if(pHDat2$DT[i] > '2020-07-16 22:00:00' & pHDat2$DT[i] < '2020-10-19 8:00') {
    pHDat2$DEP[i] = 7
  }else if(pHDat2$DT[i] > '2020-10-19 22:00:00' & pHDat2$DT[i] < '2021-07-27 8:00') {
    pHDat2$DEP[i] = 8
  }else{pHDat2$DEP[i] = NA
  }
}

#Average every day
pH <- subset(pHDat2, pHDat2$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat2, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-11-03 22:00:00' & pH$DT[i] < '2018-04-26 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-04-26 22:00:00' & pH$DT[i] < '2018-07-08 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-07-08 22:00:00' & pH$DT[i] < '2019-01-25 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2019-01-30 22:00:00' & pH$DT[i] < '2019-04-04 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-04-04 22:00:00' & pH$DT[i] < '2019-11-07 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-11-07 22:00:00' & pH$DT[i] < '2020-05-05 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2020-07-16 22:00:00' & pH$DT[i] < '2020-10-19 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-10-19 22:00:00' & pH$DT[i] < '2021-07-27 6:00') {
    pH$DEP[i] = 8
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat2[complete.cases(pHDat2$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat2, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-11-03 22:00:00' & DO$DT[i] < '2018-04-26 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-04-26 22:00:00' & DO$DT[i] < '2018-07-08 6:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-07-08 22:00:00' & DO$DT[i] < '2019-01-25 6:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2019-01-30 22:00:00' & DO$DT[i] < '2019-04-04 6:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-04-04 22:00:00' & DO$DT[i] < '2019-11-07 6:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-11-07 22:00:00' & DO$DT[i] < '2020-05-05 6:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2020-07-16 22:00:00' & DO$DT[i] < '2020-10-19 6:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-10-19 22:00:00' & DO$DT[i] < '2021-07-27 6:00') {
    DO$DEP[i] = 8
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]

TC <- pHDat2[complete.cases(pHDat2$pH.Temp),]

TC = TC %>%
  filter(QF_TC == "1") 

TC <- TC %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH.Temp))
TC <- TC[complete.cases(TC$avg),]
TC <- merge(x = TC, y = pHDat2, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-11-03 22:00:00' & TC$DT[i] < '2018-04-26 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-04-26 22:00:00' & TC$DT[i] < '2018-07-08 6:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-07-08 22:00:00' & TC$DT[i] < '2019-01-25 6:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2019-01-30 22:00:00' & TC$DT[i] < '2019-04-04 6:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-04-04 22:00:00' & TC$DT[i] < '2019-11-07 6:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-11-07 22:00:00' & TC$DT[i] < '2020-05-05 6:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2020-07-16 22:00:00' & TC$DT[i] < '2020-10-19 6:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-10-19 22:00:00' & TC$DT[i] < '2021-07-27 6:00') {
    TC$DEP[i] = 8 
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllVD <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllVD <- merge(x = AllVD, y = TC, by = "DT", all = TRUE)

names(AllVD) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllVD_Good <- AllVD

AllVD <- AllVD[complete.cases(AllVD$pH),]
AllVD <- AllVD[complete.cases(AllVD$DO),]
AllVD <- AllVD[complete.cases(AllVD$TC),]

#Turn all dates to same year in order to separate by month
AllVD <- AllVD %>%
  mutate(date=ymd_hm(format(AllVD$DT, "2017-%m-%d-%H:%M")))

pHVD <- ggplot(AllVD_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHVD

TCVD <- ggplot(AllVD_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCVD 

DOVD <- ggplot(AllVD_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOVD 

VD <- ggarrange(pHVD, TCVD, DOVD, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
VD

#Regressions of time series 
pHDO <- lm(AllVD$DO ~ AllVD$pH)
summary(pHDO)
pHTC <- lm(AllVD$TC ~ AllVD$pH)
summary(pHTC)
DOTC <- lm(AllVD$TC ~ AllVD$DO)
summary(DOTC)

pHmu <- mean(AllVD$pH); pHmu;
pHstd <- std(AllVD$pH); pHstd;
TCmu <- mean(AllVD$TC); TCmu;
TCstd <- std(AllVD$TC); TCstd;
DOmu <- mean(AllVD$DO); DOmu;
DOstd <- std(AllVD$DO); DOstd;

VDscat <- ggplot(AllVD, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 1.5, position = "jitter") +
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y="Dissolved Oxygen (mg/L)", x="pH") +
  theme_classic() +
  xlim(7.3, 8.3) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
VDscat <- VDscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
VDscat


######################################################################################################
pHDat3 = read.csv(paste('data/LB_all.csv'),
                 skip = 7)
pHDat3$Lat = as.numeric(pHDat3$Lat)
pHDat3$Lon = as.numeric(pHDat3$Lon)
pHDat3$Depth = as.numeric(pHDat3$Depth)
pHDat3$pH = as.numeric(pHDat3$pH)
pHDat3$pH.Temp = as.numeric(pHDat3$pH.Temp)
pHDat3$QF_pH = as.numeric(pHDat3$QF_pH )
pHDat3$QF_TC = as.numeric(pHDat3$QF_TC)
pHDat3$DO..umol.kg. = as.numeric(pHDat3$DO..umol.kg.)
pHDat3$DO.mgL <- pHDat3$DO..umol.kg./31.2512
pHDat3$DO..Sat = as.numeric(pHDat3$DO..Sat)
pHDat3$DO.Temp = as.numeric(pHDat3$DO.Temp)
pHDat3$QF_DO = as.numeric(pHDat3$QF_DO)
pHDat3$year = as.numeric(str_sub(pHDat3$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat3$month = as.numeric(str_sub(pHDat3$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat3$day = as.numeric(str_sub(pHDat3$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat3$DT = as.POSIXct(paste(pHDat3$day, pHDat3$month, pHDat3$year, 
                             pHDat3$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat3 <- pHDat3[complete.cases(pHDat3$Lat),]

for (i in 1:nrow(pHDat3)) {
  if(pHDat3$DT[i] > '2017-12-12 24:00:00' & pHDat3$DT[i] < '2018-01-09 6:00') {
    pHDat3$DEP[i] = 1
  }else if(pHDat3$DT[i] > '2018-06-12 24:00:00' & pHDat3$DT[i] < '2018-09-05 6:00') {
    pHDat3$DEP[i] = 2
  }else if(pHDat3$DT[i] > '2018-09-05 24:00:00' & pHDat3$DT[i] < '2018-12-10 6:00') {
    pHDat3$DEP[i] = 3
  }else if(pHDat3$DT[i] > '2018-12-10 24:00:00' & pHDat3$DT[i] < '2019-03-14 6:00') {
    pHDat3$DEP[i] = 4
  }else if(pHDat3$DT[i] > '2019-03-14 24:00:00' & pHDat3$DT[i] < '2019-06-10 6:00') {
    pHDat3$DEP[i] = 5
  }else if(pHDat3$DT[i] > '2019-06-10 24:00:00' & pHDat3$DT[i] < '2019-12-18 6:00') {
    pHDat3$DEP[i] = 6
  }else if(pHDat3$DT[i] > '2019-12-18 24:00:00' & pHDat3$DT[i] < '2020-03-04 6:00') {
    pHDat3$DEP[i] = 7
  }else if(pHDat3$DT[i] > '2020-03-04 24:00:00' & pHDat3$DT[i] < '2020-5-05 6:00') {
    pHDat3$DEP[i] = 8
  }else if(pHDat3$DT[i] > '2020-11-19 24:00:00' & pHDat3$DT[i] < '2021-04-02 6:00') {
    pHDat3$DEP[i] = 9
  }else if(pHDat3$DT[i] > '2021-04-02 24:00:00' & pHDat3$DT[i] < '2021-10-29 6:00') {
    pHDat3$DEP[i] = 10
  }else{pHDat3$DEP[i] = NA
  }
}


#Average every hour
pH <- subset(pHDat3, pHDat3$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat3, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-12-12 24:00:00' & pH$DT[i] < '2018-01-09 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-06-12 24:00:00' & pH$DT[i] < '2018-09-05 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-09-05 24:00:00' & pH$DT[i] < '2018-12-10 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-12-10 24:00:00' & pH$DT[i] < '2019-03-14 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-03-14 24:00:00' & pH$DT[i] < '2019-06-10 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-06-10 24:00:00' & pH$DT[i] < '2019-12-18 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2019-12-18 24:00:00' & pH$DT[i] < '2020-03-04 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-03-04 24:00:00' & pH$DT[i] < '2020-5-01 6:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2020-11-19 24:00:00' & pH$DT[i] < '2021-03-02 6:00') {
    pH$DEP[i] = 9
  }else if(pH$DT[i] > '2021-04-02 24:00:00' & pH$DT[i] < '2021-10-29 6:00') {
    pH$DEP[i] = 10
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat3[complete.cases(pHDat3$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat3, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-12-12 24:00:00' & DO$DT[i] < '2018-01-09 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-06-12 24:00:00' & DO$DT[i] < '2018-09-05 6:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-09-05 24:00:00' & DO$DT[i] < '2018-12-10 6:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2018-12-10 24:00:00' & DO$DT[i] < '2019-03-01 6:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-03-14 24:00:00' & DO$DT[i] < '2019-06-10 6:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-06-10 24:00:00' & DO$DT[i] < '2019-12-01 6:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2019-12-18 24:00:00' & DO$DT[i] < '2020-03-04 6:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-03-04 24:00:00' & DO$DT[i] < '2020-5-01 6:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2020-11-19 24:00:00' & DO$DT[i] < '2021-03-02 6:00') {
    DO$DEP[i] = 9
  }else if(DO$DT[i] > '2021-04-02 24:00:00' & DO$DT[i] < '2021-10-29 6:00') {
    DO$DEP[i] = 10
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]

TC <- pHDat3[complete.cases(pHDat3$pH.Temp),]

TC = TC %>%
  filter(QF_TC == "1") 

TC <- TC %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH.Temp))
TC <- TC[complete.cases(TC$avg),]
TC <- merge(x = TC, y = pHDat3, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-12-12 24:00:00' & TC$DT[i] < '2018-01-09 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-06-12 24:00:00' & TC$DT[i] < '2018-09-05 6:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-09-05 24:00:00' & TC$DT[i] < '2018-12-10 6:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2018-12-10 24:00:00' & TC$DT[i] < '2019-03-01 6:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-03-14 24:00:00' & TC$DT[i] < '2019-06-10 6:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-06-10 24:00:00' & TC$DT[i] < '2019-12-01 6:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2019-12-18 24:00:00' & TC$DT[i] < '2020-03-04 6:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-03-04 24:00:00' & TC$DT[i] < '2020-5-01 6:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2020-11-19 24:00:00' & TC$DT[i] < '2021-03-02 6:00') {
    TC$DEP[i] = 9
  }else if(TC$DT[i] > '2021-04-02 24:00:00' & TC$DT[i] < '2021-10-29 6:00') {
    TC$DEP[i] = 10 
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllLB <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllLB <- merge(x = AllLB, y = TC, by = "DT", all = TRUE)

names(AllLB) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllLB_Good <- AllLB

AllLB <- AllLB[complete.cases(AllLB$pH),]
AllLB <- AllLB[complete.cases(AllLB$DO),]
AllLB <- AllLB[complete.cases(AllLB$TC),]

#Turn all dates to same year in order to separate by month
AllLB <- AllLB %>%
  mutate(date=ymd_hm(format(AllLB$DT, "2017-%m-%d-%H:%M")))

pHLB <- ggplot(AllLB_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHLB

TCLB <- ggplot(AllLB_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCLB 

DOLB <- ggplot(AllLB_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOLB 

LB <- ggarrange(pHLB, TCLB, DOLB, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
LB

#Regressions of time series 
pHDO <- lm(AllLB$DO ~ AllLB$pH)
summary(pHDO)
pHTC <- lm(AllLB$TC ~ AllLB$pH)
summary(pHTC)
DOTC <- lm(AllLB$TC ~ AllLB$DO)
summary(DOTC)

pHmu <- mean(AllLB$pH); pHmu;
pHstd <- std(AllLB$pH); pHstd;
TCmu <- mean(AllLB$TC); TCmu;
TCstd <- std(AllLB$TC); TCstd;
DOmu <- mean(AllLB$DO); DOmu;
DOstd <- std(AllLB$DO); DOstd;

LBscat <- ggplot(AllLB, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 1.5, position = "jitter") +
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y="Dissolved Oxygen (mg/L)", x="pH") +
  theme_classic() +
  xlim(7.3, 8.3) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
LBscat <- LBscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
LBscat

######################################################################################################
pHDat4 = read.csv(paste('data/CI_all.csv'),
                  skip = 7)
pHDat4$Lat = as.numeric(pHDat4$Lat)
pHDat4$Lon = as.numeric(pHDat4$Lon)
pHDat4$Depth = as.numeric(pHDat4$Depth)
pHDat4$pH = as.numeric(pHDat4$pH)
pHDat4$pH.Temp = as.numeric(pHDat4$pH.Temp)
pHDat4$QF_pH = as.numeric(pHDat4$QF_pH )
pHDat4$QF_TC = as.numeric(pHDat4$QF_TC)
pHDat4$DO..umol.kg. = as.numeric(pHDat4$DO..umol.kg.)
pHDat4$DO.mgL <- pHDat4$DO..umol.kg./31.2512
pHDat4$DO..Sat = as.numeric(pHDat4$DO..Sat)
pHDat4$DO.Temp = as.numeric(pHDat4$DO.Temp)
pHDat4$QF_DO = as.numeric(pHDat4$QF_DO)
pHDat4$year = as.numeric(str_sub(pHDat4$pH_Date_UTC_ddmmyyyy,-4,-1))
pHDat4$month = as.numeric(str_sub(pHDat4$pH_Date_UTC_ddmmyyyy,-6,-5))
pHDat4$day = as.numeric(str_sub(pHDat4$pH_Date_UTC_ddmmyyyy,-8,-7))
pHDat4$DT = as.POSIXct(paste(pHDat4$day, pHDat4$month, pHDat4$year, 
                             pHDat4$pH_Time_UTC_hhmmss), format = "%d%m%Y %H:%M:%S", tz = "UTC")

pHDat4 <- pHDat4[complete.cases(pHDat4$Lat),]


for (i in 1:nrow(pHDat4)) {
  if(pHDat4$DT[i] > '2017-12-11 22:00:00' & pHDat4$DT[i] < '2018-03-01 6:00') {
    pHDat4$DEP[i] = 1
  }else if(pHDat4$DT[i] > '2018-03-09 22:00:00' & pHDat4$DT[i] < '2018-06-11 6:00') {
    pHDat4$DEP[i] = 2
  }else if(pHDat4$DT[i] > '2018-06-11 22:00:00' & pHDat4$DT[i] < '2018-09-05 6:00') {
    pHDat4$DEP[i] = 3
  }else if(pHDat4$DT[i] > '2018-12-10 22:00:00' & pHDat4$DT[i] < '2019-03-14 6:00') {
    pHDat4$DEP[i] = 4
  }else if(pHDat4$DT[i] > '2019-03-14 22:00:00' & pHDat4$DT[i] < '2019-06-10 6:00') {
    pHDat4$DEP[i] = 5
  }else if(pHDat4$DT[i] > '2019-06-10 22:00:00' & pHDat4$DT[i] < '2019-12-18 6:00') {
    pHDat4$DEP[i] = 6
  }else if(pHDat4$DT[i] > '2019-12-18 22:00:00' & pHDat4$DT[i] < '2020-03-03 6:00') {
    pHDat4$DEP[i] = 7
  }else if(pHDat4$DT[i] > '2020-03-03 22:00:00' & pHDat4$DT[i] < '2020-10-14 6:00') {
    pHDat4$DEP[i] = 8
  }else if(pHDat4$DT[i] > '2020-10-14 22:00:00' & pHDat4$DT[i] < '2021-06-26 6:00') {
    pHDat4$DEP[i] = 9
  }else{pHDat4$DEP[i] = NA
  }
}

#Average every hour
pH <- subset(pHDat4, pHDat4$pH< 8.5)
pH = pH %>%
  filter(QF_pH == "1") 

pH <- pH %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(pH))
pH <- pH[complete.cases(pH$avg),]
pH <- merge(x = pH, y = pHDat4, by = "DT", all.x = TRUE)
pH <- pH[c(1,2,20)]

for (i in 1:nrow(pH)) {
  if(pH$DT[i] > '2017-12-11 22:00:00' & pH$DT[i] < '2018-03-01 6:00') {
    pH$DEP[i] = 1
  }else if(pH$DT[i] > '2018-03-09 22:00:00' & pH$DT[i] < '2018-06-11 6:00') {
    pH$DEP[i] = 2
  }else if(pH$DT[i] > '2018-06-11 22:00:00' & pH$DT[i] < '2018-09-05 6:00') {
    pH$DEP[i] = 3
  }else if(pH$DT[i] > '2018-12-10 22:00:00' & pH$DT[i] < '2019-03-14 6:00') {
    pH$DEP[i] = 4
  }else if(pH$DT[i] > '2019-03-14 22:00:00' & pH$DT[i] < '2019-06-10 6:00') {
    pH$DEP[i] = 5
  }else if(pH$DT[i] > '2019-06-10 22:00:00' & pH$DT[i] < '2019-12-18 6:00') {
    pH$DEP[i] = 6
  }else if(pH$DT[i] > '2019-12-18 22:00:00' & pH$DT[i] < '2020-03-03 6:00') {
    pH$DEP[i] = 7
  }else if(pH$DT[i] > '2020-03-03 22:00:00' & pH$DT[i] < '2020-10-14 6:00') {
    pH$DEP[i] = 8
  }else if(pH$DT[i] > '2020-10-14 22:00:00' & pH$DT[i] < '2021-06-26 6:00') {
    pH$DEP[i] = 9
  }else{pH$DEP[i] = NA
  }
}
pH <- pH[complete.cases(pH$DEP),]

DO <- pHDat4[complete.cases(pHDat4$DO.mgL),]
DO = DO %>%
  filter(QF_DO == "1")
DO <- DO %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(DO.mgL))
DO <- DO[complete.cases(DO$avg),]
DO <- merge(x = DO, y = pHDat4, by = "DT", all.x = TRUE)
DO <- DO[c(1,2,20)]

for (i in 1:nrow(DO)) {
  if(DO$DT[i] > '2017-12-11 22:00:00' & DO$DT[i] < '2018-03-01 6:00') {
    DO$DEP[i] = 1
  }else if(DO$DT[i] > '2018-03-09 22:00:00' & DO$DT[i] < '2018-06-11 6:00') {
    DO$DEP[i] = 2
  }else if(DO$DT[i] > '2018-06-11 22:00:00' & DO$DT[i] < '2018-09-05 6:00') {
    DO$DEP[i] = 3
  }else if(DO$DT[i] > '2018-12-10 22:00:00' & DO$DT[i] < '2019-03-14 6:00') {
    DO$DEP[i] = 4
  }else if(DO$DT[i] > '2019-03-14 22:00:00' & DO$DT[i] < '2019-06-10 6:00') {
    DO$DEP[i] = 5
  }else if(DO$DT[i] > '2019-06-10 22:00:00' & DO$DT[i] < '2019-12-18 6:00') {
    DO$DEP[i] = 6
  }else if(DO$DT[i] > '2019-12-18 22:00:00' & DO$DT[i] < '2020-03-03 6:00') {
    DO$DEP[i] = 7
  }else if(DO$DT[i] > '2020-03-03 22:00:00' & DO$DT[i] < '2020-10-14 6:00') {
    DO$DEP[i] = 8
  }else if(DO$DT[i] > '2020-10-14 22:00:00' & DO$DT[i] < '2021-06-26 6:00') {
    DO$DEP[i] = 9
  }else{DO$DEP[i] = NA
  }
}
DO <- DO[complete.cases(DO$DEP),]


for (i in 1:nrow(pHDat4)) {
  if(pHDat4$QF_TC[i]==1) {
    pHDat4$TC_Both[i] = pHDat4$pH.Temp[i]
  }else{
    pHDat4$TC_Both[i] = pHDat4$DO.Temp[i]
  }
}

TC <- pHDat4[complete.cases(pHDat4$TC_Both),]

TC <- TC %>%
  mutate(DT = floor_date(DT,"1 days")) %>%
  group_by(DT) %>%
  summarize(avg = mean(TC_Both))
TC <- merge(x = TC, y = pHDat4, by = "DT", all.x = TRUE)
TC <- TC[c(1,2,20)]

for (i in 1:nrow(TC)) {
  if(TC$DT[i] > '2017-12-11 22:00:00' & TC$DT[i] < '2018-03-01 6:00') {
    TC$DEP[i] = 1
  }else if(TC$DT[i] > '2018-03-09 22:00:00' & TC$DT[i] < '2018-06-11 6:00') {
    TC$DEP[i] = 2
  }else if(TC$DT[i] > '2018-06-11 22:00:00' & TC$DT[i] < '2018-09-05 6:00') {
    TC$DEP[i] = 3
  }else if(TC$DT[i] > '2018-12-10 22:00:00' & TC$DT[i] < '2019-03-14 6:00') {
    TC$DEP[i] = 4
  }else if(TC$DT[i] > '2019-03-14 22:00:00' & TC$DT[i] < '2019-06-10 6:00') {
    TC$DEP[i] = 5
  }else if(TC$DT[i] > '2019-06-10 22:00:00' & TC$DT[i] < '2019-12-18 6:00') {
    TC$DEP[i] = 6
  }else if(TC$DT[i] > '2019-12-18 22:00:00' & TC$DT[i] < '2020-03-03 6:00') {
    TC$DEP[i] = 7
  }else if(TC$DT[i] > '2020-03-03 22:00:00' & TC$DT[i] < '2020-10-14 6:00') {
    TC$DEP[i] = 8
  }else if(TC$DT[i] > '2020-10-14 22:00:00' & TC$DT[i] < '2021-06-26 6:00') {
    TC$DEP[i] = 9
  }else{TC$DEP[i] = NA
  }
}
TC <- TC[complete.cases(TC$DEP),]

AllCI <- merge(x = pH, y = DO, by = "DT", all = TRUE)
AllCI <- merge(x = AllCI, y = TC, by = "DT", all = TRUE)

names(AllCI) <- c('DT', 'pH','Dep1','DO','Dep2','TC','Dep3')
AllCI_Good <- AllCI


AllCI <- AllCI[complete.cases(AllCI$pH),]
AllCI <- AllCI[complete.cases(AllCI$DO),]
AllCI <- AllCI[complete.cases(AllCI$TC),]

#Turn all dates to same year in order to separate by month
AllCI <- AllCI %>%
  mutate(date=ymd_hm(format(AllCI$DT, "2017-%m-%d-%H:%M")))

pHCI <- ggplot(AllCI_Good, aes(x = DT, y = pH, group = Dep1)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)", y = "pH", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
pHCI

TCCI <- ggplot(AllCI_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line() +
  theme_classic() +
  labs(x="Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0))
TCCI 

DOCI <- ggplot(AllCI_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line() +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x="Date (MM/YY)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0))
DOCI 

CI <- ggarrange(pHCI, TCCI, DOCI, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                legend = "bottom", common.legend = FALSE, align = "hv")
CI

#Regressions of time series 
pHDO <- lm(AllCI$DO ~ AllCI$pH)
summary(pHDO)
pHTC <- lm(AllCI$TC ~ AllCI$pH)
summary(pHTC)
DOTC <- lm(AllCI$TC ~ AllCI$DO)
summary(DOTC)

pHmu <- mean(AllCI$pH); pHmu;
pHstd <- std(AllCI$pH); pHstd;
TCmu <- mean(AllCI$TC); TCmu;
TCstd <- std(AllCI$TC); TCstd;
DOmu <- mean(AllCI$DO); DOmu;
DOstd <- std(AllCI$DO); DOstd;

CIscat <- ggplot(AllCI, aes(x = pH, y = DO, color = TC)) + 
  geom_point(size = 1.5, position = "jitter") +
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24)) +
  labs(y="Dissolved Oxygen (mg/L)", x="pH") +
  theme_classic() +
  xlim(7.3, 8.3) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "right", legend.key.width = unit(1,'cm'))
CIscat <- CIscat + guides(colour = guide_colourbar(title.position="top",title.hjust =0.5))+ scale_y_continuous(limits = c(0.5,11), breaks = c(2,4,6,8,10))
CIscat



Scatters <- ggarrange(PAscat, VDscat, LBscat, CIscat, labels = c("(a)", "(b)", "(c)","(d)"), ncol = 2, nrow = 2, 
                legend = "bottom", common.legend = TRUE, align = "hv")
Scatters

ggsave(plot = Scatters, file = "Env_Scatters.png", 
       type = "cairo-png",  bg = "transparent",
       width = 30, height = 30, units = "cm", dpi = 300)


################################################################################################################
DS = read.csv('data/AllDS.csv')
# Convert time
DS$DT <- force_tz(as.POSIXct(with(DS, mdy(DATE)), origin = as.POSIXct("1970-01-01",TZ = "America/Los_Angeles")), tzone = "America/Los_Angeles")

#Choose just samples from bins
DS_Bin <- DS[ which(DS$SH=='S'),]

pHave <- summarySE(data=DS_Bin, measurevar = "pH.out",
                   groupvars = c("HEADER", "DATE"))
pHave$Trt <- c(rep('SC',6),rep('NC',6),rep('SC',6),rep('NC',6),
               rep('NF',6),rep('SC',6),rep('NF',6),rep('SF',6),rep('NC',6),
               rep('SF',6),rep('NF',6),rep('SF',6))

pHave2 <- summarySE(data=pHave, measurevar = "pH.out", groupvars = c("HEADER"))
pHave2$Trt <- c(rep('SC',1),rep('NC',1),rep('SC',1),rep('NC',1),
                rep('NF',1),rep('SC',1),rep('NF',1),rep('SF',1),rep('NC',1),
                rep('SF',1),rep('NF',1),rep('SF',1))

pHave3 <- summarySE(data=pHave2, measurevar = "pH.out", groupvars = c("Trt"))

TCave <- summarySE(data=DS_Bin, measurevar = "Tout",
                   groupvars = c("HEADER", "DATE"))
TCave$Trt <- c(rep('SC',6),rep('NC',6),rep('SC',6),rep('NC',6),
               rep('NF',6),rep('SC',6),rep('NF',6),rep('SF',6),rep('NC',6),
               rep('SF',6),rep('NF',6),rep('SF',6))

TCave2 <- summarySE(data=TCave, measurevar = "Tout", groupvars = c("HEADER"))
TCave2$Trt <- c(rep('SC',1),rep('NC',1),rep('SC',1),rep('NC',1),
                rep('NF',1),rep('SC',1),rep('NF',1),rep('SF',1),rep('NC',1),
                rep('SF',1),rep('NF',1),rep('SF',1))
TCave3 <- summarySE(data=TCave2, measurevar = "Tout", groupvars = c("Trt"))


DOave <- summarySE(data=DS_Bin, measurevar = "DO",
                   groupvars = c("HEADER", "DATE"))
DOave$Trt <- c(rep('SC',6),rep('NC',6),rep('SC',6),rep('NC',6),
               rep('NF',6),rep('SC',6),rep('NF',6),rep('SF',6),rep('NC',6),
               rep('SF',6),rep('NF',6),rep('SF',6))

DOave2 <- summarySE(data=DOave, measurevar = "DO", groupvars = c("HEADER"))
DOave2$Trt <- c(rep('SC',1),rep('NC',1),rep('SC',1),rep('NC',1),
                rep('NF',1),rep('SC',1),rep('NF',1),rep('SF',1),rep('NC',1),
                rep('SF',1),rep('NF',1),rep('SF',1))
DOave3 <- summarySE(data=DOave2, measurevar = "DO", groupvars = c("Trt"))

Discrete <- merge(pHave, TCave, by = c("HEADER", "DATE", "Trt"))
Discrete <- merge(Discrete, DOave, by = c("HEADER", "DATE", "Trt"))
Discrete2 = subset(Discrete, select = c("DATE","HEADER","Trt","pH.out","Tout","DO") )
names(Discrete2) <- c('Date','HEADER','Trt','pH','TC','DO')
Discrete2$Date <- as.POSIXct(with(Discrete2, mdy(Date)), origin = as.POSIXct("1970-01-01",TZ = "America/Los_Angeles"))


PA_day <- AllPA[c(1,2,4,6)]
PA_day$Header <- "NC"
PA_day$Trt <- "PA"

VD_day <- AllVD[c(1,2,4,6)]
VD_day$Header <- "NC"
VD_day$Trt <- "VD"

CI_day <- AllCI[c(1,2,4,6)]
CI_day$Header <- "SC"
CI_day$Trt <- "CI"

LB_day <- AllLB[c(1,2,4,6)]
LB_day$Header <- "SC"
LB_day$Trt <- "LB"

LB_day <- LB_day[c(1,5,6,2,4,3)]
CI_day <- CI_day[c(1,5,6,2,4,3)]
VD_day <- VD_day[c(1,5,6,2,4,3)]
PA_day <- PA_day[c(1,5,6,2,4,3)]
SensData <- rbind(LB_day, CI_day, VD_day, PA_day)

names(SensData) <- c('Date','HEADER','Trt','pH','TC','DO')

EnvData <- rbind(Discrete2,SensData)

EnvDataSub <- subset(EnvData, select = c("pH","TC","DO") )

DisPCA <- merge(pHave3, TCave3, by = c( "Trt"))
DisPCA2 <- merge(DisPCA, DOave3, by = c("Trt"))
DisPCA3 = subset(Discrete, select = c("HEADER","Trt","pH.out","Tout","DO") )
names(DisPCA3) <- c('HEADER','Trt','pH','TC','DO')

 TS.PCA <- prcomp(EnvData[4:6]) #Run PCA of pH, DO, tempDO/temppH
 summary(TS.PCA) #Look at the amount of variance explained
 Load <- data.frame(TS.PCA$x) #Turn loadings into dataframe
 TrtPCA <- data.frame(c(EnvData, Load)) #Combine PCA loadings with OG data
 Ar <- data.frame(TS.PCA$rotation) #Extract vectors
 TrtPCA$Trt <- factor(TrtPCA$Trt, c("PA","VD","NC","NF","CI","LB","SC","SF"))
 for (i in 1:nrow(TrtPCA)) {
   if(TrtPCA$Trt[i] == "PA") {
     TrtPCA$Shape[i] = 'S'
   }else if(TrtPCA$Trt[i] == "VD") {
     TrtPCA$Shape[i] = 'S'
   }else if(TrtPCA$Trt[i] == "NC") {
     TrtPCA$Shape[i] = 'E'
   }else if(TrtPCA$Trt[i] == "NF") {
     TrtPCA$Shape[i] = 'E'
   }else if(TrtPCA$Trt[i] == "CI") {
     TrtPCA$Shape[i] = 'S'
   }else if(TrtPCA$Trt[i] == "LB") {
     TrtPCA$Shape[i] = 'S'
   }else if(TrtPCA$Trt[i] == "SC") {
     TrtPCA$Shape[i] = 'E'
   }else if(TrtPCA$Trt[i] == "SF") {
     TrtPCA$Shape[i] = 'E'
   }
 }
 
 for (i in 1:nrow(TrtPCA)) {
   if(TrtPCA$Trt[i] == "PA") {
     TrtPCA$Color[i] = 'NC'
   }else if(TrtPCA$Trt[i] == "VD") {
     TrtPCA$Color[i] = 'NC'
   }else if(TrtPCA$Trt[i] == "NC") {
     TrtPCA$Color[i] = 'NC'
   }else if(TrtPCA$Trt[i] == "NF") {
     TrtPCA$Color[i] = 'NF'
   }else if(TrtPCA$Trt[i] == "CI") {
     TrtPCA$Color[i] = 'SC'
   }else if(TrtPCA$Trt[i] == "LB") {
     TrtPCA$Color[i] = 'SC'
   }else if(TrtPCA$Trt[i] == "SC") {
     TrtPCA$Color[i] = 'SC'
   }else if(TrtPCA$Trt[i] == "SF") {
     TrtPCA$Color[i] = 'SF'
   }
 }
 
 for (i in 1:nrow(TrtPCA)) {
   if(TrtPCA$Trt[i] == "PA") {
     TrtPCA$CF[i] = 'Fi'
   }else if(TrtPCA$Trt[i] == "VD") {
     TrtPCA$CF[i] = 'Fi'
   }else if(TrtPCA$Trt[i] == "NC") {
     TrtPCA$CF[i] = 'C'
   }else if(TrtPCA$Trt[i] == "NF") {
     TrtPCA$CF[i] = 'F'
   }else if(TrtPCA$Trt[i] == "CI") {
     TrtPCA$CF[i] = 'Fi'
   }else if(TrtPCA$Trt[i] == "LB") {
     TrtPCA$CF[i] = 'Fi'
   }else if(TrtPCA$Trt[i] == "SC") {
     TrtPCA$CF[i] = 'C'
   }else if(TrtPCA$Trt[i] == "SF") {
     TrtPCA$CF[i] = 'F'
   }
 }
 
 for (i in 1:nrow(TrtPCA)) {
   if(TrtPCA$Trt[i] == "PA") {
     TrtPCA$WS[i] = 'S'
   }else if(TrtPCA$Trt[i] == "VD") {
     TrtPCA$WS[i] = 'S'
   }else if(TrtPCA$Trt[i] == "NC") {
     TrtPCA$WS[i] = 'S'
   }else if(TrtPCA$Trt[i] == "NF") {
     TrtPCA$WS[i] = 'S'
   }else if(TrtPCA$Trt[i] == "CI") {
     TrtPCA$WS[i] = 'W'
   }else if(TrtPCA$Trt[i] == "LB") {
     TrtPCA$WS[i] = 'W'
   }else if(TrtPCA$Trt[i] == "SC") {
     TrtPCA$WS[i] = 'W'
   }else if(TrtPCA$Trt[i] == "SF") {
     TrtPCA$WS[i] = 'W'
   }
 }
 
PCA <-ggplot(data=TrtPCA, aes(x=PC1, y=PC2, color = Color)) +
   geom_point(aes(fill = Color, shape = Shape),color = "black", size = 4, alpha = 1/3) + 
   stat_ellipse(aes(PC1, PC2, group = Trt, color = Color), size = 1, type = "norm") +
   theme(legend.key = element_rect(colour = NA, fill = NA),)  +
   scale_shape_manual(values=c(21, 22, 23, 24, 25))+
   scale_fill_brewer(palette = "Spectral", direction=-1, labels = c("North Current","North Future",
                                                                    "South Current","South Future")) +
   scale_color_brewer(palette = "Spectral", direction=-1, labels =  c("North Current","North Future",
                                                                     "South Current","South Future"))
PCA 

#Plotting Average PC scores per group
TrtPCA$Group <- str_c(TrtPCA$CF, TrtPCA$WS)

PC1Ave <- summarySE(data=TrtPCA, measurevar = "PC1", groupvars = c("Group"))
PC2Ave <- summarySE(data=TrtPCA, measurevar = "PC2", groupvars = c("Group"))
PCAAve <- merge(PC1Ave, PC2Ave, by = c( "Group"))
PCAAve$Group <- factor(PCAAve$Group, c("CS","CW","FS","FW","FiS","FiW"))
PCAAve$CF <- c("C","C","Fi","Fi","F","F")
PCAAve$WS <- c("S","W","S","W","S","W")
PCAAve <- subset(PCAAve, PCAAve$CF != 'Fi')

SensorAve <- summarySE(data=SensData, measurevar = c("DO"), groupvars = c("Trt"))



PCA <-ggplot(data=PCAAve, aes(x=PC1, y=PC2, fill = WS, shape = CF)) +
  geom_point(data = TrtPCA, aes(fill = WS, shape = CF), color = "black", size = 4, alpha = 1/4) + 
  geom_errorbar(aes(ymin = PC2-sd.y, ymax = PC2+sd.y),color = "black", width = 0.2, alpha = 2/3)+
  geom_errorbar(aes(xmin = PC1-sd.x, xmax = PC1+sd.x),color = "black", width = 0.2, alpha = 2/3)+
  geom_point(color = "black", size = 6) + 
  annotate("segment", x = 0, xend = Ar$PC1[1]*43, y = 0, yend = Ar$PC2[1]*43, colour = "black", size=0.5, alpha=0.9, arrow=arrow())+
  annotate("segment", x = 0, xend = Ar$PC1[2]*6.05, y = 0, yend = Ar$PC2[2]*6.05, colour = "black", size=0.5, alpha=0.9, arrow=arrow())+
  annotate("segment", x = 0, xend = Ar$PC1[3]*3.5, y = 0, yend = Ar$PC2[3]*3.5, colour = "black", size=0.5, alpha=0.9, arrow=arrow())+
  annotate("text", x = -1.42, y = 3.02, label = "pH")+
  annotate("text", x = -6.01, y = -1.7, label = "Temperature")+
  annotate("text", x = -.95, y = 3.50, label = "Dissolved Oxygen")+
  theme(legend.key = element_rect(colour = NA, fill = NA),)  +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("PC1 (92%)") +
  ylab("PC2 (8%)") +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12)) +
  scale_fill_brewer(name = "Region", palette = "Dark2", labels = c("Strong Upwelling","Weak Upwelling")) +
  scale_color_brewer(name = "Region",palette = "Dark2", labels = c("Strong Upwelling","Weak Upwelling")) +
  scale_shape_manual(values = c(23, 24, 21), name = "Mesocosm/Field", drop = FALSE,
                     labels = c("Current - Mesocosm","Future - Mesocosm","Field")) +
  guides(fill=guide_legend(override.aes=list(shape=21)))
PCA 

ggsave(plot = PCA, file = "Fig2_rev.png", 
       type = "cairo-png",  bg = "transparent",
       width = 30, height = 20, units = "cm", dpi = 600)


#################################################################################################################
#################################################################################################################
SoCA <- rbind(LB_day, CI_day)
SoDis <- subset(Discrete2, Trt != 'NC')
SoDis <- subset(SoDis, Trt != 'NF')
SoDis[nrow(SoDis) +1,] <- c("2021-02-21 16:00", "NA","N",NA,NA,NA)
SoDis$pH <- as.numeric(SoDis$pH)
SoDis$DO <- as.numeric(SoDis$DO)
SoDis$TC <- as.numeric(SoDis$TC)

NoCA <- rbind(VD_day, PA_day)
NoDis <- subset(Discrete2, Trt != 'SC')
NoDis <- subset(NoDis, Trt != 'SF')
NoDis[nrow(NoDis) +1,] <- c("2021-02-21 16:00", "NA","N",NA,NA,NA)
NoDis$pH <- as.numeric(NoDis$pH)
NoDis$DO <- as.numeric(NoDis$DO)
NoDis$TC <- as.numeric(NoDis$TC)

SoDis$Trt <- factor(SoDis$Trt, levels = c("SC","SF","N"))
SoCAFig <-ggplot(data=SoCA, aes(x=pH, y=DO, color = TC, fill = TC)) +
  geom_point(size = 3, position = "jitter", pch = 21, color = "#525252") +
  geom_point(size = 5, data = SoDis, aes(x = pH, y = DO, fill = TC, shape = Trt), color = "black") +
  scale_shape_manual(values = c(23, 24, 21), name = "Mesocosm/Field", drop = FALSE,
                    labels = c("Current - Mesocosm","Future - Mesocosm","Field"))+
  scale_fill_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24), guide = "none") +
  labs(y="Dissolved Oxygen (mg/L)", x="pH", title = "Weak Upwelling") +
  theme_classic() +
  xlim(7.4, 8.3) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") 
SoCAFig <- SoCAFig + #guides(fill = guide_colourbar(title.position="top",title.hjust =0.5)) + 
  scale_y_continuous(limits = c(2,10), breaks = c(2,4,6,8,10)) + 
  theme(text = element_text(size = 16))
SoCAFig

NoDis$Trt <- factor(NoDis$Trt, levels = c("NC","NF","N"))
NoCAFig <-ggplot(data=NoCA, aes(x=pH, y=DO, color = TC, fill = TC)) +
  geom_point(size = 3, position = "jitter", pch = 21, color = "#525252") +
  geom_point(size = 5, data = NoDis, aes(x = pH, y = DO, fill = TC, shape = Trt), color = "black") +
  scale_shape_manual(values = c(23, 24, 21), name = "Mesocosm/Field", drop = FALSE,
                     labels = c("Current - Mesocosm","Future - Mesocosm","Field"))+
  scale_color_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24), guide = "none") +
  scale_fill_gradientn(colours = pal, name = "Temp (\u00B0C)", limits = c(7, 24), guide = "none") +
  labs(y="Dissolved Oxygen (mg/L)", x=element_blank(), title = "Strong Upwelling") +
  theme_classic() +
  xlim(7.4, 8.3) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2.5)) +
  theme(legend.position = "none") 
NoCAFig <- NoCAFig + #guides(shape = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_y_continuous(limits = c(2,10), breaks = c(2,4,6,8,10)) + 
  theme(text = element_text(size = 16))
NoCAFig

Fig_scatter <- ggarrange(NoCAFig, SoCAFig, labels = c("(B)", "(C)"), ncol = 1, nrow = 2, 
                       common.legend = FALSE, align = "hv")
Fig_scatter

ggsave(plot = Fig1_rev, file = "POPUrchin_NoExp2.png", 
       type = "cairo-png",  bg = "white",
       width = 15, height = 25, units = "cm", dpi = 600)

##################################################################################################################
AllLB_Good$Site <- "Laguna Beach"
AllCI_Good$Site <- "Catalina Island"
AllVD_Good$Site <- "Van Damme"
AllPA_Good$Site <- "Point Arena"
AllLB_Good$Pop <- "A"
AllCI_Good$Pop <- "B"
AllVD_Good$Pop <- "A"
AllPA_Good$Pop <- "B"
All <- rbind(AllLB_Good, AllCI_Good, AllPA_Good, AllVD_Good)
All$Group1 <- paste(All$Dep1, All$Site)
All$Site <- factor(All$Site, levels=c('Point Arena','Van Damme','Catalina Island','Laguna Beach'))

pHSo <- ggplot(All, aes(x = DT, y = pH, group = Group1, alpha = Site)) + 
  geom_line(aes(color = Site)) +
  theme_classic() +
  labs(x= element_blank(), y = "pH", element_text(size = 20)) +
  scale_color_manual(values=c("#1b9e77","#1b9e77","#d95f02","#d95f02"), name = "Site") +
  scale_alpha_manual(values = c(1, 0.5, 1, 0.5), name = "Site") +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(),
        legend.position = c(0.9,0.85), legend.title = element_blank()) +
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  scale_y_continuous(limits = c(7.3, 8.4), breaks = c(7.4,7.6,7.8,8.0,8.2))
pHSo

TCSo <- ggplot(AllPA_Good, aes(x = DT, y = TC, group = Dep3)) + 
  geom_line(color = "#1b9e77") +
  geom_line(data = AllVD_Good, aes(x = DT, y = TC, group = Dep3), color = "#1b9e77", alpha = 0.5) +
  geom_line(data = AllCI_Good, aes(x = DT, y = TC, group = Dep3), color = "#d95f02") +
  geom_line(data = AllLB_Good, aes(x = DT, y = TC, group = Dep3), color = "#d95f02", alpha = 0.5) +
  theme_classic() +
  scale_color_manual(values=c("#1b9e77","#1b9e77","#d95f02","#d95f02"), name = "Site") +
  scale_alpha_manual(values = c(1, 0.5, 1, 0.5), name = "Site") +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "none")+
  scale_x_datetime(date_breaks = "6 months", expand = c(0, 0)) +
  ylim(7, 23)
TCSo 

DOSo <- ggplot(AllPA_Good, aes(x = DT, y = DO, group = Dep2)) + 
  geom_line(color = "#1b9e77") +
  geom_line(data = AllVD_Good, aes(x = DT, y = DO, group = Dep2), color = "#1b9e77", alpha = 0.5) +
  geom_line(data = AllCI_Good, aes(x = DT, y = DO, group = Dep2), color = "#d95f02") +
  geom_line(data = AllLB_Good, aes(x = DT, y = DO, group = Dep2), color = "#d95f02", alpha = 0.5) +
  scale_color_manual(values=c("#1b9e77","#1b9e77","#d95f02","#d95f02"), name = "Site") +
  scale_alpha_manual(values = c(1, 0.5, 1, 0.5), name = "Site") +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  theme_classic() +
  labs(y="Dissolved Oxygen (mg/L)", x=element_blank(), element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 10), breaks = c(3,5,7,9))+
  theme(legend.position = "none")
DOSo

So <- ggarrange(pHSo, TCSo, DOSo, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, 
                common.legend = FALSE, align = "hv")
So



All_FT <- annotate_figure(So,bottom = text_grob("Date (MM/YY)", size = 10)) 
All_FT

ggsave(plot = All_FT, file = "FigS1.png", 
       type = "cairo-png",  bg = "white",
       width = 30, height = 25, units = "cm", dpi = 300)

###################################################################################################################
sohobo = read.csv(paste('data/SoCal_dayavg_temp.csv'))
sohobo$DT <- as.POSIXct(sohobo$Time, format = "%d-%b-%Y", tz = "UTC")
sohobo.long <- gather(sohobo, site, temp, CTR:WPT, factor_key=TRUE)

CI_temp <- merge(AllCI_Good, sohobo, by.x = "DT", by.y = "DT",
                 all.y = TRUE)
CI_temp <- data.frame(CI_temp$Time, CI_temp$DT, "CI", CI_temp$TC)
colnames(CI_temp) <- c('Time','DT','site', 'temp')
sohobo.long <- rbind(sohobo.long, CI_temp)

LB_temp <- merge(AllLB_Good, sohobo, by.x = "DT", by.y = "DT",
                 all.y = TRUE)
LB_temp <- data.frame(LB_temp$Time, LB_temp$DT, "LB", LB_temp$TC)
colnames(LB_temp) <- c('Time','DT','site', 'temp')
sohobo.long <- rbind(sohobo.long, LB_temp)


Main.TC <- ggplot(sohobo.long, aes(x = DT, y = temp, group = site, color = site, alpha = site)) + 
  geom_line(size = 0.5) +
  scale_color_manual(labels = c("Christmas Tree Reef", "Point Vicente East", "White Point", "Catalina Island", "Laguna Beach"), values=c("#000000","#000000","#000000","#d95f02","#d95f02"), name = "Site") +
  scale_alpha_manual(labels = c("Christmas Tree Reef", "Point Vicente East", "White Point", "Catalina Island", "Laguna Beach"), values = c(1, 0.66, 0.33, 1, 0.5), name = "Site") +
  theme_classic() +
  theme(legend.position = c(0.9,0.9), legend.title = element_blank()) +
  labs(x= "Date (MM/YY)",element_text(size = 20)) +
  labs(y="Temperature (\u00B0C)", element_text(size = 20)) +
  scale_x_datetime(date_breaks = "6 months", date_labels = "%m/%y", expand = c(0, 0)) +
  ylim(10, 24)
Main.TC 

ggsave(plot = Main.TC, file = "FigS4.png", 
       type = "cairo-png",  bg = "white",
       width = 18, height = 12, units = "cm", dpi = 600)




