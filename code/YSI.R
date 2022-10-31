#Written by Emily Donham 05/04/21

# LoLinR is not on CRAN, so install with: 
#install_github('colin-olito/LoLinR')
###############################################################################################
###############################################################################################
library(plyr); library(dplyr);library(broom);library(ggplot2); library(lubridate);
library(LoLinR); library(stringr); library(lme4); library(lmerTest); library(multcomp); 
library(phytotools); library(googledrive); library(Rmisc);
library(tibble); library(ggpubr); library(wesanderson); library(tidyverse);library(vegan);
library(lsmeans); library(RLRsim); library(gridExtra); library(RColorBrewer);
###############################################################################################
###############################################################################################

rm(list = ls())

pH = read.csv('data/pH_processed.csv')
pH$DT <- force_tz(as.POSIXct((pH$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")
pHsubset = pH %>%
  filter(pH$DT > '2021-02-12 12:00')
pHsubset = pHsubset %>%
  filter(pHsubset$DT < '2021-05-15 12:00')

idx = which(pHsubset$QF_pH_1 == 1)
pH1 = pHsubset$pH_1[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_1[idx]
pH1_1 <- data.frame(cbind(rep('H1',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_2 == 1)
pH1 = pHsubset$pH_2[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_2[idx]
pH2 <- data.frame(cbind(rep('H2',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_3 == 1)
pH1 = pHsubset$pH_3[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_3[idx]
pH3 <- data.frame(cbind(rep('H3',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_4 == 1)
pH1 = pHsubset$pH_4[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_4[idx]
pH4 <- data.frame(cbind(rep('H4',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_5 == 1)
pH1 = pHsubset$pH_5[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_5[idx]
pH5 <- data.frame(cbind(rep('H5',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_6 == 1)
pH1 = pHsubset$pH_6[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_6[idx]
pH6 <- data.frame(cbind(rep('H6',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_7 == 1)
pH1 = pHsubset$pH_7[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_7[idx]
pH7 <- data.frame(cbind(rep('H7',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_8 == 1)
pH1 = pHsubset$pH_8[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_8[idx]
pH8 <- data.frame(cbind(rep('H8',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_9 == 1)
pH1 = pHsubset$pH_9[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_9[idx]
pH9 <- data.frame(cbind(rep('H9',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_10 == 1)
pH1 = pHsubset$pH_10[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_10[idx]
pH10 <- data.frame(cbind(rep('H10',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_11 == 1)
pH1 = pHsubset$pH_11[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_11[idx]
pH11 <- data.frame(cbind(rep('H11',length(pH1)),pH1SDN,pH1,pH1TC))

idx = which(pHsubset$QF_pH_12 == 1)
pH1 = pHsubset$pH_12[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_12[idx]
pH12 <- data.frame(cbind(rep('H12',length(pH1)),pH1SDN,pH1,pH1TC))

pH_long <- data.frame(rbind(pH1_1,pH2,pH3,pH4,pH5,pH6,pH7,pH8,pH9,pH10,pH11,pH12))
pH_long <- rename(pH_long, c("H" = "V1", "Date" = "pH1SDN", "pH" = "pH1", "Temp" = "pH1TC")) #Renmae columns
pH_long$Date <- as.numeric(pH_long$Date)
pH_long$Date <- force_tz(as.POSIXct((pH_long$Date - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")
options(digits = 4)
pH_long$pH <- as.numeric(pH_long$pH) # Converting to numeric from factor is odd...
pH_long$Temp <- as.numeric(pH_long$Temp)
pH_long$Trt <- "NA"

for (i in 1:nrow(pH_long)) {
  if (pH_long$H[i] == "H1") {
    pH_long$Trt[i] = "SC" } 
  else if (pH_long$H[i] == "H2") {
    pH_long$Trt[i] = "NF" }
  else if (pH_long$H[i] == "H3") {
    pH_long$Trt[i] = "SC" }
  else if (pH_long$H[i] == "H4") {
    pH_long$Trt[i] = "NF" }
  else if (pH_long$H[i] == "H5") {
    pH_long$Trt[i] = "SF" }
  else if (pH_long$H[i] == "H6") {
    pH_long$Trt[i] = "NC" }
  else if (pH_long$H[i] == "H7") {
    pH_long$Trt[i] = "SF" }
  else if (pH_long$H[i] == "H8") {
    pH_long$Trt[i] = "NF" }
  else if (pH_long$H[i] == "H9") {
    pH_long$Trt[i] = "SF" }
  else if (pH_long$H[i] == "H10") {
    pH_long$Trt[i] = "NC" }
  else if (pH_long$H[i] == "H11") {
    pH_long$Trt[i] = "SC" }
  else if (pH_long$H[i] == "H12") {
    pH_long$Trt[i] = "NC" }
}

chem = read.csv('data/AllDS.csv',
                colClasses=c('character','character','character','character',
                             'integer', 'character','character','character',
                             'double', 'double','double','double','double',
                             'double','double','double','double','double',
                             'double','double','double','double','double',
                             'double','double','double','double','double',
                             'double'))
chem$DT <- mdy_hm(paste(chem$DATE, chem$TIME)) # Creating R date/time stamp
chem$DT <- force_tz(as.POSIXct(chem$DT, origin = as.POSIXct("1970-01-01", TZ = "America/Los_Angeles"), 
                               TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles") #Super weird, but couldn't get it to not be in UTC, so had to do this
chem$Tco2 <- chem$HCO3 + chem$CO3 + chem$CO2 # Calculating total CO2

# Summary stats ran through to calc all for Bins
Bins = chem %>%
  filter(chem$SH == 'S')
pHave <- summarySE(data=Bins, measurevar = "War",
                   groupvars = c("TRT","TIMEPOINT","HEADER"))
pH2ave <- summarySE(data=pHave, measurevar = "War",
                    groupvars = c("TRT","TIMEPOINT"))
pH3ave <- summarySE(data=pH2ave, measurevar = "War",
                    groupvars = c("TRT"))

Headers = chem %>%    # Subset only header samples 
  filter(SH == "H") 
Headers$ID <- factor(Headers$ID)    #Need to reset the levels of the factors, for looping
u <- data.frame(unique(Headers$ID))    #This will be used to loop through all headers
u <- rename(u, c("ID" = "unique.Headers.ID."))
u$ID <- factor(u$ID)
offset <- data.frame(H1 = NA, H2 = NA, H3 = NA, H4 = NA, H5 = NA, H6 = NA, H7 = NA,
                     H8 = NA, H9 = NA, H10 = NA, H11 = NA, H12 = NA)
for (i in 1:nrow(u)) {
  temp = Headers %>%
    filter(ID == u$ID[i])
  for (j in 1:nrow(temp)) {
    temp2 = pH_long %>%
      filter(H == u$ID[i])
    ind <- which.min(abs(temp$DT[j]-temp2$Date))
    offset[j, i] <- temp2$pH[ind] - temp$pH.out[j]
  }
}   
offset <- cbind(offset, unique(Headers$DATE), unique(Headers$TIME)) # combine dates to time series of discrete samples
offset <-rename(offset, c("Date" = "unique(Headers$DATE)"))
offset <-rename(offset, c("DT" = "unique(Headers$TIME)"))

# Need to restructure to plot
offsetP <- data.frame(rbind(cbind(rep('H1',nrow(offset)), offset$H1, offset$Date, offset$DT), cbind(rep('H2',nrow(offset)), offset$H2, offset$Date, offset$DT),
                            cbind(rep('H3',nrow(offset)), offset$H3, offset$Date, offset$DT), cbind(rep('H4',nrow(offset)), offset$H4, offset$Date, offset$DT),
                            cbind(rep('H5',nrow(offset)), offset$H5, offset$Date, offset$DT), cbind(rep('H6',nrow(offset)), offset$H6, offset$Date, offset$DT),
                            cbind(rep('H7',nrow(offset)), offset$H7, offset$Date, offset$DT), cbind(rep('H8',nrow(offset)), offset$H8, offset$Date, offset$DT),
                            cbind(rep('H9',nrow(offset)), offset$H9, offset$Date, offset$DT), cbind(rep('H10',nrow(offset)), offset$H10, offset$Date, offset$DT),
                            cbind(rep('H11',nrow(offset)), offset$H11, offset$Date, offset$DT), cbind(rep('H12',nrow(offset)), offset$H12, offset$Date, offset$DT)))
offsetP$X2 <- as.numeric(offsetP$X2)
offsetP$DT <- mdy_hm(paste(offsetP$X3, offsetP$X4)) # Creating R date/time stamp
offsetP$DT <- force_tz(as.POSIXct(offsetP$DT, origin = as.POSIXct("1970-01-01", TZ = "America/Los_Angeles"), 
                               TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")

offsetM <- summarySE(data=offsetP, measurevar = "X2",
                     groupvars = "X1", na.rm = TRUE)
for(i in 1:nrow(pH_long)) {
  if(pH_long$H[i] =='H1') {
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[1]
  }else if (pH_long$H[i] =='H10'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[2]
  }else if (pH_long$H[i] =='H11'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[3]
  }else if (pH_long$H[i] =='H12'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[4]
  }else if (pH_long$H[i] =='H2'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[5]
  }else if (pH_long$H[i] =='H3'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[6]
  }else if (pH_long$H[i] =='H4'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[7]
  }else if (pH_long$H[i] =='H5'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[7]
  }else if (pH_long$H[i] =='H6'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[9]
  }else if (pH_long$H[i] =='H7'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[10]
  }else if (pH_long$H[i] =='H8'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[11]
  }else if (pH_long$H[i] =='H9'){
    pH_long$pHcorr2[i] = pH_long$pH[i] - offsetM$X2[12]
  }else {
    pH_long$pHcorr2[i] = pH_long$pH[i]
  }
}   

offsetP$SDN <- as.numeric(offsetP$DT)

for (i in 1:nrow(u)) {
  temp = offsetP %>%
    filter(X1 == u$ID[i])
  lmU <- lm(temp$X2 ~ temp$SDN)
  temp <- na.omit(temp)
  polyU <- lm(temp$X2 ~ poly(temp$SDN, 2, raw = TRUE))
  u$intLM[i] <- coef(lmU)["(Intercept)"]
  u$slope[i] <- coef(lmU)["temp$SDN"]
  u$intPOLY[i] <- coef(polyU)["(Intercept)"]
  u$poly1[i] <- coef(polyU)["poly(temp$SDN, 2, raw = TRUE)1"]
  u$poly2[i] <- coef(polyU)["poly(temp$SDN, 2, raw = TRUE)2"]
}   

pH_long$DateNum <- as.numeric(pH_long$Date)
for (i in 1:nrow(u)) {
    t <- which(pH_long$H == u$ID[i])
  for (j in 1:length(t)){
    pH_long$corrP[t[j]] <-  (pH_long$DateNum[t[j]]^2 * u$poly1[i] + pH_long$DateNum[t[j]] * u$poly2[i] + u$intPOLY[i])
} }


# Plot offsets over time for all headers
i <- ggplot(offsetP, aes(x = SDN, y = X2, color = X1)) +
  geom_point() +
  theme_classic() +
  geom_smooth(data = offsetP, method = lm, formula = y ~ poly(x, 2), se = FALSE) +
  facet_wrap(~offsetP$X1) +
  scale_color_manual(values = wes_palette("Zissou1", 12, "continuous")) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title="Durafet Calibration",y="Offset", x = "Time Point")
i 
ggsave(plot = i, file = "DurafetCalWinter2021.png", 
       type = "cairo-png",  bg = "transparent",
       width = 10, height = 10, units = "cm", dpi = 300)

DO = read.csv('data/DO_processed_edit.csv')
# Create R date/time
DO$DT <- force_tz(as.POSIXct((DO$SDN - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")
DO$DT <- format(DO$DT, format = '%Y-%m-%d %H:%M')
DOsubset = DO %>%
  filter(DO$DT > '2021-02-12 12:00')
DOsubset = DOsubset %>%
  filter(DOsubset$DT < '2021-05-15 12:00')

DO_long <- data.frame(rbind(cbind(rep('H1',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_1, DOsubset$TC_1, DOsubset$QF_DO_1, DOsubset$QF_DO2_1), cbind(rep('H2',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_2, DOsubset$TC_2,DOsubset$QF_DO_2, DOsubset$QF_DO2_2),
                            cbind(rep('H3',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_3, DOsubset$TC_3, DOsubset$QF_DO_3, DOsubset$QF_DO2_3), cbind(rep('H4',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_4, DOsubset$TC_4, DOsubset$QF_DO_4, DOsubset$QF_DO2_4),
                            cbind(rep('H5',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_5, DOsubset$TC_5, DOsubset$QF_DO_5, DOsubset$QF_DO2_5), cbind(rep('H6',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_6, DOsubset$TC_6,DOsubset$QF_DO_6, DOsubset$QF_DO2_6),
                            cbind(rep('H7',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_7, DOsubset$TC_7, DOsubset$QF_DO_7, DOsubset$QF_DO2_7), cbind(rep('H8',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_8, DOsubset$TC_8,DOsubset$QF_DO_8, DOsubset$QF_DO2_8),
                            cbind(rep('H9',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_9, DOsubset$TC_9,DOsubset$QF_DO_9, DOsubset$QF_DO2_9),cbind(rep('H10',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_10, DOsubset$TC_10,DOsubset$QF_DO_10, DOsubset$QF_DO2_10),
                            cbind(rep('H11',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_11, DOsubset$TC_11, DOsubset$QF_DO_11, DOsubset$QF_DO2_11), cbind(rep('H12',nrow(DOsubset)), DOsubset$SDN, DOsubset$mgL_12, DOsubset$TC_12, DOsubset$QF_DO_12, DOsubset$QF_DO2_12)))


DO_long <- rename(DO_long, c("H" = "X1", "Date" = "X2", "DO" = "X3", "Temp" = "X4")) #Rename columns
DO_long$Date <- as.numeric(DO_long$Date)
DO_long$Date <- force_tz(as.POSIXct((DO_long$Date - 719529)*86400, origin = "1970-01-01", TZ = "America/Los_Angeles"), tzone = "America/Los_Angeles")
DO_long$DO <- as.numeric(DO_long$DO) # Converting to numeric from factor is odd...
DO_long$Temp <- as.numeric(DO_long$Temp)

DO_long = DO_long %>%
  filter(DO_long$X5 == 1)
DO_long = DO_long %>%
  filter(DO_long$X6 == 0)

for (i in 1:nrow(DO_long)) {
  if (DO_long$H[i] == "H1") {
    DO_long$Trt[i] = "SC" } 
  else if (DO_long$H[i] == "H2") {
    DO_long$Trt[i] = "NF" }
  else if (DO_long$H[i] == "H3") {
    DO_long$Trt[i] = "SC" }
  else if (DO_long$H[i] == "H4") {
    DO_long$Trt[i] = "NF" }
  else if (DO_long$H[i] == "H5") {
    DO_long$Trt[i] = "SF" }
  else if (DO_long$H[i] == "H6") {
    DO_long$Trt[i] = "NC" }
  else if (DO_long$H[i] == "H7") {
    DO_long$Trt[i] = "SF" }
  else if (DO_long$H[i] == "H8") {
    DO_long$Trt[i] = "NF" }
  else if (DO_long$H[i] == "H9") {
    DO_long$Trt[i] = "SF" }
  else if (DO_long$H[i] == "H10") {
    DO_long$Trt[i] = "NC" }
  else if (DO_long$H[i] == "H11") {
    DO_long$Trt[i] = "SC" }
  else if (DO_long$H[i] == "H12") {
    DO_long$Trt[i] = "NC" }
}

idx = which(pHsubset$QF_pH_1 == 1)
pH1 = pHsubset$pH_1[idx]
pH1SDN = pHsubset$SDN[idx]
pH1TC = pHsubset$TC_1[idx]
pH1_1 <- data.frame(cbind(rep('H1',length(pH1)),pH1SDN,pH1,pH1TC))



YSI = read.csv('data/YSI_Winter2021.csv')
# Convert time
YSI$DT <- force_tz(as.POSIXct(with(YSI, mdy(Date)), origin = as.POSIXct("1970-01-01",TZ = "America/Los_Angeles")), tzone = "America/Los_Angeles")
YSI$Trt <- rep(c('SC','SC','NF','NF','SC','SC','NF','NF','SF','SF','NC','NC',
                 'SF','SF','NF','NF','SF','SF','NC','NC','SC','SC','NC','NC'), 97 )
YSI$pH <- YSI$pH - 0.2
YSI$Tub <- as.character((YSI$Tub))
YSI$Header <- substr(YSI$Tub, 1, nchar(YSI$Tub) -1) 

pHave <- summarySE(data=YSI, measurevar = "pH",
                   groupvars = c("Header","DT"))

pHave$Trt <- c(rep('SC',93),rep('NC',93),rep('SC',93),rep('NC',93),
               rep('NF',93),rep('SC',93),rep('NF',93),rep('SF',93),rep('NC',93),
               rep('SF',93),rep('NF',93),rep('SF',93))

TCave <- summarySE(data=YSI, measurevar = "Temp",
                   groupvars = c("Header","DT"))

TCave$Trt <- c(rep('SC',93),rep('NC',93),rep('SC',93),rep('NC',93),
               rep('NF',93),rep('SC',93),rep('NF',93),rep('SF',93),rep('NC',93),
               rep('SF',93),rep('NF',93),rep('SF',93))

TC2ave <- summarySE(data=TCave, measurevar = "Temp",
                    groupvars = c("Trt", "DT"))

TC3ave <- summarySE(data=TC2ave, measurevar = "Temp",
                    groupvars = c("Trt"))

DOave <- summarySE(data=YSI, measurevar = "DO..mg.L.",
                   groupvars = c("Header","DT"))

DOave$Trt <- c(rep('SC',93),rep('NC',93),rep('SC',93),rep('NC',93),
               rep('NF',93),rep('SC',93),rep('NF',93),rep('SF',93),rep('NC',93),
               rep('SF',93),rep('NF',93),rep('SF',93))

DO2ave <- summarySE(data=DOave, measurevar = "DO..mg.L.",
                    groupvars = c("Trt", "DT"))

DO3ave <- summarySE(data=DO2ave, measurevar = "DO..mg.L.",
                    groupvars = c("Trt"))

Salave <- summarySE(data=YSI, measurevar = "Salinity..ppt.",
                   groupvars = c("Header","DT"))

Salave$Trt <- c(rep('SC',93),rep('NC',93),rep('SC',93),rep('NC',93),
               rep('NF',93),rep('SC',93),rep('NF',93),rep('SF',93),rep('NC',93),
               rep('SF',93),rep('NF',93),rep('SF',93))

Sal2ave <- summarySE(data=Salave, measurevar = "Salinity..ppt.",
                    groupvars = c("Trt", "DT"))

Sal3ave <- summarySE(data=Sal2ave, measurevar = "Salinity..ppt.",
                    groupvars = c("Trt"))



i <- ggplot(TCave, aes(x = DT, y = Temp, color = Trt)) +
  geom_point(size = 2) +
  geom_line(pH_long, mapping = aes(x = Date, y = Temp, color = Trt), alpha = 1/3) +
  theme_classic() +
  scale_color_brewer(palette = "Dark2", labels=c("Strong Upwelling Current", "Strong Upwelling Future", 
                                                 "Weak Upwelling Current", "Weak Upwelling Future")) +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Temp (C)", x = "Date") +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "bottom")+
  ylim(9.5, 19.5) + 
  scale_y_continuous(breaks = c(10,13,16,19))
i
ggsave(plot = i, file = "YSI_MesoWinter2021Temp.png", 
       type = "cairo-png",  bg = "transparent",
       width = 18, height = 5, units = "cm", dpi = 300)

chem2 <- subset(chem, chem$TRT != 'CRM')
chem2 <- subset(chem2, chem2$SH != 'S')

chem2$TRT <- factor(chem2$TRT, levels = c("NC","NF","SC","SF"))

k <- ggplot(chem2, aes(x = DT, y = pH.out, color = TRT)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_line(pH_long, mapping = aes(x = Date, y = pHcorr2, color = Trt), alpha = 1/3) +
  scale_color_brewer(palette = "Dark2", labels=c("Strong Upwelling Current", "Strong Upwelling Future", 
                                                 "Weak Upwelling Current", "Weak Upwelling Future")) +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="pH", x = "Date") +
  theme(axis.title.x=element_blank(), axis.text.x = element_blank()) +
  theme(legend.position = "bottom")+
  ylim(7.4, 8.1)
k
ggsave(plot = k, file = "YSI_MesoWinter2021AllDS.png", 
       type = "cairo-png",  bg = "transparent",
       width = 18, height = 5, units = "cm", dpi = 300)

j <- ggplot(DOave, aes(x = DT, y = DO..mg.L., color = Trt)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_line(DO_long, mapping = aes(x = Date, y = DO, color = Trt), alpha = 1/3) +
  scale_color_brewer(palette = "Dark2", labels=c("Strong Upwelling Current", "Strong Upwelling Future", 
                                                 "Weak Upwelling Current", "Weak Upwelling Future")) +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(y="Dissolved Oxygen (mg/L)", x = "Date") +
  ylim(3, 13) +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = c(3,7,10,13))

j
ggsave(plot = j, file = "YSI_MesoWinter2021DO.png", 
       type = "cairo-png",  bg = "transparent",
       width = 18, height = 5, units = "cm", dpi = 300)

l <- ggarrange(k, i, j, labels = c("(a)", "(b)", "(c)"), ncol = 1, nrow = 3, align = "hv",
               common.legend = TRUE, legend = "bottom")
l

ggsave(plot = l, file = "YSI_MesoWinter.png", 
       type = "cairo-png",  bg = "white",
       width = 23, height = 21, units = "cm", dpi = 300)

DS = read.csv('AllDS.csv')
# Convert time
DS$DT <- force_tz(as.POSIXct(with(DS, mdy(DATE)), origin = as.POSIXct("1970-01-01",TZ = "America/Los_Angeles")), tzone = "America/Los_Angeles")

# based on variable values

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

Discrete <- merge(pHave3, TCave3, by = "Trt")
Discrete <- merge(Discrete, DOave3, by = "Trt")




