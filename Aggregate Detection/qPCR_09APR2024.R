library(ggplot2)
library(forcats)
library(tidyverse)
library(xml2)
library(rsvg)
library(svglite)
library(reshape2)
library(grid)
library(gridExtra)
library(ggpubr)
library(xlsx)
library(dplyr)
library(stringr)
library(readr)

melt_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-04-SS/Aggregate Probe SS 04APR2024_FULL_Melt_09APR2024.csv")
sampIDs = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-04-SS/Aggregate Probe SS 04APR2024_FULL_IDs.csv")
standardCurve = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-04-SS/Aggregate Probe SS 04APR2024_FULL_StandardCurve.csv")
melt_file_A1 <- melt_file[which(melt_file$Well.Position == "A1"),]

melt_file_parse <- melt_file[which(melt_file$Well.Position == c("B2","B4","B6")),]
min_temp <- round(min(melt_file_parse$Temperature))
max_temp <- round(max(melt_file_parse$Temperature))
temps <- seq(min_temp,max_temp,by=1)

melt_file_parse$Fluorescence <- gsub(",","",melt_file_parse$Fluorescence)
melt_file_parse$Fluorescence <- as.numeric(as.character(melt_file_parse$Fluorescence))

min_fluor <- min(melt_file_parse$Fluorescence)
max_fluor <- max(melt_file_parse$Fluorescence)
fluor <- round(seq(min_fluor,max_fluor,length.out=10))

plot_melt <- ggplot(melt_file_parse) + 
  geom_line(aes(x = Temperature, y=Fluorescence, group=Well.Position,color=Well.Position),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=Fluorescence, color=Well.Position)) + 
  scale_y_continuous(breaks = fluor, labels = as.character(fluor)) +
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
plot_melt

##################################################################################
##################################################################################
##################################################################################

melt_file_parse <- melt_file[which(melt_file$Well.Position == c("B1")),]
min_temp <- round(min(melt_file_parse$Temperature))
max_temp <- round(max(melt_file_parse$Temperature))
temps <- seq(min_temp,max_temp,by=1)

melt_file_parse$Fluorescence <- gsub(",","",melt_file_parse$Fluorescence)
melt_file_parse$Fluorescence <- as.numeric(as.character(melt_file_parse$Fluorescence))

min_fluor <- min(melt_file_parse$Fluorescence)
max_fluor <- max(melt_file_parse$Fluorescence)
fluor <- round(seq(min_fluor,max_fluor,length.out=10))

plot_melt <- ggplot(melt_file_parse) + 
  geom_line(aes(x = Temperature, y=Fluorescence, group=Target.Name,color=Target.Name),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=Fluorescence, color=Target.Name)) + 
  scale_y_continuous(breaks = fluor, labels = as.character(fluor)) +
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
plot_melt

##################################################################################
############################AGGREGATE PERCENTAGE##################################
##################################################################################

melt_file_A1 <- melt_file[which(melt_file$Well.Position == c("A1")),]
melt_file_A1$Fluorescence <- gsub(",","",melt_file_A1$Fluorescence)
melt_file_A1$Fluorescence <- as.numeric(as.character(melt_file_A1$Fluorescence))
melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("J18")),]
melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

percAGG <- melt_file_A1$Fluorescence / melt_file_AGG$Fluorescence * 100

AGG <- data.frame(Temperature = melt_file_A1$Temperature, 
                  percAGG = percAGG,
                  Target.Name = melt_file_A1$Target.Name)

min_temp <- round(min(melt_file_A1$Temperature))
max_temp <- round(max(melt_file_A1$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(AGG) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Target.Name),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Target.Name)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                  subtitle = "SYPRO")
plot_melt





# AGGREGATE PERCENTAGE FROM FREEZE
##############################DP 1.9E11#############################################
sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])
# A = which(inx %in% c("B1","B2","B3","B4","B5","B6"))
# A = which(inx %in% c("B7","B8","B9","B10","B11","B12"))

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep

# melt_file_B <- melt_file[which(melt_file$Well.Position %in% c("B2","B4","B6")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B8","B10","B12")),]
# melt_file_B <- melt_file[which(melt_file$Well.Position == c("B14","B16","B18")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B20","B22","B24")),]

melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B20","B22","B24")),]



sampIDs_B <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_B_rep <- rep(sampIDs_B$Sample[1:8], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_B$Well[1:8])
# A = which(inx %in% c("B1","B2","B3","B4","B5","B6"))
# A = which(inx %in% c("B7","B8","B9","B10","B11","B12"))

melt_file_B <- melt_file[A,]
melt_file_B$Sample <- sampIDs_B_rep
  
# melt_file_B <- melt_file[which(melt_file$Well.Position %in% c("B2","B4","B6")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B8","B10","B12")),]
# melt_file_B <- melt_file[which(melt_file$Well.Position == c("B14","B16","B18")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B20","B22","B24")),]

melt_file_B$Fluorescence <- gsub(",","",melt_file_B$Fluorescence)
melt_file_B$Fluorescence <- as.numeric(as.character(melt_file_B$Fluorescence))

melt_file_B <- melt_file_B %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_B_plot <- melt_file_B[which(melt_file_B$Well.Position %in% c("B2","B4","B6")),]

#subtract background
melt_file_B_plot$Fluorescence <- melt_file_B_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

B = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(B,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_B_plot$Fluorescence[i],
                                    standardCurve_slope[i],
                                    standardCurve_intercept[i])
}
melt_file_B_plot$percAGG <- percAGG


#melt_file_B_plot$percAGG <- melt_file_B_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_B_plot$Sample,
#                   Temperature = melt_file_B_plot$Temperature, 
#                   percAGG = melt_file_B_plot$percAGG,
#                   Well.Position = melt_file_B_plot$Well.Position,
#                   Target.Name = melt_file_B_plot$Target.Name)

min_temp <- round(min(melt_file_B_plot$Temperature))
max_temp <- round(max(melt_file_B_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_B_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 1.9E11 | PROTEOSTAT | FREEZE")+
  labs(color="Sample")
  # scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
  #                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt

##############################DP 3.5E11######################################
sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])
# A = which(inx %in% c("B1","B2","B3","B4","B5","B6"))
# A = which(inx %in% c("B7","B8","B9","B10","B11","B12"))

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep

# melt_file_B <- melt_file[which(melt_file$Well.Position %in% c("B2","B4","B6")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B8","B10","B12")),]
# melt_file_B <- melt_file[which(melt_file$Well.Position == c("B14","B16","B18")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B20","B22","B24")),]

melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B20","B22","B24")),]

sampIDs_B <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_B_rep <- rep(sampIDs_B$Sample[7:12], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_B$Well[7:12])
# A = which(inx %in% c("B1","B2","B3","B4","B5","B6"))
# A = which(inx %in% c("B7","B8","B9","B10","B11","B12"))

melt_file_B <- melt_file[A,]
melt_file_B$Sample <- sampIDs_B_rep

# melt_file_B <- melt_file[which(melt_file$Well.Position %in% c("B2","B4","B6")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B8","B10","B12")),]
# melt_file_B <- melt_file[which(melt_file$Well.Position == c("B14","B16","B18")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B20","B22","B24")),]

melt_file_B$Fluorescence <- gsub(",","",melt_file_B$Fluorescence)
melt_file_B$Fluorescence <- as.numeric(as.character(melt_file_B$Fluorescence))

melt_file_B <- melt_file_B %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_B_plot <- melt_file_B[which(melt_file_B$Well.Position %in% c("B8","B10","B12")),]

#subtract background
melt_file_B_plot$Fluorescence <- melt_file_B_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

B = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(B,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_B_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_B_plot$percAGG <- percAGG


#melt_file_B_plot$percAGG <- melt_file_B_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_B_plot$Sample,
#                   Temperature = melt_file_B_plot$Temperature, 
#                   percAGG = melt_file_B_plot$percAGG,
#                   Well.Position = melt_file_B_plot$Well.Position,
#                   Target.Name = melt_file_B_plot$Target.Name)

min_temp <- round(min(melt_file_B_plot$Temperature))
max_temp <- round(max(melt_file_B_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_B_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 3.5E11 | PROTEOSTAT | FREEZE")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt

##############################DS3##################################################
sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B20","B22","B24")),]


sampIDs_B <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_B_rep <- rep(sampIDs_B$Sample[13:18], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_B$Well[13:18])
# A = which(inx %in% c("B1","B2","B3","B4","B5","B6"))
# A = which(inx %in% c("B7","B8","B9","B10","B11","B12"))

melt_file_B <- melt_file[A,]
melt_file_B$Sample <- sampIDs_B_rep

# melt_file_B <- melt_file[which(melt_file$Well.Position %in% c("B2","B4","B6")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B8","B10","B12")),]
# melt_file_B <- melt_file[which(melt_file$Well.Position == c("B14","B16","B18")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B20","B22","B24")),]

melt_file_B$Fluorescence <- gsub(",","",melt_file_B$Fluorescence)
melt_file_B$Fluorescence <- as.numeric(as.character(melt_file_B$Fluorescence))

melt_file_B <- melt_file_B %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_B_plot <- melt_file_B[which(melt_file_B$Well.Position %in% c("B14","B16","B18")),]

#subtract background
melt_file_B_plot$Fluorescence <- melt_file_B_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

B = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(B,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))


interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_B_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_B_plot$percAGG <- percAGG

#melt_file_B_plot$percAGG <- melt_file_B_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_B_plot$Sample,
#                   Temperature = melt_file_B_plot$Temperature, 
#                   percAGG = melt_file_B_plot$percAGG,
#                   Well.Position = melt_file_B_plot$Well.Position,
#                   Target.Name = melt_file_B_plot$Target.Name)

min_temp <- round(min(melt_file_B_plot$Temperature))
max_temp <- round(max(melt_file_B_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_B_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DS3 1E12 | PROTEOSTAT | FREEZE")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt
##############################BUFFER##################################################
sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])
# A = which(inx %in% c("B1","B2","B3","B4","B5","B6"))
# A = which(inx %in% c("B7","B8","B9","B10","B11","B12"))

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep

# melt_file_B <- melt_file[which(melt_file$Well.Position %in% c("B2","B4","B6")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B8","B10","B12")),]
# melt_file_B <- melt_file[which(melt_file$Well.Position == c("B14","B16","B18")),]
# #melt_file_B <- melt_file[which(melt_file$Well.Position == c("B20","B22","B24")),]

melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_plot <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B20","B22","B24")),]

B = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(B,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_B_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_B_plot$percAGG <- percAGG



# AGG <- data.frame(Sample = melt_file_B_plot$Sample,
#                   Temperature = melt_file_B_plot$Temperature, 
#                   percAGG = melt_file_B_plot$percAGG,
#                   Well.Position = melt_file_B_plot$Well.Position,
#                   Target.Name = melt_file_B_plot$Target.Name)

min_temp <- round(min(melt_file_B_plot$Temperature))
max_temp <- round(max(melt_file_B_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)

# melt_file_B_plot$Sample <- factor(melt_file_B_plot$Sample, levels=c("FFB2 0% PS20", "FFB 0.001% PS20", "FFB 0.008% PS20"))

plot_melt <- ggplot(melt_file_B_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "Sample Buffers | PROTEOSTAT")+
  labs(color="Sample")
  #scale_color_manual(values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt




# AGGREGATE PERCENTAGE FROM 4C HOLD


# AGGREGATE PERCENTAGE FROM 4C HOLD
##############################DP 1.9E11#############################################
sampIDs_BUFFER <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("G20","G22","G24")),]

sampIDs_G <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_G_rep <- rep(sampIDs_G$Sample[1:8], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_G$Well[1:8])
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_G <- melt_file[A,]
melt_file_G$Sample <- sampIDs_G_rep


melt_file_G$Fluorescence <- gsub(",","",melt_file_G$Fluorescence)
melt_file_G$Fluorescence <- as.numeric(as.character(melt_file_G$Fluorescence))

melt_file_G <- melt_file_G %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_G_plot <- melt_file_G[which(melt_file_G$Well.Position %in% c("G2","G4","G6")),]

#subtract background
melt_file_G_plot$Fluorescence <- melt_file_G_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

G = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(G,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_G_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_G_plot$percAGG <- percAGG


#melt_file_G_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_G_plot$Temperature))
max_temp <- round(max(melt_file_G_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_G_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 1.9E11 | PROTEOSTAT | 4C Hold")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt

##############################DP 3.5E11######################################

sampIDs_BUFFER <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("G20","G22","G24")),]

sampIDs_G <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_G_rep <- rep(sampIDs_G$Sample[7:12], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_G$Well[7:12])
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_G <- melt_file[A,]
melt_file_G$Sample <- sampIDs_G_rep

# melt_file_G <- melt_file[which(melt_file$Well.Position %in% c("G2","G4","G6")),]
# #melt_file_G <- melt_file[which(melt_file$Well.Position == c("G8","G10","G12")),]
# melt_file_G <- melt_file[which(melt_file$Well.Position == c("G14","G16","G18")),]
# #melt_file_G <- melt_file[which(melt_file$Well.Position == c("G20","G22","G24")),]

melt_file_G$Fluorescence <- gsub(",","",melt_file_G$Fluorescence)
melt_file_G$Fluorescence <- as.numeric(as.character(melt_file_G$Fluorescence))

melt_file_G <- melt_file_G %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_G_plot <- melt_file_G[which(melt_file_G$Well.Position %in% c("G8","G10","G12")),]

#subtract background
melt_file_G_plot$Fluorescence <- melt_file_G_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

G = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(G,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_G_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_G_plot$percAGG <- percAGG


#melt_file_G_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_G_plot$Temperature))
max_temp <- round(max(melt_file_G_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_G_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 3.5E11 | PROTEOSTAT | 4C Hold")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt

##############################DS3##################################################
sampIDs_BUFFER <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("G22","G24")),]

sampIDs_G <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_G_rep <- rep(sampIDs_G$Sample[13:18], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_G$Well[13:18])
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_G <- melt_file[A,]
melt_file_G$Sample <- sampIDs_G_rep


melt_file_G$Fluorescence <- gsub(",","",melt_file_G$Fluorescence)
melt_file_G$Fluorescence <- as.numeric(as.character(melt_file_G$Fluorescence))

melt_file_G <- melt_file_G %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_G_plot <- melt_file_G[which(melt_file_G$Well.Position %in% c("G16","G18")),] #DS3 0% not enough sample

#subtract background
melt_file_G_plot$Fluorescence <- melt_file_G_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

G = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(G,2),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))


interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],2) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],2)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_G_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_G_plot$percAGG <- percAGG

#melt_file_G_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_G_plot$Temperature))
max_temp <- round(max(melt_file_G_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_G_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DS3 1E12 | PROTEOSTAT | 4C Hold")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt
##############################BUFFER##################################################
sampIDs_G <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_G_rep <- rep(sampIDs_G$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_G$Well[19:24])
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_G <- melt_file[A,]
melt_file_G$Sample <- sampIDs_G_rep

# melt_file_G <- melt_file[which(melt_file$Well.Position %in% c("G2","G4","G6")),]
# #melt_file_G <- melt_file[which(melt_file$Well.Position == c("G8","G10","G12")),]
# melt_file_G <- melt_file[which(melt_file$Well.Position == c("G14","G16","G18")),]
# #melt_file_G <- melt_file[which(melt_file$Well.Position == c("G20","G22","G24")),]

melt_file_G$Fluorescence <- gsub(",","",melt_file_G$Fluorescence)
melt_file_G$Fluorescence <- as.numeric(as.character(melt_file_G$Fluorescence))

melt_file_G <- melt_file_G %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_G_plot <- melt_file_G[which(melt_file_G$Well.Position %in% c("G20","G22","G24")),]

G = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(G,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_G_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_G_plot$percAGG <- percAGG



# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_G_plot$Temperature))
max_temp <- round(max(melt_file_G_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)

# melt_file_G_plot$Sample <- factor(melt_file_G_plot$Sample, levels=c("FFG2 0% PS20", "FFG 0.001% PS20", "FFG 0.008% PS20"))

plot_melt <- ggplot(melt_file_G_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "Sample Buffers | PROTEOSTAT | 4C Hold")+
  labs(color="Sample")
#scale_color_manual(values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt



# AGGREGATE PERCENTAGE ROOM TEMP
##############################DP 1.9E11#############################################
sampIDs_BUFFER <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("G20","G22","G24")),]

sampIDs_K <- sampIDs[grepl('K',sampIDs$Well, fixed=T),]
sampIDs_K_rep <- rep(sampIDs_K$Sample[1:8], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_K$Well[1:8])
# A = which(inx %in% c("K1","K2","K3","K4","K5","K6"))
# A = which(inx %in% c("K7","K8","K9","K10","K11","K12"))

melt_file_K <- melt_file[A,]
melt_file_K$Sample <- sampIDs_K_rep

# melt_file_K <- melt_file[which(melt_file$Well.Position %in% c("K2","K4","K6")),]
# #melt_file_K <- melt_file[which(melt_file$Well.Position == c("K8","K10","K12")),]
# melt_file_K <- melt_file[which(melt_file$Well.Position == c("K14","K16","K18")),]
# #melt_file_K <- melt_file[which(melt_file$Well.Position == c("K20","K22","K24")),]

melt_file_K$Fluorescence <- gsub(",","",melt_file_K$Fluorescence)
melt_file_K$Fluorescence <- as.numeric(as.character(melt_file_K$Fluorescence))

melt_file_K <- melt_file_K %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_K_plot <- melt_file_K[which(melt_file_K$Well.Position %in% c("K2","K4","K6")),]

#subtract background
melt_file_K_plot$Fluorescence <- melt_file_K_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_K_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_K_plot$percAGG <- percAGG


#melt_file_G_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_K_plot$Temperature))
max_temp <- round(max(melt_file_K_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_K_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 1.9E11 | PROTEOSTAT | ROOM TEMP")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt

##############################DP 3.5E11######################################
sampIDs_BUFFER <- sampIDs[grepl('G',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[19:24], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[19:24])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("G20","G22","G24")),]

sampIDs_K <- sampIDs[grepl('K',sampIDs$Well, fixed=T),]
sampIDs_K_rep <- rep(sampIDs_K$Sample[7:12], each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_K$Well[7:12])
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_K <- melt_file[A,]
melt_file_K$Sample <- sampIDs_K_rep

# melt_file_K <- melt_file[which(melt_file$Well.Position %in% c("G2","G4","G6")),]
# #melt_file_K <- melt_file[which(melt_file$Well.Position == c("G8","G10","G12")),]
# melt_file_K <- melt_file[which(melt_file$Well.Position == c("G14","G16","G18")),]
# #melt_file_G <- melt_file[which(melt_file$Well.Position == c("K20","K22","K24")),]

melt_file_K$Fluorescence <- gsub(",","",melt_file_K$Fluorescence)
melt_file_K$Fluorescence <- as.numeric(as.character(melt_file_K$Fluorescence))

melt_file_K <- melt_file_K %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_K_plot <- melt_file_K[which(melt_file_K$Well.Position %in% c("K8","K10","K12")),]

#subtract background
melt_file_K_plot$Fluorescence <- melt_file_K_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_K_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_K_plot$percAGG <- percAGG


#melt_file_K_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_K_plot$Temperature))
max_temp <- round(max(melt_file_K_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_K_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 3.5E11 | PROTEOSTAT | ROOM TEMP")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt


##################################################################################
############################AGGREGATE PERCENTAGE##################################
##############################    HOLD COND       ################################

# AGGREGATE PERCENTAGE AT CONDS FOR 1.9E11
##############################DP 1.9E11 0%##########################################
sampIDs_BUFFER_A <- sampIDs[grepl('19',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_B <- sampIDs[grepl('20',sampIDs$Well, fixed=T),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER_A,sampIDs_BUFFER_B)
sampIDs_BUFFER <- sampIDs_BUFFER[which(sampIDs_BUFFER$Stain == "PROTEO"),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER, sampIDs_BUFFER[c(2,4),])
sampIDs_BUFFER$Experiment <- c("FT","4C","RT","FT","4C","RT")
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_BUFFER_exp <- rep(sampIDs_BUFFER$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well)

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER <- rbind(melt_file_BUFFER, melt_file_BUFFER[349:696,])
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep
melt_file_BUFFER$Experiment <- sampIDs_BUFFER_exp


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B20","G20")),]

#extract_numeric(sampIDs$Well)
sampIDs_A <- sampIDs[parse_number(sampIDs$Well) == 1,]
sampIDs_B <- sampIDs[parse_number(sampIDs$Well) == 2,]
sampIDs_sample <- rbind(sampIDs_A,sampIDs_B)
sampIDs_sample <- sampIDs_sample[which(sampIDs_sample$Stain == "PROTEO"),]
sampIDs_sample <- sampIDs_sample[order(sampIDs_sample$Well),] 
#parse_character(sampIDs$Well)
sampIDs_sample_rep <- rep(sampIDs_sample$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_sample_exp <- rep(sampIDs_sample$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_sample$Well)
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_sample <- melt_file[A,]
melt_file_sample$Sample <- sampIDs_sample_rep
melt_file_sample$Experiment <- sampIDs_sample_exp


melt_file_sample$Fluorescence <- gsub(",","",melt_file_sample$Fluorescence)
melt_file_sample$Fluorescence <- as.numeric(as.character(melt_file_sample$Fluorescence))

melt_file_sample <- melt_file_sample %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_sample_plot <- melt_file_sample[which(melt_file_sample$Well.Position %in% c("B1","G1","K1")),]

#subtract background
melt_file_sample_plot$Fluorescence <- melt_file_sample_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_sample_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_sample_plot$percAGG <- percAGG


#melt_file_K_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_sample_plot$Temperature))
max_temp <- round(max(melt_file_sample_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_sample_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Experiment),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Experiment)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 1.9E11 | PROTEOSTAT | 0% PS20")+
  labs(color="Condition")

plot_melt

##############################DP 1.9E11 0.001%####################################
sampIDs_BUFFER_A <- sampIDs[grepl('21',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_B <- sampIDs[grepl('22',sampIDs$Well, fixed=T),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER_A,sampIDs_BUFFER_B)
sampIDs_BUFFER <- sampIDs_BUFFER[which(sampIDs_BUFFER$Stain == "PROTEO"),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER, sampIDs_BUFFER[c(2,4),])
sampIDs_BUFFER$Experiment <- c("FT","4C","RT","FT","4C","RT")
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_BUFFER_exp <- rep(sampIDs_BUFFER$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well)

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER <- rbind(melt_file_BUFFER, melt_file_BUFFER[349:696,])
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep
melt_file_BUFFER$Experiment <- sampIDs_BUFFER_exp


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B21","G21")),]

#extract_numeric(sampIDs$Well)
sampIDs_A <- sampIDs[parse_number(sampIDs$Well) == 3,]
sampIDs_B <- sampIDs[parse_number(sampIDs$Well) == 4,]
sampIDs_sample <- rbind(sampIDs_A,sampIDs_B)
sampIDs_sample <- sampIDs_sample[which(sampIDs_sample$Stain == "PROTEO"),]
sampIDs_sample <- sampIDs_sample[order(sampIDs_sample$Well),] 
#parse_character(sampIDs$Well)
sampIDs_sample_rep <- rep(sampIDs_sample$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_sample_exp <- rep(sampIDs_sample$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_sample$Well)
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_sample <- melt_file[A,]
melt_file_sample$Sample <- sampIDs_sample_rep
melt_file_sample$Experiment <- sampIDs_sample_exp


melt_file_sample$Fluorescence <- gsub(",","",melt_file_sample$Fluorescence)
melt_file_sample$Fluorescence <- as.numeric(as.character(melt_file_sample$Fluorescence))

melt_file_sample <- melt_file_sample %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_sample_plot <- melt_file_sample[which(melt_file_sample$Well.Position %in% c("B3","G3","K3")),]

#subtract background
melt_file_sample_plot$Fluorescence <- melt_file_sample_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_sample_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_sample_plot$percAGG <- percAGG


#melt_file_K_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_sample_plot$Temperature))
max_temp <- round(max(melt_file_sample_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_sample_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Experiment),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Experiment)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 1.9E11 | PROTEOSTAT | 0.001% PS20")+
  labs(color="Condition")

plot_melt
##############################DP 1.9E11 0.008%####################################
sampIDs_BUFFER_A <- sampIDs[grepl('23',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_B <- sampIDs[grepl('24',sampIDs$Well, fixed=T),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER_A,sampIDs_BUFFER_B)
sampIDs_BUFFER <- sampIDs_BUFFER[which(sampIDs_BUFFER$Stain == "PROTEO"),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER, sampIDs_BUFFER[c(2,4),])
sampIDs_BUFFER$Experiment <- c("FT","4C","RT","FT","4C","RT")
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_BUFFER_exp <- rep(sampIDs_BUFFER$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well)

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER <- rbind(melt_file_BUFFER, melt_file_BUFFER[349:696,])
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep
melt_file_BUFFER$Experiment <- sampIDs_BUFFER_exp


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B23","G23")),]

#extract_numeric(sampIDs$Well)
sampIDs_A <- sampIDs[parse_number(sampIDs$Well) == 5,]
sampIDs_B <- sampIDs[parse_number(sampIDs$Well) == 6,]
sampIDs_sample <- rbind(sampIDs_A,sampIDs_B)
sampIDs_sample <- sampIDs_sample[which(sampIDs_sample$Stain == "PROTEO"),]
sampIDs_sample <- sampIDs_sample[order(sampIDs_sample$Well),] 
#parse_character(sampIDs$Well)
sampIDs_sample_rep <- rep(sampIDs_sample$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_sample_exp <- rep(sampIDs_sample$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_sample$Well)
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_sample <- melt_file[A,]
melt_file_sample$Sample <- sampIDs_sample_rep
melt_file_sample$Experiment <- sampIDs_sample_exp


melt_file_sample$Fluorescence <- gsub(",","",melt_file_sample$Fluorescence)
melt_file_sample$Fluorescence <- as.numeric(as.character(melt_file_sample$Fluorescence))

melt_file_sample <- melt_file_sample %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_sample_plot <- melt_file_sample[which(melt_file_sample$Well.Position %in% c("B5","G5","K5")),]

#subtract background
melt_file_sample_plot$Fluorescence <- melt_file_sample_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_sample_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_sample_plot$percAGG <- percAGG


#melt_file_K_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_sample_plot$Temperature))
max_temp <- round(max(melt_file_sample_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_sample_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Experiment),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Experiment)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 1.9E11 | PROTEOSTAT | 0.008% PS20")+
  labs(color="Condition")

plot_melt





###### AGGREGATE PERCENTAGE ACROSS CONDS FOR 3.5E11
##############################DP 3.5E11 0%#######################################
sampIDs_BUFFER_A <- sampIDs[grepl('19',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_B <- sampIDs[grepl('20',sampIDs$Well, fixed=T),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER_A,sampIDs_BUFFER_B)
sampIDs_BUFFER <- sampIDs_BUFFER[which(sampIDs_BUFFER$Stain == "PROTEO"),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER, sampIDs_BUFFER[c(2,4),])
sampIDs_BUFFER$Experiment <- c("FT","4C","RT","FT","4C","RT")
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_BUFFER_exp <- rep(sampIDs_BUFFER$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well)

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER <- rbind(melt_file_BUFFER, melt_file_BUFFER[349:696,])
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep
melt_file_BUFFER$Experiment <- sampIDs_BUFFER_exp


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B20","G20")),]

#extract_numeric(sampIDs$Well)
sampIDs_A <- sampIDs[parse_number(sampIDs$Well) == 7,]
sampIDs_B <- sampIDs[parse_number(sampIDs$Well) == 8,]
sampIDs_sample <- rbind(sampIDs_A,sampIDs_B)
sampIDs_sample <- sampIDs_sample[which(sampIDs_sample$Stain == "PROTEO"),]
sampIDs_sample <- sampIDs_sample[order(sampIDs_sample$Well),] 
#parse_character(sampIDs$Well)
sampIDs_sample_rep <- rep(sampIDs_sample$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_sample_exp <- rep(sampIDs_sample$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_sample$Well)
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_sample <- melt_file[A,]
melt_file_sample$Sample <- sampIDs_sample_rep
melt_file_sample$Experiment <- sampIDs_sample_exp


melt_file_sample$Fluorescence <- gsub(",","",melt_file_sample$Fluorescence)
melt_file_sample$Fluorescence <- as.numeric(as.character(melt_file_sample$Fluorescence))

melt_file_sample <- melt_file_sample %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_sample_plot <- melt_file_sample[which(melt_file_sample$Well.Position %in% c("B7","G7","K7")),]

#subtract background
melt_file_sample_plot$Fluorescence <- melt_file_sample_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_sample_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_sample_plot$percAGG <- percAGG


#melt_file_K_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_sample_plot$Temperature))
max_temp <- round(max(melt_file_sample_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_sample_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Experiment),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Experiment)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 3.5E11 | PROTEOSTAT | 0% PS20")+
  labs(color="Condition")

plot_melt


##############################DP 3.5E11 0.001%#######################################
sampIDs_BUFFER_A <- sampIDs[grepl('21',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_B <- sampIDs[grepl('22',sampIDs$Well, fixed=T),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER_A,sampIDs_BUFFER_B)
sampIDs_BUFFER <- sampIDs_BUFFER[which(sampIDs_BUFFER$Stain == "PROTEO"),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER, sampIDs_BUFFER[c(2,4),])
sampIDs_BUFFER$Experiment <- c("FT","4C","RT","FT","4C","RT")
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_BUFFER_exp <- rep(sampIDs_BUFFER$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well)

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER <- rbind(melt_file_BUFFER, melt_file_BUFFER[349:696,])
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep
melt_file_BUFFER$Experiment <- sampIDs_BUFFER_exp


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B21","G21")),]

#extract_numeric(sampIDs$Well)
sampIDs_A <- sampIDs[parse_number(sampIDs$Well) == 9,]
sampIDs_B <- sampIDs[parse_number(sampIDs$Well) == 10,]
sampIDs_sample <- rbind(sampIDs_A,sampIDs_B)
sampIDs_sample <- sampIDs_sample[which(sampIDs_sample$Stain == "PROTEO"),]
sampIDs_sample <- sampIDs_sample[order(sampIDs_sample$Well),] 
#parse_character(sampIDs$Well)
sampIDs_sample_rep <- rep(sampIDs_sample$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_sample_exp <- rep(sampIDs_sample$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_sample$Well)
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_sample <- melt_file[A,]
melt_file_sample$Sample <- sampIDs_sample_rep
melt_file_sample$Experiment <- sampIDs_sample_exp


melt_file_sample$Fluorescence <- gsub(",","",melt_file_sample$Fluorescence)
melt_file_sample$Fluorescence <- as.numeric(as.character(melt_file_sample$Fluorescence))

melt_file_sample <- melt_file_sample %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_sample_plot <- melt_file_sample[which(melt_file_sample$Well.Position %in% c("B9","G9","K9")),]

#subtract background
melt_file_sample_plot$Fluorescence <- melt_file_sample_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_sample_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_sample_plot$percAGG <- percAGG


#melt_file_K_plot$percAGG <- melt_file_G_plot$Fluorescence / melt_file_AGG$Fluorescence *100


# AGG <- data.frame(Sample = melt_file_G_plot$Sample,
#                   Temperature = melt_file_G_plot$Temperature, 
#                   percAGG = melt_file_G_plot$percAGG,
#                   Well.Position = melt_file_G_plot$Well.Position,
#                   Target.Name = melt_file_G_plot$Target.Name)

min_temp <- round(min(melt_file_sample_plot$Temperature))
max_temp <- round(max(melt_file_sample_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_sample_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Experiment),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Experiment)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 3.5E11 | PROTEOSTAT | 0.001% PS20")+
  labs(color="Condition")

plot_melt
##############################DP 3.5E11 0.008%#######################################
sampIDs_BUFFER_A <- sampIDs[grepl('23',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_B <- sampIDs[grepl('24',sampIDs$Well, fixed=T),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER_A,sampIDs_BUFFER_B)
sampIDs_BUFFER <- sampIDs_BUFFER[which(sampIDs_BUFFER$Stain == "PROTEO"),]
sampIDs_BUFFER <- rbind(sampIDs_BUFFER, sampIDs_BUFFER[c(2,4),])
sampIDs_BUFFER$Experiment <- c("FT","4C","RT","FT","4C","RT")
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_BUFFER_exp <- rep(sampIDs_BUFFER$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well)

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER <- rbind(melt_file_BUFFER, melt_file_BUFFER[349:696,])
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep
melt_file_BUFFER$Experiment <- sampIDs_BUFFER_exp


melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B23","G23")),]

#extract_numeric(sampIDs$Well)
sampIDs_A <- sampIDs[parse_number(sampIDs$Well) == 11,]
sampIDs_B <- sampIDs[parse_number(sampIDs$Well) == 12,]
sampIDs_sample <- rbind(sampIDs_A,sampIDs_B)
sampIDs_sample <- sampIDs_sample[which(sampIDs_sample$Stain == "PROTEO"),]
sampIDs_sample <- sampIDs_sample[order(sampIDs_sample$Well),] 
#parse_character(sampIDs$Well)
sampIDs_sample_rep <- rep(sampIDs_sample$Sample, each=174) #number of reads per well (in melt_file)
sampIDs_sample_exp <- rep(sampIDs_sample$Experiment, each=174) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_sample$Well)
# A = which(inx %in% c("G1","G2","G3","G4","G5","G6"))
# A = which(inx %in% c("G7","G8","G9","G10","G11","G12"))

melt_file_sample <- melt_file[A,]
melt_file_sample$Sample <- sampIDs_sample_rep
melt_file_sample$Experiment <- sampIDs_sample_exp


melt_file_sample$Fluorescence <- gsub(",","",melt_file_sample$Fluorescence)
melt_file_sample$Fluorescence <- as.numeric(as.character(melt_file_sample$Fluorescence))

melt_file_sample <- melt_file_sample %>%
  group_by(Experiment, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_sample_plot <- melt_file_sample[which(melt_file_sample$Well.Position %in% c("B11","G11","K11")),]

#subtract background
melt_file_sample_plot$Fluorescence <- melt_file_sample_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence

K = which(inx %in% c("K18"))
melt_file_AGG <- melt_file[rep(K,3),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:174,4],3) #number of samples
standardCurve_intercept <- rep(standardCurve[1:174,5],3)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_sample_plot$Fluorescence[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_sample_plot$percAGG <- percAGG


min_temp <- round(min(melt_file_sample_plot$Temperature))
max_temp <- round(max(melt_file_sample_plot$Temperature))
temps <- seq(min_temp,max_temp,by=1)


plot_melt <- ggplot(melt_file_sample_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Experiment),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG, color=Experiment)) + 
  scale_x_continuous(breaks=temps,
                     labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP 3.5E11 | PROTEOSTAT | 0.008% PS20")+
  labs(color="Condition")

plot_melt