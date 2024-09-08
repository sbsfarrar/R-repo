library(ggplot2)
#theme_set(theme_grey())
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
library(cowplot)
library(magick)

##############################DP STEP 1 - MELT CURVE############################################
melt_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 1_SS19APR2024.csv")
sampIDs = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 1_SS19APR2024_SampleIds.csv")
standardCurve = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 1_SS19APR2024_StandardCurve.csv")

sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[which(sampIDs_BUFFER$Sample == "FFB ")], each=211) #number of reads per well (in melt_file)

buffNum <- which(sampIDs_BUFFER$Sample == "FFB ")
inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[buffNum])

melt_file_BUFFER <- melt_file[A,]
melt_file_BUFFER$Sample <- sampIDs_BUFFER_rep

melt_file_BUFFER$Fluorescence <- gsub(",","",melt_file_BUFFER$Fluorescence)
melt_file_BUFFER$Fluorescence <- as.numeric(as.character(melt_file_BUFFER$Fluorescence))

melt_file_BUFFER <- melt_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_BUFFER_background <- melt_file_BUFFER[which(melt_file_BUFFER$Well.Position %in% c("B9")),]



sampIDs_B <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_B_rep <- rep(sampIDs_B$Sample[1:8], each=211) #number of reads per well (in melt_file)

inx <- as.matrix(melt_file["Well.Position"])
A = which(inx %in% sampIDs_B$Well[1:8])

melt_file_B <- melt_file[A,]
melt_file_B$Sample <- sampIDs_B_rep

melt_file_B$Fluorescence <- gsub(",","",melt_file_B$Fluorescence)
melt_file_B$Fluorescence <- as.numeric(as.character(melt_file_B$Fluorescence))

melt_file_B <- melt_file_B %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt_file_B_plot <- melt_file_B[which(melt_file_B$Well.Position %in% c("B2","B4","B6","B8")),]

#subtract background
melt_file_B_plot$FluorescenceSub <- melt_file_B_plot$Fluorescence - melt_file_BUFFER_background$Fluorescence #check

B = which(inx %in% c("B11","B12"))
melt_file_AGG <- melt_file[rep(B,4),]

# melt_file_AGG <- melt_file[which(melt_file$Well.Position == c("K18")),]

melt_file_AGG$Fluorescence <- gsub(",","",melt_file_AGG$Fluorescence)
melt_file_AGG$Fluorescence <- as.numeric(as.character(melt_file_AGG$Fluorescence))

melt_file_AGG <- melt_file_AGG %>%
  group_by(Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:211,4],4) #number of samples
standardCurve_intercept <- rep(standardCurve[1:211,5],4)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt_file_B_plot$FluorescenceSub[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt_file_B_plot$percAGG <- percAGG


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
                                                        subtitle = "DP17APR2024 | PROTEOSTAT | STEP 1")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_melt

##############################DP STEP 2 - 25C CYCLE HOLD#############################################
hold_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 2_SS19APR2024.csv")
#same sample IDs as before, duh
sampIDs = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 1_SS19APR2024_SampleIds.csv")
standardCurve = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 2_SS19APR2024_StandardCurve.csv")
sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[which(sampIDs_BUFFER$Sample == "FFB ")], each=120) #number of reads per well (in hold_file), can find max in read column

buffNum <- which(sampIDs_BUFFER$Sample == "FFB ")
inx <- as.matrix(hold_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[buffNum])

hold_file_BUFFER <- hold_file[A,]
hold_file_BUFFER$Sample <- sampIDs_BUFFER_rep

hold_file_BUFFER$Fluorescence <- gsub(",","",hold_file_BUFFER$Fluorescence)
hold_file_BUFFER$Fluorescence <- as.numeric(as.character(hold_file_BUFFER$Fluorescence))

hold_file_BUFFER <- hold_file_BUFFER %>%
  group_by(Sample, Cycle) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

hold_file_BUFFER_background <- hold_file_BUFFER[which(hold_file_BUFFER$Well.Position %in% c("B9")),]



sampIDs_B <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_B_rep <- rep(sampIDs_B$Sample[1:8], each=120) #number of reads per well (in hold_file)

inx <- as.matrix(hold_file["Well.Position"])
A = which(inx %in% sampIDs_B$Well[1:8])

hold_file_B <- hold_file[A,]
hold_file_B$Sample <- sampIDs_B_rep

hold_file_B$Fluorescence <- gsub(",","",hold_file_B$Fluorescence)
hold_file_B$Fluorescence <- as.numeric(as.character(hold_file_B$Fluorescence))

hold_file_B <- hold_file_B %>%
  group_by(Sample, Cycle) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

hold_file_B_plot <- hold_file_B[which(hold_file_B$Well.Position %in% c("B2","B4","B6","B8")),]

#subtract background
hold_file_B_plot$FluorescenceSub <- hold_file_B_plot$Fluorescence - hold_file_BUFFER_background$Fluorescence #check

B = which(inx %in% c("B11","B12"))
hold_file_AGG <- hold_file[rep(B,4),]

# hold_file_AGG <- hold_file[which(hold_file$Well.Position == c("K18")),]

hold_file_AGG$Fluorescence <- gsub(",","",hold_file_AGG$Fluorescence)
hold_file_AGG$Fluorescence <- as.numeric(as.character(hold_file_AGG$Fluorescence))

hold_file_AGG <- hold_file_AGG %>%
  group_by(Cycle) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:120,4],4) #number of samples
standardCurve_intercept <- rep(standardCurve[1:120,5],4)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(hold_file_B_plot$FluorescenceSub[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
hold_file_B_plot$percAGG <- percAGG


# min_temp <- round(min(hold_file_B_plot$Temperature))
# max_temp <- round(max(hold_file_B_plot$Temperature))
# temps <- seq(min_temp,max_temp,by=1)
cycs = seq(1,120,by=1)


plot_hold <- ggplot(hold_file_B_plot) + 
  geom_line(aes(x = Cycle, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Cycle, y=percAGG, color=Sample)) + 
  scale_x_continuous(breaks=cycs[c(FALSE,FALSE,FALSE,FALSE,TRUE)],
                     labels= as.character(cycs[c(FALSE,FALSE,FALSE,FALSE,TRUE)])) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Cycle")+ labs(title = "Percent Aggregate per Cycle at 25C",
                                                        subtitle = "DP17APR2024 | PROTEOSTAT | STEP 2")+
  labs(color="Sample")
# scale_color_manual(labels=c("0% PS20","0.001% PS20","0.008% PS20"), 
#                    values=c("black","blue","green", "yellow", "orange", "red"))
plot_hold

##############################DP STEP 3 - MELT CURVE############################################
melt2_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 3_SS19APR2024.csv")
sampIDs = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 1_SS19APR2024_SampleIds.csv")
standardCurve = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/2024-04-19-SS/Aggregate Detection Step 3_SS19APR2024_StandardCurve.csv")

sampIDs_BUFFER <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_BUFFER_rep <- rep(sampIDs_BUFFER$Sample[which(sampIDs_BUFFER$Sample == "FFB ")], each=196) #number of reads per well (in melt2_file)

buffNum <- which(sampIDs_BUFFER$Sample == "FFB ")
inx <- as.matrix(melt2_file["Well.Position"])
A = which(inx %in% sampIDs_BUFFER$Well[buffNum])

melt2_file_BUFFER <- melt2_file[A,]
melt2_file_BUFFER$Sample <- sampIDs_BUFFER_rep

melt2_file_BUFFER$Fluorescence <- gsub(",","",melt2_file_BUFFER$Fluorescence)
melt2_file_BUFFER$Fluorescence <- as.numeric(as.character(melt2_file_BUFFER$Fluorescence))

melt2_file_BUFFER <- melt2_file_BUFFER %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt2_file_BUFFER_background <- melt2_file_BUFFER[which(melt2_file_BUFFER$Well.Position %in% c("B9")),]



sampIDs_B <- sampIDs[grepl('B',sampIDs$Well, fixed=T),]
sampIDs_B_rep <- rep(sampIDs_B$Sample[1:8], each=196) #number of reads per well (in melt2_file)

inx <- as.matrix(melt2_file["Well.Position"])
A = which(inx %in% sampIDs_B$Well[1:8])

melt2_file_B <- melt2_file[A,]
melt2_file_B$Sample <- sampIDs_B_rep

melt2_file_B$Fluorescence <- gsub(",","",melt2_file_B$Fluorescence)
melt2_file_B$Fluorescence <- as.numeric(as.character(melt2_file_B$Fluorescence))

melt2_file_B <- melt2_file_B %>%
  group_by(Sample, Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

melt2_file_B_plot <- melt2_file_B[which(melt2_file_B$Well.Position %in% c("B2","B4","B6","B8")),]

#subtract background
melt2_file_B_plot$FluorescenceSub <- melt2_file_B_plot$Fluorescence - melt2_file_BUFFER_background$Fluorescence #check

B = which(inx %in% c("B11","B12"))
melt2_file_AGG <- melt2_file[rep(B,4),]

# melt2_file_AGG <- melt2_file[which(melt2_file$Well.Position == c("K18")),]

melt2_file_AGG$Fluorescence <- gsub(",","",melt2_file_AGG$Fluorescence)
melt2_file_AGG$Fluorescence <- as.numeric(as.character(melt2_file_AGG$Fluorescence))

melt2_file_AGG <- melt2_file_AGG %>%
  group_by(Reading) %>%
  mutate(Mean = mean(Fluorescence)) %>%
  ungroup()

interp_fluorescence <- function(x, slope, intercept){
  return ((x)*slope + intercept)
}


standardCurve_slope <- rep(standardCurve[1:196,4],4) #number of samples
standardCurve_intercept <- rep(standardCurve[1:196,5],4)
count = length(standardCurve_slope)
percAGG <- matrix(data=0, nrow=length(standardCurve_slope))
for (i in 1:count){
  percAGG[i,1] <- interp_fluorescence(melt2_file_B_plot$FluorescenceSub[i],
                                      standardCurve_slope[i],
                                      standardCurve_intercept[i])
}
melt2_file_B_plot$percAGG <- percAGG


min_temp <- round(min(melt2_file_B_plot$Temperature))
max_temp <- round(max(melt2_file_B_plot$Temperature))
temps <- seq(max_temp,min_temp,by=-1)


plot_melt2 <- ggplot(melt2_file_B_plot) + 
  geom_line(aes(x = Temperature, y=percAGG,color=Sample),size=0.8,alpha = 0.5) + 
  geom_point(aes(x= Temperature, y=percAGG,color=Sample)) + 
  scale_x_reverse(breaks=temps,
                  labels = as.character(temps))+
  # scale_x_continuous(breaks=temps,
  #                    labels= as.character(temps)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) + 
  ylab("Percent Aggregate") + xlab("Temperature")+ labs(title = "Percent Aggregate by Temperature",
                                                        subtitle = "DP17APR2024 | PROTEOSTAT | STEP 3")+
  labs(color="Sample")
plot_melt2

















##########################method image overlay###################

logo_file <- "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/qPCR/Aggregate Detection Method Step 2.png"
theme_set(theme_cowplot())
p1 <- ggdraw() +
  draw_image(logo_file, x=0.3,y=0.4,scale=0.5)+
  draw_plot(plot_melt) +
  theme_cowplot()
p1

image = image_read(logo_file)
p2 <- image_ggplot(image, interpolate = FALSE)
p2