#230-280 data
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

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

fc_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/FC_Counts.csv")
#fc_counts = fc_file[which(fc_file$Sample != "FFB2+PS20"),1:4]
fc_counts <- fc_file
fc_counts$Counts[which(fc_counts$Sample == "DI2")] <- fc_counts$Counts[which(fc_counts$Sample == "DI2")]*16 #adjust for dilution factor

plot_fc <- ggplot(fc_counts, aes(x=Temp, y=Counts, group=Sample,color=Sample))
plot_fc <- plot_fc + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc <- plot_fc + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Particulates/mL") + xlab("Temperature") + labs(title = "Particle Counts at Increasing Temperatures, by Sample",
                                                    subtitle = "19FEB2024")
plot_fc

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

turb_jmp = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity_JMP.csv")
#fc_counts = fc_file[which(fc_file$Sample != "FFB2+PS20"),1:4]

plot_fc <- ggplot(fc_counts, aes(x=Temp, y=Counts, group=Sample,color=Sample))
plot_fc <- plot_fc + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc <- plot_fc + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Particles/mL") + xlab("Temperature") + labs(title = "Particle Counts at Increasing Temperatures, by Sample",
                                                    subtitle = "15FEB2024")
plot_fc







turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity_19FEB2024.csv")

turb_230 = turb[which(turb$Wavelength == "230"),setdiff(1:5,3)]
turb_230 = turb_230[which(turb_230$Dilution == "1"),-which(colnames(turb_230) == "Dilution")]
turb_230_melt <- melt(turb_230,id= c("Sample","Temp"))

turb_230_melt_samp <- turb_230_melt[grepl('DP',turb_230_melt$Sample, fixed=T),]

plot_230 <- ggplot(turb_230_melt_samp, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_230 <- plot_230 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_230 <- plot_230 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 230nm") + xlab("Temperature")+ labs(title = "OD Values at Increasing Temperatures, by Sample",
  subtitle = "230nm")
plot_230

turb_280 = turb[which(turb$Wavelength == "280"),setdiff(1:5,3)]
turb_280 = turb_280[which(turb_280$Dilution == "1"),-which(colnames(turb_280) == "Dilution")]
turb_280_melt <- melt(turb_280,id= c("Sample","Temp"))

turb_280_melt_samp <- turb_280_melt[grepl('DP',turb_280_melt$Sample, fixed=T),]

plot_280 <- ggplot(turb_280_melt_samp, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_280 <- plot_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_280 <- plot_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 280nm") + xlab("Temperature")+ labs(title = "OD Values at Increasing Temperatures, by Sample",
  subtitle = "280nm")
plot_280


turb_ratio <- cbind(turb_230_melt, data.frame(x280 = turb_280_melt$value), data.frame(ratio=turb_230_melt$value/turb_280_melt$value))
#filter
turb_ratio_samp <- turb_ratio[grepl('DP',turb_ratio$Sample, fixed=T),]

plot_230_280 <- ggplot(turb_ratio_samp, aes(x=Temp, y=ratio, group=Sample,color=Sample))
plot_230_280 <- plot_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_230_280 <- plot_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Temperature")+ labs(title = "OD Values at Increasing Temperatures, by Sample",
 subtitle = "Aggregate Percentage, 230/280")
plot_230_280


turb_650 = turb[which(turb$Wavelength == "650"),setdiff(1:5,3)]
turb_650 = turb_650[which(turb_650$Dilution == "1"),-which(colnames(turb_650) == "Dilution")]
turb_650_melt <- melt(turb_650,id= c("Sample","Temp"))

turb_650_melt_samp <- turb_650_melt[grepl('DP',turb_650_melt$Sample, fixed=T),]

plot_650 <- ggplot(turb_650_melt_samp, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_650 <- plot_650 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_650 <- plot_650 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 650nm") + xlab("Temperature")+ labs(title = "OD Values at Increasing Temperatures, by Sample",
  subtitle = "650nm")
plot_650



turb_980 = turb[which(turb$Wavelength == "980"),setdiff(1:5,3)]
turb_980 = turb_980[which(turb_980$Dilution == "1"),-which(colnames(turb_980) == "Dilution")]
turb_980_melt <- melt(turb_980,id= c("Sample","Temp"))

turb_980_melt_samp <- turb_980_melt[grepl('DP',turb_980_melt$Sample, fixed=T),]

plot_980 <- ggplot(turb_980_melt_samp, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_980 <- plot_980 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_980 <- plot_980 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 980nm") + xlab("Temperature")+ labs(title = "OD Values at Increasing Temperatures, by Sample",
  subtitle = "980nm")
plot_980

























temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_280_full = turb[which(turb$Wavelength == "280"),setdiff(1:5,3)]


turb_280_full_DP0 <- turb_280_full[which(turb_280_full$Sample == "DP 0%"),]
plot_280_full <- ggplot(turb_280_full_DP0, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_280_full$labels$colour <- "Temp(C)"
plot_280_full