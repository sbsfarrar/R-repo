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
#library(wesanderson)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

fc_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/FC_Counts_22FEB2024.csv")
#fc_counts = fc_file[which(fc_file$Sample != "FFB2+PS20"),1:4]
fc_counts <- fc_file[,1:6]
fc_counts$Total.Particulate[which(fc_counts$Sample == "DI2")] <- fc_counts$Total.Particulate[which(fc_counts$Sample == "DI2")]*10 #adjust for dilution factor
fc_counts$Small.Particulate[which(fc_counts$Sample == "DI2")] <- fc_counts$Small.Particulate[which(fc_counts$Sample == "DI2")]*10 #adjust for dilution factor
fc_counts$Large.Particulate[which(fc_counts$Sample == "DI2")] <- fc_counts$Large.Particulate[which(fc_counts$Sample == "DI2")]*10 #adjust for dilution factor
#remove spin
fc_counts <- fc_counts[which(fc_counts$Sample != "DP0 SPIN"),1:6]
plot_fc <- ggplot(fc_counts, aes(x=Temp, y=Total.Particulate, group=Sample,color=Sample))
plot_fc <- plot_fc + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc <- plot_fc + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Particles/mL") + xlab("Temperature") + labs(title = "Total Particulate Counts at Increasing Temperatures, by Sample",
  subtitle = "22FEB2024")
plot_fc

plot_fc_small <- ggplot(fc_counts, aes(x=Temp, y=Small.Particulate, group=Sample,color=Sample))
plot_fc_small <- plot_fc_small + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc_small <- plot_fc_small + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Particles/mL") + xlab("Temperature") + labs(title = "Small Particulate Counts at Increasing Temperatures, by Sample",
                                                    subtitle = "22FEB2024")
plot_fc_small

plot_fc_large <- ggplot(fc_counts, aes(x=Temp, y=Large.Particulate, group=Sample,color=Sample))
plot_fc_large <- plot_fc_large + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc_large <- plot_fc_large + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Particles/mL") + xlab("Temperature") + labs(title = "Large Particulate Counts at Increasing Temperatures, by Sample",
                                                    subtitle = "22FEB2024")
plot_fc_large

fc_counts$PercentSmall <- fc_counts$Small.Particulate/(fc_counts$Small.Particulate + 
                                                         fc_counts$Large.Particulate)
fc_counts$PercentLarge <- fc_counts$Large.Particulate/(fc_counts$Small.Particulate + 
                                                         fc_counts$Large.Particulate)

plot_fc_smallperc <- ggplot(fc_counts, aes(x=Temp, y=PercentSmall, group=Sample,color=Sample))
plot_fc_smallperc <- plot_fc_smallperc + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc_smallperc <- plot_fc_smallperc + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Percent of Particulates Counted") + xlab("Temperature") + labs(title = "Percent 'Small' Particulate at Increasing Temperatures, by Sample",
                                                    subtitle = "22FEB2024")
plot_fc_smallperc

plot_fc_largeperc <- ggplot(fc_counts, aes(x=Temp, y=PercentLarge, group=Sample,color=Sample))
plot_fc_largeperc <- plot_fc_largeperc + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_fc_largeperc <- plot_fc_largeperc + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("Percent of Particulates Counted") + xlab("Temperature") + labs(title = "Percent 'Large' Particulate at Increasing Temperatures, by Sample",
                                                    subtitle = "22FEB2024")
plot_fc_largeperc
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
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity_22FEB2024.csv")
turb["OD"][turb["OD"] < 0] <- 0
turb_DP <- turb[grepl('DP',turb$Sample, fixed=T),]
turb_DP <- turb_DP[which(turb_DP$Sample != "DP 0% SPIN"),]
turb_DI2 <-turb[grepl('DI2',turb$Sample, fixed=T),]
turb_buff <- turb[grepl('FFB2',turb$Sample, fixed=T),]
fc_counts_DP <- fc_counts[which(fc_counts$Sample != "DI2"),1:6]

turb_230 = turb_DP[which(turb_DP$Wavelength == "230"),setdiff(1:5,3)]
turb_230 = turb_230[which(turb_230$Dilution == "2"),-which(colnames(turb_230) == "Dilution")]
turb_230_melt <- melt(turb_230,id= c("Sample","Temp"))

plot_230 <- ggplot(turb_230_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_230 <- plot_230 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_230 <- plot_230 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 230nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
#  subtitle = "230nm")
plot_230
#from https://finchstudio.io/blog/ggplot-dual-y-axes/
max_first  <- max(mseResultsDF$MSE)   # Specify max of first y axis
max_second <- max(mseResultsDF$Accuracy) # Specify max of second y axis
min_first  <- min(mseResultsDF$MSE)   # Specify min of first y axis
min_second <- min(mseResultsDF$Accuracy) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x, scale, shift){
  return ((x)*scale - shift)
}

# Function to scale secondary variable values
inv_scale_function <- function(x, scale, shift){
  return ((x + shift)/scale)
}

p2 <- ggplot(mseResultsDF,aes(x=k)) + 
  geom_line(aes(y=MSE,alpha=0.5,color="MSE")) + geom_point(aes(y=MSE, color="MSE"),alpha=0.5) + 
  geom_line(aes(y=inv_scale_function(Accuracy,scale,shift),alpha=0.5,color="Accuracy")) + 
  geom_point(aes(y=inv_scale_function(Accuracy,scale,shift), color="Accuracy"),alpha=0.5) + 
  scale_color_brewer(palette="Set1") + guides(alpha="none") + 
  ggtitle("K Parameter Accuracy and Mean Squared Error") + 
  theme(plot.title=element_text(hjust = 0.5)) + xlab("C-parameter Value") +
  ylab("Model Accuracy") + labs(colour="") +
  scale_x_continuous(breaks=kVals[c(FALSE,TRUE,FALSE)],labels=xLabel[c(FALSE,TRUE,FALSE)]) +
  scale_y_continuous("Mean Squared Error",  # limits=c(min_first,max_first)
                     sec.axis = sec_axis(~scale_function(.,scale,shift),name="Model Accuracy")) +
  theme(axis.text.x = element_text(angle=90,hjust=0.5))
p2


turb_260 = turb_DP[which(turb_DP$Wavelength == "260"),setdiff(1:5,3)]
turb_260 = turb_260[which(turb_260$Dilution == "2"),-which(colnames(turb_260) == "Dilution")]
turb_260_melt <- melt(turb_260,id= c("Sample","Temp"))

plot_260 <- ggplot(turb_260_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_260 <- plot_260 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_260 <- plot_260 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 260nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
#  subtitle = "260nm")
plot_260


turb_280 = turb_DP[which(turb_DP$Wavelength == "280"),setdiff(1:5,3)]
turb_280 = turb_280[which(turb_280$Dilution == "2"),-which(colnames(turb_280) == "Dilution")]
turb_280_melt <- melt(turb_280,id= c("Sample","Temp"))

plot_280 <- ggplot(turb_280_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_280 <- plot_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_280 <- plot_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 280nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
#  subtitle = "280nm")
plot_280

turb_340 = turb_DP[which(turb_DP$Wavelength == "340"),setdiff(1:5,3)]
turb_340 = turb_340[which(turb_340$Dilution == "2"),-which(colnames(turb_340) == "Dilution")]
turb_340_melt <- melt(turb_340,id= c("Sample","Temp"))

plot_340 <- ggplot(turb_340_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_340 <- plot_340 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_340 <- plot_340 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 340nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "340nm")
plot_340

turb_650 = turb_DP[which(turb_DP$Wavelength == "650"),setdiff(1:5,3)]
turb_650 = turb_650[which(turb_650$Dilution == "2"),-which(colnames(turb_650) == "Dilution")]
turb_650_melt <- melt(turb_650,id= c("Sample","Temp"))

plot_650 <- ggplot(turb_650_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_650 <- plot_650 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_650 <- plot_650 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 650nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "650nm")
plot_650




turb_ratio <- cbind(turb_230_melt, data.frame(ratio=turb_230_melt$value/turb_280_melt$value))

plot_230_280 <- ggplot(turb_ratio, aes(x=Temp, y=ratio, group=Sample,color=Sample))
plot_230_280 <- plot_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_230_280 <- plot_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_230_280

grid_plot <- ggarrange(plot_230,plot_260,plot_280,plot_340, plot_650, plot_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at Increasing Temperatures, by Sample", 
                                           face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

turb_230 = turb_buff[which(turb_buff$Wavelength == "230"),setdiff(1:5,3)]
turb_230 = turb_230[which(turb_230$Dilution == "1"),-which(colnames(turb_230) == "Dilution")]
turb_230_melt <- melt(turb_230,id= c("Sample","Temp"))

plot_230 <- ggplot(turb_230_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_230 <- plot_230 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_230 <- plot_230 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 230nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
#  subtitle = "230nm")
plot_230


turb_260 = turb_buff[which(turb_buff$Wavelength == "260"),setdiff(1:5,3)]
turb_260 = turb_260[which(turb_260$Dilution == "1"),-which(colnames(turb_260) == "Dilution")]
turb_260_melt <- melt(turb_260,id= c("Sample","Temp"))

plot_260 <- ggplot(turb_260_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_260 <- plot_260 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_260 <- plot_260 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 260nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
#  subtitle = "260nm")
plot_260


turb_280 = turb_buff[which(turb_buff$Wavelength == "280"),setdiff(1:5,3)]
turb_280 = turb_280[which(turb_280$Dilution == "1"),-which(colnames(turb_280) == "Dilution")]
turb_280_melt <- melt(turb_280,id= c("Sample","Temp"))

plot_280 <- ggplot(turb_280_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_280 <- plot_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_280 <- plot_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 280nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
#  subtitle = "280nm")
plot_280

turb_340 = turb_buff[which(turb_buff$Wavelength == "340"),setdiff(1:5,3)]
turb_340 = turb_340[which(turb_340$Dilution == "1"),-which(colnames(turb_340) == "Dilution")]
turb_340_melt <- melt(turb_340,id= c("Sample","Temp"))

plot_340 <- ggplot(turb_340_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_340 <- plot_340 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_340 <- plot_340 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 340nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "340nm")
plot_340

turb_650 = turb_buff[which(turb_buff$Wavelength == "650"),setdiff(1:5,3)]
turb_650 = turb_650[which(turb_650$Dilution == "1"),-which(colnames(turb_650) == "Dilution")]
turb_650_melt <- melt(turb_650,id= c("Sample","Temp"))

plot_650 <- ggplot(turb_650_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_650 <- plot_650 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_650 <- plot_650 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 650nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "650nm")
plot_650




turb_ratio <- cbind(turb_230_melt, data.frame(ratio=turb_230_melt$value/turb_280_melt$value))

plot_230_280 <- ggplot(turb_ratio, aes(x=Temp, y=ratio, group=Sample,color=Sample))
plot_230_280 <- plot_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_230_280 <- plot_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_230_280

grid_plot <- ggarrange(plot_230,plot_260,plot_280,plot_340, plot_650, plot_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at Increasing Temperatures, by Sample", 
                                           face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#











turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity_3.csv")
turb_25 = turb[which(turb$Temp == "25"),setdiff(1:5,2)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_25_230 <- turb_25[which(turb_25$Wavelength == "230"),]
plot_25_1 <- ggplot(turb_25_230, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_25_1

turb_25_260 <- turb_25[which(turb_25$Wavelength == "260"),]
plot_25_2 <- ggplot(turb_25_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_25_2

turb_25_280 <- turb_25[which(turb_25$Wavelength == "280"),]
plot_25_3 <- ggplot(turb_25_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_25_3

turb_25_650 <- turb_25[which(turb_25$Wavelength == "650"),]
plot_25_4 <- ggplot(turb_25_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_25_4

turb_25_340 <- turb_25[which(turb_25$Wavelength == "340"),]
plot_25_5 <- ggplot(turb_25_340, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_25_5

turb_ratio_25 <- cbind(turb_25_230, data.frame(ratio=turb_25_230$OD/turb_25_280$OD))

plot_25_230_280 <- ggplot(turb_ratio_25, aes(x=1/Dilution, y=ratio, group=Sample,color=Sample))
plot_25_230_280 <- plot_25_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_25_230_280 <- plot_25_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_25_230_280

grid_plot <- ggarrange(plot_25_1, plot_25_2, plot_25_3, plot_25_4,plot_25_5,
                       plot_25_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at 25C across Selected Wavelengths, by Sample", 
                                           face = "bold", size = 14))  

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb_35 = turb[which(turb$Temp == "35"),setdiff(1:5,2)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_35_230 <- turb_35[which(turb_35$Wavelength == "230"),]
plot_35_1 <- ggplot(turb_35_230, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_35_1

turb_35_260 <- turb_35[which(turb_35$Wavelength == "260"),]
plot_35_2 <- ggplot(turb_35_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_35_2

turb_35_280 <- turb_35[which(turb_35$Wavelength == "280"),]
plot_35_3 <- ggplot(turb_35_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_35_3

turb_35_650 <- turb_35[which(turb_35$Wavelength == "650"),]
plot_35_4 <- ggplot(turb_35_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_35_4

turb_35_340 <- turb_35[which(turb_35$Wavelength == "340"),]
plot_35_5 <- ggplot(turb_35_340, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_35_5

turb_ratio_35 <- cbind(turb_35_230, data.frame(ratio=turb_35_230$OD/turb_35_280$OD))

plot_35_230_280 <- ggplot(turb_ratio_35, aes(x=1/Dilution, y=ratio, group=Sample,color=Sample))
plot_35_230_280 <- plot_35_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_35_230_280 <- plot_35_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_35_230_280

grid_plot <- ggarrange(plot_35_1, plot_35_2, plot_35_3, plot_35_4,plot_35_5,
                       plot_35_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at 35C across Selected Wavelengths, by Sample", 
                                           face = "bold", size = 14))  


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb_45 = turb[which(turb$Temp == "45"),setdiff(1:5,2)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_45_230 <- turb_45[which(turb_45$Wavelength == "230"),]
plot_45_1 <- ggplot(turb_45_230, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_45_1

turb_45_260 <- turb_45[which(turb_45$Wavelength == "260"),]
plot_45_2 <- ggplot(turb_45_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_45_2

turb_45_280 <- turb_45[which(turb_45$Wavelength == "280"),]
plot_45_3 <- ggplot(turb_45_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_45_3

turb_45_650 <- turb_45[which(turb_45$Wavelength == "650"),]
plot_45_4 <- ggplot(turb_45_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_45_4

turb_45_340 <- turb_45[which(turb_45$Wavelength == "340"),]
plot_45_5 <- ggplot(turb_45_340, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_45_5

turb_ratio_45 <- cbind(turb_45_230, data.frame(ratio=turb_45_230$OD/turb_45_280$OD))

plot_45_230_280 <- ggplot(turb_ratio_45, aes(x=1/Dilution, y=ratio, group=Sample,color=Sample))
plot_45_230_280 <- plot_45_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_45_230_280 <- plot_45_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_45_230_280

grid_plot <- ggarrange(plot_45_1, plot_45_2, plot_45_3, plot_45_4,plot_45_5,
                       plot_45_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at 45C across Selected Wavelengths, by Sample", 
                                           face = "bold", size = 14))  

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb_55 = turb[which(turb$Temp == "55"),setdiff(1:5,2)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_55_230 <- turb_55[which(turb_55$Wavelength == "230"),]
plot_55_1 <- ggplot(turb_55_230, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_55_1

turb_55_260 <- turb_55[which(turb_55$Wavelength == "260"),]
plot_55_2 <- ggplot(turb_55_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_55_2

turb_55_280 <- turb_55[which(turb_55$Wavelength == "280"),]
plot_55_3 <- ggplot(turb_55_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_55_3

turb_55_650 <- turb_55[which(turb_55$Wavelength == "650"),]
plot_55_4 <- ggplot(turb_55_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_55_4

turb_55_340 <- turb_55[which(turb_55$Wavelength == "340"),]
plot_55_5 <- ggplot(turb_55_340, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_55_5

turb_ratio_55 <- cbind(turb_55_230, data.frame(ratio=turb_55_230$OD/turb_55_280$OD))

plot_55_230_280 <- ggplot(turb_ratio_55, aes(x=1/Dilution, y=ratio, group=Sample,color=Sample))
plot_55_230_280 <- plot_55_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_55_230_280 <- plot_55_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_55_230_280

grid_plot <- ggarrange(plot_55_1, plot_55_2, plot_55_3, plot_55_4,plot_55_5,
                       plot_55_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at 55C across Selected Wavelengths, by Sample", 
                                           face = "bold", size = 14)) 

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb_65 = turb[which(turb$Temp == "65"),setdiff(1:5,2)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_65_230 <- turb_65[which(turb_65$Wavelength == "230"),]
plot_65_1 <- ggplot(turb_65_230, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_65_1

turb_65_260 <- turb_65[which(turb_65$Wavelength == "260"),]
plot_65_2 <- ggplot(turb_65_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_65_2

turb_65_280 <- turb_65[which(turb_65$Wavelength == "280"),]
plot_65_3 <- ggplot(turb_65_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_65_3

turb_65_650 <- turb_65[which(turb_65$Wavelength == "650"),]
plot_65_4 <- ggplot(turb_65_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_65_4

turb_65_340 <- turb_65[which(turb_65$Wavelength == "340"),]
plot_65_5 <- ggplot(turb_65_340, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_65_5

turb_ratio_65 <- cbind(turb_65_230, data.frame(ratio=turb_65_230$OD/turb_65_280$OD))

plot_65_230_280 <- ggplot(turb_ratio_65, aes(x=1/Dilution, y=ratio, group=Sample,color=Sample))
plot_65_230_280 <- plot_65_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_65_230_280 <- plot_65_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_65_230_280

grid_plot <- ggarrange(plot_65_1, plot_65_2, plot_65_3, plot_65_4,plot_65_5,
                       plot_65_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at 65C across Selected Wavelengths, by Sample", 
                                           face = "bold", size = 14)) 


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity_3.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_H2O = turb[which(turb$Sample == "H2O"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_H2O_230 <- turb_H2O[which(turb_H2O$Wavelength == "230"),]
plot_H2O_0 <- ggplot(turb_H2O_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_H2O_0$labels$colour <- "Temp(C)"
plot_H2O_0

turb_H2O_260 <- turb_H2O[which(turb_H2O$Wavelength == "260"),]
plot_H2O_1 <- ggplot(turb_H2O_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_H2O_1$labels$colour <- "Temp(C)"
plot_H2O_1

turb_H2O_280 <- turb_H2O[which(turb_H2O$Wavelength == "280"),]
plot_H2O_2 <- ggplot(turb_H2O_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_H2O_2$labels$colour <- "Temp(C)"
plot_H2O_2

turb_H2O_650 <- turb_H2O[which(turb_H2O$Wavelength == "650"),]
plot_H2O_3 <- ggplot(turb_H2O_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_H2O_3$labels$colour <- "Temp(C)"
plot_H2O_3

turb_H2O_340 <- turb_H2O[which(turb_H2O$Wavelength == "340"),]
plot_H2O_4 <- ggplot(turb_H2O_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_H2O_4$labels$colour <- "Temp(C)"
plot_H2O_4

turb_ratio_H2O <- cbind(turb_H2O_230, data.frame(ratio=turb_H2O_230$OD/turb_H2O_280$OD))

plot_H2O_230_280 <- ggplot(turb_ratio_H2O, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_H2O_230_280 <- plot_H2O_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_H2O_230_280 <- plot_H2O_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_H2O_230_280$labels$colour <- "Temp(C)"
plot_H2O_230_280

grid_plot <- ggarrange(plot_H2O_0, plot_H2O_1, plot_H2O_2, plot_H2O_3, plot_H2O_4, 
                       plot_H2O_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for H2O across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_FFB2 = turb[which(turb$Sample == "FFB2"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_FFB2_230 <- turb_FFB2[which(turb_FFB2$Wavelength == "230"),]
plot_FFB2_0 <- ggplot(turb_FFB2_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_0$labels$colour <- "Temp(C)"
plot_FFB2_0

turb_FFB2_260 <- turb_FFB2[which(turb_FFB2$Wavelength == "260"),]
plot_FFB2_1 <- ggplot(turb_FFB2_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_1$labels$colour <- "Temp(C)"
plot_FFB2_1

turb_FFB2_280 <- turb_FFB2[which(turb_FFB2$Wavelength == "280"),]
plot_FFB2_2 <- ggplot(turb_FFB2_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_2$labels$colour <- "Temp(C)"
plot_FFB2_2

turb_FFB2_650 <- turb_FFB2[which(turb_FFB2$Wavelength == "650"),]
plot_FFB2_3 <- ggplot(turb_FFB2_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_3$labels$colour <- "Temp(C)"
plot_FFB2_3

turb_FFB2_340 <- turb_FFB2[which(turb_FFB2$Wavelength == "340"),]
plot_FFB2_4 <- ggplot(turb_FFB2_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_4$labels$colour <- "Temp(C)"
plot_FFB2_4

turb_ratio_FFB2 <- cbind(turb_FFB2_230, data.frame(ratio=turb_FFB2_230$OD/turb_FFB2_280$OD))

plot_FFB2_230_280 <- ggplot(turb_ratio_FFB2, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_FFB2_230_280 <- plot_FFB2_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_FFB2_230_280 <- plot_FFB2_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_FFB2_230_280$labels$colour <- "Temp(C)"
plot_FFB2_230_280

grid_plot <- ggarrange(plot_FFB2_0, plot_FFB2_1, plot_FFB2_2, plot_FFB2_3, plot_FFB2_4, 
                       plot_FFB2_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for FFB2 across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_FFB2_PS20 = turb[which(turb$Sample == "FFB2+PS20"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_FFB2_PS20_230 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "230"),]
plot_FFB2_PS20_0 <- ggplot(turb_FFB2_PS20_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_PS20_0$labels$colour <- "Temp(C)"
plot_FFB2_PS20_0

turb_FFB2_PS20_260 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "260"),]
plot_FFB2_PS20_1 <- ggplot(turb_FFB2_PS20_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_PS20_1$labels$colour <- "Temp(C)"
plot_FFB2_PS20_1

turb_FFB2_PS20_280 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "280"),]
plot_FFB2_PS20_2 <- ggplot(turb_FFB2_PS20_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_PS20_2$labels$colour <- "Temp(C)"
plot_FFB2_PS20_2

turb_FFB2_PS20_650 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "650"),]
plot_FFB2_PS20_3 <- ggplot(turb_FFB2_PS20_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_PS20_3$labels$colour <- "Temp(C)"
plot_FFB2_PS20_3

turb_FFB2_PS20_340 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "340"),]
plot_FFB2_PS20_4 <- ggplot(turb_FFB2_PS20_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_PS20_4$labels$colour <- "Temp(C)"
plot_FFB2_PS20_4

turb_ratio_FFB2_PS20 <- cbind(turb_FFB2_PS20_230, data.frame(ratio=turb_FFB2_PS20_230$OD/turb_FFB2_PS20_280$OD))

plot_FFB2_PS20_230_280 <- ggplot(turb_ratio_FFB2_PS20, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_FFB2_PS20_230_280 <- plot_FFB2_PS20_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_FFB2_PS20_230_280 <- plot_FFB2_PS20_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_FFB2_PS20_230_280$labels$colour <- "Temp(C)"
plot_FFB2_PS20_230_280

grid_plot <- ggarrange(plot_FFB2_PS20_0, plot_FFB2_PS20_1, plot_FFB2_PS20_2, plot_FFB2_PS20_3, plot_FFB2_PS20_4, 
                       plot_FFB2_PS20_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for FFB2+PS20 across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_DP_0 = turb[which(turb$Sample == "DP 0%"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_DP_0_230 <- turb_DP_0[which(turb_DP_0$Wavelength == "230"),]
plot_DP_0_0 <- ggplot(turb_DP_0_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_0_0$labels$colour <- "Temp(C)"
plot_DP_0_0

turb_DP_0_260 <- turb_DP_0[which(turb_DP_0$Wavelength == "260"),]
plot_DP_0_1 <- ggplot(turb_DP_0_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_0_1$labels$colour <- "Temp(C)"
plot_DP_0_1

turb_DP_0_280 <- turb_DP_0[which(turb_DP_0$Wavelength == "280"),]
plot_DP_0_2 <- ggplot(turb_DP_0_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_0_2$labels$colour <- "Temp(C)"
plot_DP_0_2

turb_DP_0_650 <- turb_DP_0[which(turb_DP_0$Wavelength == "650"),]
plot_DP_0_3 <- ggplot(turb_DP_0_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_0_3$labels$colour <- "Temp(C)"
plot_DP_0_3

turb_DP_0_340 <- turb_DP_0[which(turb_DP_0$Wavelength == "340"),]
plot_DP_0_4 <- ggplot(turb_DP_0_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_0_4$labels$colour <- "Temp(C)"
plot_DP_0_4

turb_ratio_DP_0 <- cbind(turb_DP_0_230, data.frame(ratio=turb_DP_0_230$OD/turb_DP_0_280$OD))

plot_DP_0_230_280 <- ggplot(turb_ratio_DP_0, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_DP_0_230_280 <- plot_DP_0_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_DP_0_230_280 <- plot_DP_0_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_DP_0_230_280$labels$colour <- "Temp(C)"
plot_DP_0_230_280

grid_plot <- ggarrange(plot_DP_0_0, plot_DP_0_1, plot_DP_0_2, plot_DP_0_3, plot_DP_0_4, 
                       plot_DP_0_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for DP 0% across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_DP_001 = turb[which(turb$Sample == "DP 0.001%"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_DP_001_230 <- turb_DP_001[which(turb_DP_001$Wavelength == "230"),]
plot_DP_001_0 <- ggplot(turb_DP_001_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_001_0$labels$colour <- "Temp(C)"
plot_DP_001_0

turb_DP_001_260 <- turb_DP_001[which(turb_DP_001$Wavelength == "260"),]
plot_DP_001_1 <- ggplot(turb_DP_001_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_001_1$labels$colour <- "Temp(C)"
plot_DP_001_1

turb_DP_001_280 <- turb_DP_001[which(turb_DP_001$Wavelength == "280"),]
plot_DP_001_2 <- ggplot(turb_DP_001_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_001_2$labels$colour <- "Temp(C)"
plot_DP_001_2

turb_DP_001_650 <- turb_DP_001[which(turb_DP_001$Wavelength == "650"),]
plot_DP_001_3 <- ggplot(turb_DP_001_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_001_3$labels$colour <- "Temp(C)"
plot_DP_001_3

turb_DP_001_340 <- turb_DP_001[which(turb_DP_001$Wavelength == "340"),]
plot_DP_001_4 <- ggplot(turb_DP_001_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_001_4$labels$colour <- "Temp(C)"
plot_DP_001_4

turb_ratio_DP_001 <- cbind(turb_DP_001_230, data.frame(ratio=turb_DP_001_230$OD/turb_DP_001_280$OD))

plot_DP_001_230_280 <- ggplot(turb_ratio_DP_001, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_DP_001_230_280 <- plot_DP_001_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_DP_001_230_280 <- plot_DP_001_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_DP_001_230_280$labels$colour <- "Temp(C)"
plot_DP_001_230_280

grid_plot <- ggarrange(plot_DP_001_0, plot_DP_001_1, plot_DP_001_2, plot_DP_001_3, plot_DP_001_4, 
                       plot_DP_001_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for DP 0.001% across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_DP_005 = turb[which(turb$Sample == "DP 0.005%"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_DP_005_230 <- turb_DP_005[which(turb_DP_005$Wavelength == "230"),]
plot_DP_005_0 <- ggplot(turb_DP_005_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_005_0$labels$colour <- "Temp(C)"
plot_DP_005_0

turb_DP_005_260 <- turb_DP_005[which(turb_DP_005$Wavelength == "260"),]
plot_DP_005_1 <- ggplot(turb_DP_005_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_005_1$labels$colour <- "Temp(C)"
plot_DP_005_1

turb_DP_005_280 <- turb_DP_005[which(turb_DP_005$Wavelength == "280"),]
plot_DP_005_2 <- ggplot(turb_DP_005_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_005_2$labels$colour <- "Temp(C)"
plot_DP_005_2

turb_DP_005_650 <- turb_DP_005[which(turb_DP_005$Wavelength == "650"),]
plot_DP_005_3 <- ggplot(turb_DP_005_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_005_3$labels$colour <- "Temp(C)"
plot_DP_005_3

turb_DP_005_340 <- turb_DP_005[which(turb_DP_005$Wavelength == "340"),]
plot_DP_005_4 <- ggplot(turb_DP_005_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_005_4$labels$colour <- "Temp(C)"
plot_DP_005_4

turb_ratio_DP_005 <- cbind(turb_DP_005_230, data.frame(ratio=turb_DP_005_230$OD/turb_DP_005_280$OD))

plot_DP_005_230_280 <- ggplot(turb_ratio_DP_005, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_DP_005_230_280 <- plot_DP_005_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_DP_005_230_280 <- plot_DP_005_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_DP_005_230_280$labels$colour <- "Temp(C)"
plot_DP_005_230_280

grid_plot <- ggarrange(plot_DP_005_0, plot_DP_005_1, plot_DP_005_2, plot_DP_005_3, plot_DP_005_4, 
                       plot_DP_005_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for DP 0.005% across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_DI2 = turb[which(turb$Sample == "DI2"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_DI2_230 <- turb_DI2[which(turb_DI2$Wavelength == "230"),]
plot_DI2_0 <- ggplot(turb_DI2_230, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 230nm") +
  xlab("Concentration (1/Dilution)")
plot_DI2_0$labels$colour <- "Temp(C)"
plot_DI2_0

turb_DI2_260 <- turb_DI2[which(turb_DI2$Wavelength == "260"),]
plot_DI2_1 <- ggplot(turb_DI2_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_DI2_1$labels$colour <- "Temp(C)"
plot_DI2_1

turb_DI2_280 <- turb_DI2[which(turb_DI2$Wavelength == "280"),]
plot_DI2_2 <- ggplot(turb_DI2_280, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_DI2_2$labels$colour <- "Temp(C)"
plot_DI2_2

turb_DI2_650 <- turb_DI2[which(turb_DI2$Wavelength == "650"),]
plot_DI2_3 <- ggplot(turb_DI2_650, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_DI2_3$labels$colour <- "Temp(C)"
plot_DI2_3

turb_DI2_340 <- turb_DI2[which(turb_DI2$Wavelength == "340"),]
plot_DI2_4 <- ggplot(turb_DI2_340, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 340nm") +
  xlab("Concentration (1/Dilution)")
plot_DI2_4$labels$colour <- "Temp(C)"
plot_DI2_4

turb_ratio_DI2 <- cbind(turb_DI2_230, data.frame(ratio=turb_DI2_230$OD/turb_DI2_280$OD))

plot_DI2_230_280 <- ggplot(turb_ratio_DI2, aes(x=1/Dilution, y=ratio, group=Temp,color=as.factor(Temp)))
plot_DI2_230_280 <- plot_DI2_230_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_DI2_230_280 <- plot_DI2_230_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("230/280 Ratio, Aggregate %") + xlab("Concentration (1/Dilution)") #+ 
# labs(title = "OD Values at Increasing Temperatures, by Sample",
# subtitle = "280nm")
plot_DI2_230_280$labels$colour <- "Temp(C)"
plot_DI2_230_280

grid_plot <- ggarrange(plot_DI2_0, plot_DI2_1, plot_DI2_2, plot_DI2_3, plot_DI2_4, 
                       plot_DI2_230_280, ncol=2,nrow=3,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for DI2 across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#

turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity_22FEB2024.csv")
turb["OD"][turb["OD"] < 0] <- 0
turb_DP <- turb[grepl('DP',turb$Sample, fixed=T),]
turb_DP <- turb_DP[which(turb_DP$Sample != "DP 0% SPIN"),1:6]
turb_DP <- turb_DP[which(turb_DP$Dilution == "1"),1:6]
turb_DI2 <-turb[grepl('DI2',turb$Sample, fixed=T),]
turb_buff <- turb[grepl('FFB2',turb$Sample, fixed=T),]

fc_file = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/FC_Counts_22FEB2024.csv")
#fc_counts = fc_file[which(fc_file$Sample != "FFB2+PS20"),1:4]
fc_counts <- fc_file[,1:6]
fc_counts$Total.Particulate[which(fc_counts$Sample == "DI2")] <- fc_counts$Total.Particulate[which(fc_counts$Sample == "DI2")]*10 #adjust for dilution factor
fc_counts$Small.Particulate[which(fc_counts$Sample == "DI2")] <- fc_counts$Small.Particulate[which(fc_counts$Sample == "DI2")]*10 #adjust for dilution factor
fc_counts$Large.Particulate[which(fc_counts$Sample == "DI2")] <- fc_counts$Large.Particulate[which(fc_counts$Sample == "DI2")]*10 #adjust for dilution factor
#remove spin
fc_counts <- fc_counts[which(fc_counts$Sample != "DP0 SPIN"),1:6]
fc_counts_DP <- fc_counts[which(fc_counts$Sample != "DI2"),1:6]


turb_230 = turb_DP[which(turb_DP$Wavelength == "230"),setdiff(1:6,3)]
turb_230 = turb_230[which(turb_230$Dilution == "1"),-which(colnames(turb_230) == "Dilution")]
comb_turb_230 <- cbind(turb_230, totalParticulate = fc_counts_DP$Total.Particulate, 
                       smallParticulate = fc_counts_DP$Small.Particulate,
                       largeParticulate = fc_counts_DP$Large.Particulate)
#turb_230_melt <- melt(turb_230,id= c("Sample","Temp"))
# 
# plot_230 <- ggplot(turb_230_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
# plot_230 <- plot_230 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
# plot_230 <- plot_230 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
#   ylab("OD at 230nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
# #  subtitle = "230nm")
# plot_230
#from https://finchstudio.io/blog/ggplot-dual-y-axes/
max_first  <- max(comb_turb_230$RawOD)   # Specify max of first y axis
max_second <- max(comb_turb_230$totalParticulate) # Specify max of second y axis
min_first  <- min(comb_turb_230$RawOD)   # Specify min of first y axis
min_second <- min(comb_turb_230$totalParticulate) # Specify min of second y axis

# scale and shift variables calculated based on desired mins and maxes
scale = (max_second - min_second)/(max_first - min_first)
shift = min_first - min_second

# Function to scale secondary axis
scale_function <- function(x, scale, shift){
  return ((x)*scale - shift)
}

# Function to scale secondary variable values
inv_scale_function <- function(x, scale, shift){
  return ((x + shift)/scale)
}

p2 <- ggplot(comb_turb_230,aes(x=Temp)) + 
  geom_line(aes(y=RawOD,alpha=0.5,group=Sample,color=Sample)) + geom_point(aes(y=RawOD),alpha=0.5) + 
  geom_line(aes(y=inv_scale_function(totalParticulate,scale,shift),alpha=0.5,group=Sample,color=Sample)) + 
  geom_point(aes(y=inv_scale_function(totalParticulate,scale,shift)),color="white",alpha=0.5) + 
  scale_color_brewer(palette="Set1") + guides(alpha="none") + 
  scale_y_continuous("Raw OD", #limits=c(min_second,max_second),
                     sec.axis = sec_axis(~scale_function(.,scale,shift),name="Total Particulate")) +
  ggtitle("Total Particulate and Raw OD") + 
  theme(plot.title=element_text(hjust = 0.5)) + xlab("Temperature") +
  ylab("Raw OD") + labs(colour="")
  #scale_x_continuous(breaks=kVals[c(FALSE,TRUE,FALSE)],labels=xLabel[c(FALSE,TRUE,FALSE)]) +
  
  #theme(axis.text.x = element_text(angle=90,hjust=0.5))
p2

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
library(protti)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggbreak)
library(ggforce)
pDiff_df <- read_protti("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/pDiff_plotTable.csv")
pDiff_report <- read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/pDiff_reportTable.csv",
                         header = TRUE)
colnames(pDiff_report) <- c("Measurement","DP 0%", "DP 0.001%", "DP 0.005%", "%CV")
# Convert into long format
pDiff_long <- pDiff_df %>%
  pivot_longer(
    cols = starts_with("p_diff_"),
    names_to = "Sample",
    names_prefix = "p_diff_",
    values_to = "pDiff"
  )

plot <- ggplot(pDiff_long, aes(x=fct_rev(fct_reorder(category,pDiff)),y = pDiff*100, fill=Sample))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', width = 0.6) + theme_dark() 
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ylab("Percent Difference") +
  labs(title = "Percent Difference between PS20 Concentrations, by Measurement",
       subtitle = "Filtered DP Samples")
plot

plot <- ggplot(pDiff_long, aes(x=fct_rev(fct_reorder(category,pDiff)),y = pDiff*100, fill=Sample))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', width = 0.6) + theme_dark() 
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ ylab("Percent Difference") +
  labs(title = "Percent Difference between PS20 Concentrations, by Measurement",
       subtitle = "Filtered DP Samples") +
  facet_zoom(ylim = c(0,10))
plot
