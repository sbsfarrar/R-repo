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
#library(wesanderson)

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
turb_260 = turb[which(turb$Wavelength == "260"),setdiff(1:5,3)]
turb_260 = turb_260[which(turb_260$Dilution == "1"),-which(colnames(turb_260) == "Dilution")]
turb_260_melt <- melt(turb_260,id= c("Sample","Temp"))

plot_260 <- ggplot(turb_260_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_260 <- plot_260 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_260 <- plot_260 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 260nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
                        #  subtitle = "260nm")
plot_260




turb_280 = turb[which(turb$Wavelength == "280"),setdiff(1:5,3)]
turb_280 = turb_280[which(turb_280$Dilution == "1"),-which(colnames(turb_280) == "Dilution")]
turb_280_melt <- melt(turb_280,id= c("Sample","Temp"))

plot_280 <- ggplot(turb_280_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_280 <- plot_280 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_280 <- plot_280 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 280nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
                         # subtitle = "280nm")
plot_280


turb_650 = turb[which(turb$Wavelength == "650"),setdiff(1:5,3)]
turb_650 = turb_650[which(turb_650$Dilution == "1"),-which(colnames(turb_650) == "Dilution")]
turb_650_melt <- melt(turb_650,id= c("Sample","Temp"))

plot_650 <- ggplot(turb_650_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_650 <- plot_650 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_650 <- plot_650 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 650nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
                         # subtitle = "650nm")
plot_650


turb_980 = turb[which(turb$Wavelength == "980"),setdiff(1:5,3)]
turb_980 = turb_980[which(turb_980$Dilution == "1"),-which(colnames(turb_980) == "Dilution")]
turb_980_melt <- melt(turb_980,id= c("Sample","Temp"))

plot_980 <- ggplot(turb_980_melt, aes(x=Temp, y=value, group=Sample,color=Sample))
plot_980 <- plot_980 + geom_line(size = 0.8, alpha = 0.5) + geom_point(size = 1.5)
plot_980 <- plot_980 + theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  ylab("OD at 980nm") + xlab("Temperature")#+ labs(title = "OD Values at Increasing Temperatures, by Sample",
                         # subtitle = "980nm")
plot_980



grid_plot <- ggarrange(plot_260, plot_280, plot_650, plot_980, ncol=2,nrow=2,
          common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at Increasing Temperatures, by Sample", 
                                       face = "bold", size = 14))

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
turb_25 = turb[which(turb$Temp == "25"),setdiff(1:5,2)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

turb_25_260 <- turb_25[which(turb_25$Wavelength == "260"),]
plot_25_1 <- ggplot(turb_25_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_25_1

turb_25_280 <- turb_25[which(turb_25$Wavelength == "280"),]
plot_25_2 <- ggplot(turb_25_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_25_2

turb_25_650 <- turb_25[which(turb_25$Wavelength == "650"),]
plot_25_3 <- ggplot(turb_25_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_25_3

turb_25_980 <- turb_25[which(turb_25$Wavelength == "980"),]
plot_25_4 <- ggplot(turb_25_980, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_25_4
  

grid_plot <- ggarrange(plot_25_1, plot_25_2, plot_25_3, plot_25_4, ncol=2,nrow=2,
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

turb_35_260 <- turb_35[which(turb_35$Wavelength == "260"),]
plot_35_1 <- ggplot(turb_35_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_35_1

turb_35_280 <- turb_35[which(turb_35$Wavelength == "280"),]
plot_35_2 <- ggplot(turb_35_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_35_2

turb_35_650 <- turb_35[which(turb_35$Wavelength == "650"),]
plot_35_3 <- ggplot(turb_35_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_35_3

turb_35_980 <- turb_35[which(turb_35$Wavelength == "980"),]
plot_35_4 <- ggplot(turb_35_980, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_35_4


grid_plot <- ggarrange(plot_35_1, plot_35_2, plot_35_3, plot_35_4, ncol=2,nrow=2,
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

turb_45_260 <- turb_45[which(turb_45$Wavelength == "260"),]
plot_45_1 <- ggplot(turb_45_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_45_1

turb_45_280 <- turb_45[which(turb_45$Wavelength == "280"),]
plot_45_2 <- ggplot(turb_45_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_45_2

turb_45_650 <- turb_45[which(turb_45$Wavelength == "650"),]
plot_45_3 <- ggplot(turb_45_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_45_3

turb_45_980 <- turb_45[which(turb_45$Wavelength == "980"),]
plot_45_4 <- ggplot(turb_45_980, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_45_4


grid_plot <- ggarrange(plot_45_1, plot_45_2, plot_45_3, plot_45_4, ncol=2,nrow=2,
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

turb_55_260 <- turb_55[which(turb_55$Wavelength == "260"),]
plot_55_1 <- ggplot(turb_55_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_55_1

turb_55_280 <- turb_55[which(turb_55$Wavelength == "280"),]
plot_55_2 <- ggplot(turb_55_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_55_2

turb_55_650 <- turb_55[which(turb_55$Wavelength == "650"),]
plot_55_3 <- ggplot(turb_55_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_55_3

turb_55_980 <- turb_55[which(turb_55$Wavelength == "980"),]
plot_55_4 <- ggplot(turb_55_980, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_55_4


grid_plot <- ggarrange(plot_55_1, plot_55_2, plot_55_3, plot_55_4, ncol=2,nrow=2,
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

turb_65_260 <- turb_65[which(turb_65$Wavelength == "260"),]
plot_65_1 <- ggplot(turb_65_260, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 260nm") +
  xlab("Concentration (1/Dilution)")
plot_65_1

turb_65_280 <- turb_65[which(turb_65$Wavelength == "280"),]
plot_65_2 <- ggplot(turb_65_280, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 280nm") +
  xlab("Concentration (1/Dilution)")
plot_65_2

turb_65_650 <- turb_65[which(turb_65$Wavelength == "650"),]
plot_65_3 <- ggplot(turb_65_650, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 650nm") +
  xlab("Concentration (1/Dilution)")
plot_65_3

turb_65_980 <- turb_65[which(turb_65$Wavelength == "980"),]
plot_65_4 <- ggplot(turb_65_980, aes(x=1/Dilution, y=OD, group=Sample, color=Sample)) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) + ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_65_4


grid_plot <- ggarrange(plot_65_1, plot_65_2, plot_65_3, plot_65_4, ncol=2,nrow=2,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values at 65C across Selected Wavelengths, by Sample", 
                                           face = "bold", size = 14))


#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
turb = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Turbidity/Turbidity.csv")
temp_colors <- c("blue","green", "yellow", "orange", "red")
turb_H2O = turb[which(turb$Sample == "H2O"),setdiff(1:5,1)]
#turb_25 = turb_25[which(turb_25$Wavelength == "1"),-which(colnames(turb_25) == "Dilution")]
#turb_25_melt <- melt(turb_25,id= c("Sample","Temp",))

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

turb_H2O_980 <- turb_H2O[which(turb_H2O$Wavelength == "980"),]
plot_H2O_4 <- ggplot(turb_H2O_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_H2O_4$labels$colour <- "Temp(C)"
plot_H2O_4


grid_plot <- ggarrange(plot_H2O_1, plot_H2O_2, plot_H2O_3, plot_H2O_4, ncol=2,nrow=2,
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

turb_FFB2_260 <- turb_FFB2[which(turb_FFB2$Wavelength == "260"),]
plot_FFB2_1 <- ggplot(turb_FFB2_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=1) + geom_point(size=1.5) +
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

turb_FFB2_980 <- turb_FFB2[which(turb_FFB2$Wavelength == "980"),]
plot_FFB2_4 <- ggplot(turb_FFB2_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_4$labels$colour <- "Temp(C)"
plot_FFB2_4


grid_plot <- ggarrange(plot_FFB2_1, plot_FFB2_2, plot_FFB2_3, plot_FFB2_4, ncol=2,nrow=2,
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

turb_FFB2_PS20_260 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "260"),]
plot_FFB2_PS20_1 <- ggplot(turb_FFB2_PS20_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=1) + geom_point(size=1.5) +
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

turb_FFB2_PS20_980 <- turb_FFB2_PS20[which(turb_FFB2_PS20$Wavelength == "980"),]
plot_FFB2_PS20_4 <- ggplot(turb_FFB2_PS20_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_FFB2_PS20_4$labels$colour <- "Temp(C)"
plot_FFB2_PS20_4


grid_plot <- ggarrange(plot_FFB2_PS20_1, plot_FFB2_PS20_2, plot_FFB2_PS20_3, plot_FFB2_PS20_4, ncol=2,nrow=2,
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

turb_DP_0_260 <- turb_DP_0[which(turb_DP_0$Wavelength == "260"),]
plot_DP_0_1 <- ggplot(turb_DP_0_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=1) + geom_point(size=1.5) +
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

turb_DP_0_980 <- turb_DP_0[which(turb_DP_0$Wavelength == "980"),]
plot_DP_0_4 <- ggplot(turb_DP_0_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_0_4$labels$colour <- "Temp(C)"
plot_DP_0_4


grid_plot <- ggarrange(plot_DP_0_1, plot_DP_0_2, plot_DP_0_3, plot_DP_0_4, ncol=2,nrow=2,
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

turb_DP_001_260 <- turb_DP_001[which(turb_DP_001$Wavelength == "260"),]
plot_DP_001_1 <- ggplot(turb_DP_001_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=1) + geom_point(size=1.5) +
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

turb_DP_001_980 <- turb_DP_001[which(turb_DP_001$Wavelength == "980"),]
plot_DP_001_4 <- ggplot(turb_DP_001_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_001_4$labels$colour <- "Temp(C)"
plot_DP_001_4


grid_plot <- ggarrange(plot_DP_001_1, plot_DP_001_2, plot_DP_001_3, plot_DP_001_4, ncol=2,nrow=2,
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

turb_DP_005_260 <- turb_DP_005[which(turb_DP_005$Wavelength == "260"),]
plot_DP_005_1 <- ggplot(turb_DP_005_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=1) + geom_point(size=1.5) +
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

turb_DP_005_980 <- turb_DP_005[which(turb_DP_005$Wavelength == "980"),]
plot_DP_005_4 <- ggplot(turb_DP_005_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_DP_005_4$labels$colour <- "Temp(C)"
plot_DP_005_4


grid_plot <- ggarrange(plot_DP_005_1, plot_DP_005_2, plot_DP_005_3, plot_DP_005_4, ncol=2,nrow=2,
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

turb_DI2_260 <- turb_DI2[which(turb_DI2$Wavelength == "260"),]
plot_DI2_1 <- ggplot(turb_DI2_260, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=1) + geom_point(size=1.5) +
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

turb_DI2_980 <- turb_DI2[which(turb_DI2$Wavelength == "980"),]
plot_DI2_4 <- ggplot(turb_DI2_980, aes(x=1/Dilution, y=OD, group=Temp, color=as.factor(Temp))) + 
  geom_line(size = 0.8, alpha=0.5) + geom_point(size=1.5) +
  scale_color_manual(values = temp_colors,aesthetics = "color") +
  ylab("OD at 980nm") +
  xlab("Concentration (1/Dilution)")
plot_DI2_4$labels$colour <- "Temp(C)"
plot_DI2_4


grid_plot <- ggarrange(plot_DI2_1, plot_DI2_2, plot_DI2_3, plot_DI2_4, ncol=2,nrow=2,
                       common.legend = TRUE, legend="bottom")
annotate_figure(grid_plot, top = text_grob("OD Values for DI2 across Selected Wavelengths, by Temperature", 
                                           face = "bold", size = 14))
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
