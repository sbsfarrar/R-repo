##############################INITIALIZE###########################################
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyFiles)
library(htmlwidgets)
library(R.matlab)
library(ggplot2)
library(ggforce)
# library(sp)
library(tidyr)
library(dplyr)
# library(plot3D)
library(viridis)
# library(oce)
library(reshape2)
library(xml2)
library(plotly)
library(rgl)
library(data.table)
library(rvest)
library(rsvg)
library(profvis)
library(scales)
library(plotwidgets)
open3d(useNULL = TRUE)
scene <- scene3d()
rgl.close()
library(ggpubr)
library(RColorBrewer)
library(ggridges)
library(svgtools)
library(gridExtra)
library(ggpubr)
### DATA INIT
sliceZ <- 0
voxel_volume <- 20 * 20 * 20
factor_names <-
  c(
    "Akt",
    "AMPK",
    "cMyc",
    "HIF",
    "mTOR",
    "NOX",
    "p53",
    "PDK",
    "PI3K",
    "PTEN",
    "RAS",
    "SOD",
    "VEGF",
    "GluT1",
    "HK",
    "G6PD",
    "GPI",
    "PFKFB2",
    "PFK1",
    "ALD",
    "TPI",
    "GAPDH",
    "PGK",
    "PHGDH",
    "PGAM",
    "ENO",
    "PKM2",
    "PDH",
    "ACC",
    "LDH",
    "Glucose",
    "G6P",
    "F6P",
    "FBP",
    "G3P",
    "DHAP",
    "13BPG",
    "3PG",
    "2PG",
    "PEP",
    "Pyruvate",
    "Lactate",
    "R5P",
    "F26BP",
    "Serine",
    "Citrate",
    "AMP",
    "ADP",
    "ATP",
    "NAD",
    "NADH",
    "complex2",
    "ROS"
  )
setwd("C:/Users/StevenSummey/Documents/CancerSimulationPhysiCell/SHINY_APP")
load("C:/Users/StevenSummey/Documents/CancerSimulationPhysiCell/SHINY_APP/GammaNamesList.RData")
folder <- 'C:/Users/StevenSummey/Desktop/VB_Share/May_29_2024_1519/home/vboxuser/PhysiCell/output/'
measured_data <- c()
cell_measured_data_name <- c()
envir_data <- c()
NanToZero <- function(x) {
  if (is.nan(x)) {
    return(0)
  }
  return(x)
}
## R PLOT FUNCTIONS
GET_TIMESTEP <- function(folder_path = folder) {
  return(as.numeric(gsub(".xml", "", gsub(
    "output", "", dir(folder_path, pattern = "output.*.xml")
  ))))
}
GET_MEDIUM_OBJECT <- function(TIME , folder_path = folder) {
  TIME <- TIME + 1
  files.names <- dir(folder_path, pattern = "output.*.xml")
  path_file <-
    paste(folder_path,
          gsub(".xml", "_microenvironment0.mat", files.names[TIME]),
          sep = "")
  MEDIUM <- readMat(path_file)
  MEDIUM <- MEDIUM$multiscale.microenvironment
  MEDIUM <- t(MEDIUM)
  MEDIUM <- data.frame(MEDIUM, TIME)
  if ("H" %in% envir_data) {
    MEDIUM$pH <- MEDIUM[, 0]
  }
  colnames(MEDIUM) <-
    c("x", "y", "z", "constant", envir_data, "TIME")
  MEDIUM$pH <-
    -log10(MEDIUM$H / (1e-15 * 1000)) * 1e-15 # multiplied back by the volume of voxel
  MEDIUM$Distance <-
    sqrt((MEDIUM$x - 0) ^ 2 + (MEDIUM$y - 0) ^ 2 + (MEDIUM$z - 0) ^
           2 * 1.0)
  MEDIUM$oxygen <- MEDIUM$oxygen / (1.39e-3)
  
  return(MEDIUM)
}
GET_CELL_OBJECT <- function(TIME, path = folder, firstTime=FALSE) {
  TIME <- TIME + 1
  print(paste0("CELL_OBJECT FUNCTION ACCESSED AT TIME: ", TIME-1))
  files.names2 <-
    dir(path, pattern = "output.*_cells_physicell.mat")
  
  path_file <-
    paste(path, files.names2[TIME], sep = "")
  cells_temp <- readMat(path_file)
  cells_temp <- cells_temp$cells
  cells_temp <- t(cells_temp)
  cells_temp <- as.data.frame(cells_temp)
  cells_temp <-
    cells_temp %>% mutate(PHENOTYPE = cells_temp[1]) # temp column
  cells_temp <-
    cells_temp %>% mutate(PHENOTYPE2 = cells_temp[1]) # temp column
  colnames(cells_temp) <- cell_measured_data_name
  cells_temp$distance <-
    with(cells_temp, sqrt((cells_temp$position.x ^ 2) + (position.y ^ 2) + (position.z ^
                                                                              2)))
  # cells_temp$distance <-with(cells_temp, (cells_temp$position.y)-min(cells_temp$position.y))
  
  cells_temp$simulation_time <- TIME
  if(firstTime==FALSE){
    cells_temp <- cells_temp %>% filter(cell_type == 0)
    
    if ("metabolites_parameters6" %in% colnames(cells_temp)) {
      cells_temp$metabolites_parameters6 <-
        cells_temp$metabolites_parameters6 / (1.39e-3) #from mM back to mmHg
    }
    if ("metabolites_parameters118" %in% colnames(cells_temp)) {
      cells_temp$metabolites_parameters118 <- -log10(cells_temp$metabolites_parameters118)
    }
    
    
    if ("r21" %in% colnames(cells_temp)) {
      cells_temp <- cells_temp %>% mutate(
        PHENOTYPE = case_when(
          cycle_model != 9999 ~ "DEAD",
          r12 <= 0 ~ "RESPIRATION (++)",
          r12 > 0 &
            r12 > r18 & r12 <= 10 * r18 ~ "FERMENTATION (+) Respiration (-)",
          r12 > 0 & r12 > r18 & r12 > 10 * r18  ~ "FERMENTATION (++)",
          r12 > 0 &
            r18 > r12 & r18 <= 10 * r12 ~ "RESPIRATION (+) Fermentation (-)",
          r12 > 0 & r18 > r12 & r18 > 10 * r12  ~ "RESPIRATION (++)"
        )
      )
      
      cells_temp <- cells_temp %>% mutate(
        PHENOTYPE2 = case_when(
          cycle_model != 9999 ~ "DEAD",
          r21 <= 0  ~ "REVERSE WARBURG",
          metabolites_parameters6 > 4 & r21 > 0  ~ "WARBURG MODERN",
          metabolites_parameters6 <= 4 & r21 > 0  ~ "WARBURG ORIGINAL",
          metabolites_parameters6 <= 4 & r21 > 0 ~ PHENOTYPE
        )
      )
      
      ##STD DEFINITION OF WARBURG
      # cells_temp <- cells_temp %>% mutate(PHENOTYPE = case_when(
      #   r21 > 0 & metabolites_parameters6 > 8 ~ "WARBURG",
      #   r21 <= 0 & metabolites_parameters6 > 8 ~ "REVERSE WARBURG",
      #   r21 > 0 & metabolites_parameters6 <= 8 ~ "STD GLYCOLYTIC - HYPOXIC",
      #   r21 <= 0 & metabolites_parameters6 <= 8 ~ "STD MYTO - HYPOXIC"
      # ))
    }
  }
  return(cells_temp)
}
GET_ENV_ACROSS_TIME <- function(PARAMETER, sample, radius) {
  cols <- c("group", "mean","median", "q1", "q3", "n", "std","radius", "time")
  values <- data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(values) <- cols
  for (step in 1:length(sample)) {
    cut_levels <- radius
    medium_object <- GET_MEDIUM_OBJECT(sample[step])
    object <-
      medium_object %>%
      group_by(group = cut(
        Distance,
        breaks = seq(0, ceiling(max(Distance)), ceiling(max(Distance)) / cut_levels),
        include.lowest = TRUE
      )) %>%
      summarise(
        mean=mean(!!sym(PARAMETER)),
        median = median(!!sym(PARAMETER)),
        q1 = quantile(!!sym(PARAMETER), 0.25),
        q3 = quantile(!!sym(PARAMETER), 0.75),
        n = n(),
        std=sd(!!sym(PARAMETER))
      )
    
    object$radius <- as.numeric(rownames(object))
    object$time <- sample[step]
    
    values <- rbind(values, as.data.frame(object))
  }
  colnames(values) <- cols
  return(values)
}
GET_ALL_CELL_SIM <- function(from, to, by, path = folder) {
  datacollect <- GET_CELL_OBJECT(from, path = path)
  datacollect$simtime <- from
  for (timestep in seq(from + by, to, by)) {
    tempdata <- GET_CELL_OBJECT(timestep, path = path)
    tempdata$simtime <- timestep
    datacollect <- rbind(datacollect, as.data.frame(tempdata))
  }
  return(datacollect)
}
GET_ALL_SIM <- function(simlist) {
  N_EVAL <- 6
  datacollect <-
    GET_ALL_CELL_SIM(0, last(GET_TIMESTEP(simlist[1])), floor(last(GET_TIMESTEP(simlist[1])) /
                                                                (N_EVAL - 1)), simlist[1])
  datacollect$sim <- simlist[1]
  for (sim in simlist[c(-1)]) {
    tempdata <-
      GET_ALL_CELL_SIM(0, last(GET_TIMESTEP(sim)), floor(last(GET_TIMESTEP(sim)) /
                                                           (N_EVAL - 1)), sim)
    tempdata$sim <- sim
    datacollect <- bind_rows(datacollect, as.data.frame(tempdata))
  }
  return(datacollect)
}
GET_ALL_SIM_SAMPLE <- function(simlist, sampling) {
  datacollect <-
    GET_CELL_OBJECT(last(GET_TIMESTEP(simlist[1])), simlist[1])
  datacollect <- datacollect %>% filter(cycle_model == 9999)
  datacollect <-
    datacollect %>% sample_n(min(sampling, nrow(datacollect)))
  datacollect$sim <- simlist[1]
  for (sim in simlist[c(-1)]) {
    tempdata <- GET_CELL_OBJECT(last(GET_TIMESTEP(sim)), sim)
    tempdata <- tempdata %>% filter(cycle_model == 9999)
    tempdata <- tempdata %>% sample_n(min(sampling, nrow(tempdata)))
    tempdata$sim <- sim
    datacollect <- bind_rows(datacollect, as.data.frame(tempdata))
  }
  return(datacollect)
}
GET_ALL_CELL_ACROSS_TIME <- function(PARAMETER, sample, radius,RIDGE=F,path=folder) {
  cols <- c("group", "mean","median", "q1", "q3", "sd","n", "radius", "time")
  if(RIDGE==T){
    cols <- c("value","radius", "time")
  }
  values <- data.frame(matrix(ncol = length(cols), nrow = 0))
  colnames(values) <- cols
  for (step in 1:length(sample)) {
    cut_levels <- radius
    cell_object <- GET_CELL_OBJECT(sample[step],path)
    
    if(RIDGE==T){
      object<-as.data.frame(cell_object[[PARAMETER]])
      colnames(object)<-c("value")
    }else{
      object <-
        cell_object %>%
        group_by(group = cut(
          distance,
          breaks = seq(0, ceiling(max(distance)), ceiling(max(distance)) / cut_levels),
          include.lowest = TRUE
        )) %>%
        summarise(
          mean = mean(!!sym(PARAMETER)),
          median=median(!!sym(PARAMETER)),
          q1 = quantile(!!sym(PARAMETER), 0.25),
          q3 = quantile(!!sym(PARAMETER), 0.75),
          sd=sd(!!sym(PARAMETER)),
          n = n()
        )
    }
    
    object$radius <- as.numeric(rownames(object))
    object$time <- sample[step]
    
    values <- rbind(values, as.data.frame(object))
  }
  colnames(values) <- cols
  return(values)
}
# across<-GET_ALL_CELL_ACROSS_TIME("Lactate",seq(from=0,to=max(GET_TIMESTEP()),by=5))
EXPORT_CELL_OBJECT_CSV <- function() {
  steps <- GET_TIMESTEP()
  for (value in steps) {
    object <- GET_CELL_OBJECT(value)
    write.table(
      x = object,
      sep = ";",
      file = paste0(folder, paste0(
        paste0("output", sprintf("%07d", value + 1)), "_cell_physicell.csv"
      )),
      row.names = FALSE,
      col.names = FALSE,
      append = FALSE
    )
  }
}
PLOT_MEDIUM_TIME <- function(PARAMETER, SAMPLE, RADIUS) {
  across <-
    GET_ENV_ACROSS_TIME(PARAMETER, seq(
      from = 0,
      to = max(GET_TIMESTEP()),
      by = SAMPLE
    ), RADIUS)
  #
  # fig <- plot_ly(across, x = ~time, y = across$q1, split=~radius, color=~radius, type = 'scatter', mode = 'lines',
  #                line = list(color = 'transparent'),
  #                showlegend = FALSE, name = 'High')
  #
  # fig <- fig %>% add_trace(y = across$q3, split=~radius, color=~radius, type = 'scatter', mode = 'lines',
  #                          fill = 'tonexty', line = list(color = 'transparent'),
  #                          showlegend = FALSE, name = 'Low')
  
  fig <-
    plot_ly(
      across,
      x = ~ time,
      y = ~ mean / (1e-15),
      split = ~ radius,
      color = ~ radius,
      type = "scatter",
      mode = "lines",
      name = ~ paste0("Couche ", formatC(
        radius,
        width = 2, flag = "0"
      ))
    )
  
  fig <-
    fig %>% layout(
      title = "Median, High and Low of selected parameter",
      paper_bgcolor = "rgb(255,255,255)",
      plot_bgcolor = "rgb(229,229,229)",
      yaxis = list(type = "linear")
    )
  fig
  return(fig)
}
PLOT_CELL_TIME <- function(PARAMETER, SAMPLE, RADIUS,path=folder,force=FALSE) {
  
  across <-
    GET_ALL_CELL_ACROSS_TIME(PARAMETER, seq(
      from = 0,
      to = max(GET_TIMESTEP(path)),
      by = SAMPLE
    ), RADIUS,RADIUS==1,path)
  print(across)
  
  #
  # fig <- plot_ly(across, x = ~time, y = across$q1, split=~radius, color=~radius, type = 'scatter', mode = 'lines',
  #                line = list(color = 'transparent'),
  #                showlegend = FALSE, name = 'High')
  #
  # fig <- fig %>% add_trace(y = across$q3, split=~radius, color=~radius, type = 'scatter', mode = 'lines',
  #                          fill = 'tonexty', line = list(color = 'transparent'),
  #                          showlegend = FALSE, name = 'Low')
  
  if(RADIUS==1 & force==TRUE){
    column<-across$value
    xmin<-min(column[!column %in% boxplot.stats(column,coef = 0)$out])
    xmax<-max(column[!column %in% boxplot.stats(column,coef = 1.5)$out])
    fig<-  ggplot(across, aes(x = value, y = as.factor(time), fill = ..x..)) +
      geom_density_ridges_gradient(scale = 3) +
      scale_fill_viridis(option = "C",direction = -1) +
      theme_ridges() +
      xlim(xmin,xmax)+
      theme(legend.position = "none")
  }
  else{
    fig <-
      plot_ly(
        across,
        x = ~ time,
        y = ~ median,
        split = ~ radius,
        color = ~ radius,
        type = "scatter",
        mode = "lines",
        name = ~ paste0("Couche ", formatC(
          radius,
          width = 2, flag = "0"
        ))
      )
    
    fig <-
      fig %>% layout(
        title = "Median, High and Low of selected parameter",
        paper_bgcolor = "rgb(255,255,255)",
        plot_bgcolor = "rgb(229,229,229)",
        yaxis = list(type = "linear")
      )
    fig
  }
  
  return(fig)
}
plotCellPopAcrossTime <- function(sample,path=folder) {
  living_cells <- c()
  quiescent_cells <- c()
  dead_cells <- c()
  all_cells <- c()
  for (value in 1:length(sample)) {
    all_cells_at_step <- GET_CELL_OBJECT(sample[value],path = path)$current_phase
    living_cell_at_step <-
      length(all_cells_at_step[(all_cells_at_step %in% c(106, 107, 108, 109))])
    living_cells <- c(living_cells, living_cell_at_step)
    
    quiescent_cells_at_step <-length(all_cells_at_step[(all_cells_at_step %in% c(105))])
    quiescent_cells<- c(quiescent_cells, quiescent_cells_at_step)
    
    all_cells <- c(all_cells, length(all_cells_at_step))
    dead_cells <-
      c(dead_cells, length(all_cells_at_step) - living_cell_at_step - quiescent_cells_at_step)
  }
  fig <-
    plot_ly(
      x = sample,
      y = all_cells,
      type = "scatter",
      mode = "lines",
      name = "All cells"
    )
  fig <- fig %>% add_trace(y = living_cells, name = "Proliferating cells")
  fig <- fig %>% add_trace(y = quiescent_cells, name = "Quiescent cells")
  fig <- fig %>% add_trace(y = dead_cells, name = "Dead cells")
  fig
}
plotCellPhenotypeAcrossTime <- function(sample,path=folder,PARAMETER) {
  tempobject<-GET_CELL_OBJECT(sample[1],path = path)
  all_cells <- as.data.frame(table(tempobject[[PARAMETER]]))
  DEADs<-as.numeric(table(tempobject$PHENOTYPE)["DEAD"])
  if(is.na(DEADs)==TRUE){
    DEADs<-0
  }
  all_cells$TIME<-sample[1]
  all_cells$total<-length(tempobject$ID)-DEADs
  for (value in sample[-1]) {
    tempobject<-GET_CELL_OBJECT(value,path = path)
    DEADs<-as.numeric(table(tempobject$PHENOTYPE)["DEAD"])
    if(is.na(DEADs)==TRUE){
      DEADs<-0
    }
    all_cells_at_step <- as.data.frame(table(tempobject[[PARAMETER]]))
    all_cells_at_step$TIME<-value
    all_cells_at_step$total<-length(tempobject$ID)-DEADs
    all_cells <- bind_rows(all_cells, all_cells_at_step)
  }
  print(all_cells)
  fig <-
    plot_ly(all_cells,
            x = ~TIME,
            y = ~(Freq/total),
            color=~as.factor(Var1),
            name=~as.factor(Var1),
            type = "scatter",
            mode = "lines"
    )
  fig
}
EnvironmentMAPGG <-
  function(PARAMETER,
           TIME,
           GLOBAL_CURRENT_CELL_OBJECT) {
    MEDIUM_OBJECT <- GET_MEDIUM_OBJECT(TIME)
    #cells_temp <- GLOBAL_CURRENT_CELL_OBJECT
    #cells_temp$distance <-with(cells_temp, sqrt((cells_temp$position.x ^ 2) + (position.y ^ 2) + (position.z ^2)))
    #max_radius <- max(cells_temp$distance)
    MEDIUM_OBJECT <-MEDIUM_OBJECT[MEDIUM_OBJECT$z == -min(abs(MEDIUM_OBJECT$z)),]
    MEDIUM_OBJECT[, PARAMETER] <- MEDIUM_OBJECT[, PARAMETER] / (1e-15)
    rangePARAMETER <- range(MEDIUM_OBJECT[, PARAMETER])
    
    scalePARAMETER <-
      scale_fill_viridis_c(PARAMETER,
                           option = "plasma",
                           #limits = c(rangePARAMETER[1], rangePARAMETER[2]),
                           #breaks=c(rangePARAMETER[1],rangePARAMETER[1] +(rangePARAMETER[2]-rangePARAMETER[1])/2 , rangePARAMETER[2]),
      )
    
    p <- ggplot(data = MEDIUM_OBJECT,
                aes(x = MEDIUM_OBJECT$x,
                    y = MEDIUM_OBJECT$y,
                    fill = MEDIUM_OBJECT[, PARAMETER])) + geom_tile() + scalePARAMETER +
      #annotate(
      #   "path",
      #   x = 0 + max_radius * cos(seq(0, 2 * pi, length.out = 100)),
      #   y = 0 + max_radius * sin(seq(0, 2 * pi, length.out = 100)),
      #   color = "white",
      #   alpha = 0.5
      # ) +
      coord_equal() +
      labs(colour = "", x = NULL, y = NULL) +
      theme(
        legend.position = c(1, 0),
        legend.justification = c(1, 0),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
    return(p)
  }
EnvironmentMAP <-
  function(PARAMETER,
           TIME,
           GLOBAL_CURRENT_CELL_OBJECT) {
    MEDIUM_OBJECT <- GET_MEDIUM_OBJECT(TIME)
    cells_temp <- GLOBAL_CURRENT_CELL_OBJECT
    cells_temp$distance <-
      with(cells_temp, sqrt((cells_temp$position.x ^ 2) + (position.y ^ 2) + (position.z ^
                                                                                2)))
    max_radius <- max(cells_temp$distance)
    MEDIUM_OBJECT <-
      MEDIUM_OBJECT[MEDIUM_OBJECT$z == -min(abs(MEDIUM_OBJECT$z)),]
    MEDIUM_OBJECT[, PARAMETER] <- MEDIUM_OBJECT[, PARAMETER] / (1e-15)
    rangePARAMETER <- range(MEDIUM_OBJECT[, PARAMETER])
    
    scalePARAMETER <-
      scale_fill_viridis_c(
        PARAMETER,
        option = "plasma",
        limits = c(rangePARAMETER[1], rangePARAMETER[2])
      )
    themeNoLegends <-
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.justification = c(1, 0),
        legend.position = c(1, 0)
      )
    
    p <- plot_ly(
      data = MEDIUM_OBJECT,
      x = MEDIUM_OBJECT$x,
      y = MEDIUM_OBJECT$y,
      colors = "plasma",
      z = MEDIUM_OBJECT[, PARAMETER],
      type = "heatmap"
    ) %>%
      layout(shapes = list(
        list(
          type = "circle",
          xref = "x",
          x0 = -max_radius,
          x1 = max_radius,
          yref = "y",
          y0 = -max_radius,
          y1 = max_radius,
          fillcolor = "rgba(50, 20, 90,0.0)",
          line = list(color = "rgb(255, 255, 255)"),
          opacity = 0.5
        )
      )) %>%
      config(displayModeBar = F) %>%
      layout(
        margin = list(l = 0, b = 0, t = 0),
        xaxis = list(
          showticklabels = FALSE,
          showgrid = FALSE,
          zeroline = FALSE,
          showline = FALSE,
          ticks = ""
        ),
        yaxis = list(
          showticklabels = FALSE,
          showgrid = FALSE,
          zeroline = FALSE,
          showline = FALSE,
          ticks = "",
          scaleanchor = "x"
        ),
        plot_bgcolor = "rgba(0,0,0,0)",
        colorbar = list()
      )
    return(p)
  }
plotCellatTime <- function(PARAMETER, GLOBAL_CURRENT_CELL_OBJECT) {
  CELL_OBJECT <- GLOBAL_CURRENT_CELL_OBJECT
  CELL_OBJECT <-
    CELL_OBJECT[, c(
      "position.x",
      "position.y",
      "position.z",
      "distance",
      "ATP",
      PARAMETER,
      "simulation_time"
    )]
  
  oldw <- getOption("warn")
  options(warn = -1)
  plot <- ggplot(data = CELL_OBJECT,
                 aes_string(x = "distance", y = CELL_OBJECT[[PARAMETER]])) +
    geom_point(alpha = 0.04, color = "#4575b4") +
    scale_color_viridis() +
    #ggtitle(PARAMETER, "/Distance") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 13)) +
    xlab("Distance (µm)") +
    ylab(PARAMETER) +
    #stat_summary(fun.data=mean_sd, geom="ribbon", alpha=0.25)+
    stat_smooth(color = "orange",
                geom = "smooth",
                level = 0.99) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),
          text = element_text(size = 21))
  plot <- ggplotly(plot) %>%
    layout(font = list(size = 22)) %>%
    config(displayModeBar = F) %>%
    toWebGL()
  options(warn = oldw)
  return(plot)
}
plotSpheroidAtTIME <-
  function(PARAMETER,
           GLOBAL_CURRENT_CELL_OBJECT,
           CHOSENMODE = 1,
           FORCE = FALSE) {
    MODE <- 1  #FAST MODE = 1 else = 2
    
    CELL_OBJECT <- GLOBAL_CURRENT_CELL_OBJECT
    CELL_OBJECT <- CELL_OBJECT[abs(CELL_OBJECT$position.z) <= 10,]
    column <- CELL_OBJECT[[PARAMETER]]
    discr <- FALSE
    if (PARAMETER == "current_phase") {
      column[column == 100] <- "Dead"#0
      column[column == 101] <- "Dead"#0
      column[column == 102] <- "Dead"#0
      column[column == 103] <- "Dead"#0
      column[column == 104] <- "Dead"#0
      column[column == 105] <- "Quiescent"#1
      column[column == 106] <- "Proliferating"#2
      column[column == 107] <- "Proliferating"#2
      column[column == 108] <- "Proliferating"#2
      column[column == 109] <- "Proliferating"#2
      column <- as.character(column)
      column <- as.factor(column)
      discr <- TRUE
      MODE <- 4
    }
    
    if (PARAMETER == "PHENOTYPE" | PARAMETER == "PHENOTYPE2") {
      column <- as.character(column)
      column <- as.factor(column)
      discr <- TRUE
      MODE <- 3
    }
    colorpal2 <- scale_colour_viridis(
      option = "plasma",
      direction = 1,
      discrete = discr,
      name = PARAMETER
    )
    colorpal <- scale_fill_viridis(
      option = "plasma",
      direction = 1,
      discrete = discr,
      name = PARAMETER
    )
    #  colorpal<-scale_fill_manual(values = c("#c15578","#f2f951","#0d0882"),labels=c("Proliferating","Quiescent","Dead"),drop=FALSE)
    
    if (PARAMETER != "PHENOTYPE" &
        PARAMETER != "PHENOTYPE2" & PARAMETER != "current_phase") {
      if (min(column) < 0) {
        #positive<-column[column>0]
        #  negative<-column[column<0]
        q1 <- min(column[!column %in% boxplot.stats(column)$out])
        q3 <- max(column[!column %in% boxplot.stats(column)$out])
        print(q1)
        print(q3)
        
        colorpal <- scale_fill_gradient2(
          low = "#f0150a",
          mid = "white",
          high = "#005bd1",
          midpoint = 0,
          limits = c(q1, q3),
          na.value = "lightgreen"
          #oob=squish
        )
        
      }
    }
    
    alphaDead <-
      1 - 0.8 * as.integer(as.logical(CELL_OBJECT$PHENOTYPE2 == "DEAD"))
    radius <- 0.7 * (3 * CELL_OBJECT$total_volume / (4 * pi)) ^ (1 / 3)
    # alpha<-pmax(0.01, 1-0.3*(CELL_OBJECT$Density)/max(CELL_OBJECT$Density))
    
    if (FORCE == TRUE) {
      MODE <- CHOSENMODE
    }
    
    
    if (MODE == 1) {
      #FOR SHINY
      ggplot(data = CELL_OBJECT, aes(x = position.x, y = position.y, z = column)) +
        geom_point(
          aes(x = position.x, y = position.y, colour = column),
          show.legend = TRUE,
          alpha = alphaDead,
          size = 1
        ) +
        #stat_summary_hex(drop = TRUE)+
        colorpal + colorpal2 +
        coord_equal() +
        labs(colour = "", x = NULL, y = NULL) +
        theme(legend.position = "right")
    }
    else if (MODE == 2) {
      #FOR THESIS
      ggplot() +
        geom_circle(
          aes(
            x0 = position.x,
            y0 = position.y,
            r = radius,
            fill = column,
            alpha = alphaDead
          ),
          n = 20,
          show.legend = TRUE,
          data = CELL_OBJECT
        ) +
        colorpal +
        coord_equal() +
        labs(colour = "", x = NULL, y = NULL) +
        theme(legend.position = "right")
    }
    else if (MODE == 3) {
      ggplot() +
        geom_point(
          aes(x = position.x, y = position.y, colour = column),
          show.legend = TRUE,
          data = CELL_OBJECT,
          alpha = alphaDead,
          size = 2
        ) +
        scale_color_manual(
          name = "Phenotype",
          values = c(
            "FERMENTATION (+) Respiration (-)" = "#F2375C",
            "FERMENTATION (++)" = "#C50027",
            "RESPIRATION (+) Fermentation (-)" = "#6D3BC4",
            "RESPIRATION (++)" = "#3D0D90",
            "REVERSE WARBURG" = "#29B7B2",
            "WARBURG MODERN" = "#FF7B12",
            "WARBURG ORIGINAL" = "#ba5200",
            "DEAD" = "#000000"
          )
        ) +
        coord_equal() +
        labs(colour = "", x = NULL, y = NULL) +
        theme(legend.position = "right")
    }
    else{
      ggplot() +
        geom_point(
          aes(x = position.x, y = position.y, colour = column),
          show.legend = TRUE,
          data = CELL_OBJECT,
          size = 0.7
        ) +
        scale_color_manual(
          name = "Current phase",
          values = c(
            "Quiescent" = "#f2f94f",
            "Proliferating" = "#c05377",
            "Dead" = "#080584"
          )
        ) +
        coord_equal() +
        labs(colour = "", x = NULL, y = NULL) +
        theme(legend.position = "right")
    }
  }
plotSpheroid3DAtTIME <-
  function(PARAMETER, GLOBAL_CURRENT_CELL_OBJECT) {
    CELL_OBJECT <- GLOBAL_CURRENT_CELL_OBJECT
    range01 <- function(x) {
      (x - min(x)) / (max(x) - min(x))
    }
    open3d(useNULL = TRUE)
    p <-
      spheres3d(
        x = CELL_OBJECT[, "position.x"],
        y = CELL_OBJECT[, "position.y"],
        z = CELL_OBJECT[, "position.z"],
        radius = ((3 * CELL_OBJECT[, "total_volume"]) / (4 * pi)) ^ (1 / 3),
        col = viridis(501, option = "C")[as.integer(range01(CELL_OBJECT[[PARAMETER]]) *
                                                      500) + 1],
        alpha = 1
      )
    # clipplanes3d(a = 0,b = 0,c = 1,d = 0)
    rglwidget()
  }
removeOutliers <- function(x) {
  return(x[!x %in% boxplot.stats(x)$out])
}
min_max_norm_data <- function(data, minus = F,outliers=T) {
  data_filtered<-data
  if (outliers == F) {
    data_filtered<-removeOutliers(data)
  }
  if (minus == T) {
    return(sign(data) * (abs(data) - min(abs(data_filtered))) / (max(data_filtered) - min(abs(data_filtered))))
  } else{
    return((data - min(data_filtered)) / (max(data_filtered) - min(data_filtered)))
  }
}
min_max_norm <- function(x, data) {
  (x - min(data)) / (max(data) - min(data))
}
plotsvg <- function(SUBSET, GLOBAL_CURRENT_CELL_OBJECT) {
  #fileXML <- read_xml(folder, )
  files.names <- dir(folder, pattern = "snapshot.*.svg")
  path_xml <- paste0(folder, files.names[TIME])
  # files.names <- dir(folder, pattern = "output.*.xml")
  # path <- paste0(folder, files.names[TIME])
  fileXML <- read_xml(path_xml)
  flux_list <-
    c(
      "GluT1",
      "HK",
      "GPI",
      "PFK1",
      "ALD",
      "ALD2",
      "TPI",
      "GAPDH",
      "PGK",
      "PGAM",
      "ENO",
      "PKM2",
      "LDH",
      "G6PD",
      "ATPases",
      "AK",
      "PFKFB",
      "PHGDH",
      "PDH",
      "PDH2",
      "ACC",
      "SOD",
      "SOD2",
      "SOD3",
      "NUCLEOTIDEBIOSYNTHESIS",
      "SERINECONSUMPTION",
      "GPDH",
      "GPDH2",
      "GPDH3",
      "NOX",
      "NOX2"
    )
  flux_reaction <-
    c(
      "r1",
      "r2",
      "r3",
      "r4",
      "r5",
      "r5",
      "r6",
      "r7",
      "r8",
      "r9",
      "r10",
      "r11",
      "r12",
      "r13",
      "r14",
      "r15",
      "r16",
      "r17",
      "r18",
      "r18",
      "r19",
      "r20",
      "r21",
      "r22",
      "r23",
      "r24",
      "r25",
      "r27",
      "r28",
      "r29",
      "r30"
    )
  WHOLE_OBJECT <- GLOBAL_CURRENT_CELL_OBJECT
  CELL_OBJECT <-
    GLOBAL_CURRENT_CELL_OBJECT %>% filter(eval(parse(text = SUBSET)))
  
  for (flux in seq(1:length(flux_list))) {
    path <- paste0("[id='", flux_list[flux])
    path <- paste0(path, "']")
    object <- CELL_OBJECT[, flux_reaction[flux]]
    object <- removeOutliers(object)
    group_median <- mean(object)
    
    # removing outliers
    COLLECTION <-
      (WHOLE_OBJECT[flux_reaction]) %>% select(-c(r20, r29, r14))
    # COLLECTION<-(WHOLE_OBJECT[flux_reaction])%>% select(c(r1))
    # COLLECTION<-colMeans(COLLECTION)
    COLLECTION <- as.vector(as.matrix(COLLECTION))
    # COLLECTION<-removeOutliers(COLLECTION)
    
    flux_scaled <- min_max_norm(group_median, COLLECTION)
    # flux_scaled=NanToZero((abs(mean_column)-mean(as.matrix(abs(column)))) / (max(as.matrix(abs(column)))-mean(as.matrix(abs(column)))))
    flux_sign <- sign(group_median)
    flux_value <- 0.2 + flux_scaled * 50
    col_positive <-
      colorRampPalette(c("white", "blue"))(255)[round(1 + (flux_scaled) * 255)]
    col_negative <-
      colorRampPalette(c("white", "red"))(255)[round(1 + (flux_scaled) * 255)]
    flux_color <- ifelse(flux_sign >= 0, col_positive, col_negative)
    
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 2) %>%
      xml_set_attr("stroke-width", toString(flux_value))
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 2) %>%
      xml_set_attr("stroke", flux_color)
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 3) %>%
      xml_set_attr("fill", flux_color)
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 4) %>%
      xml_set_attr("fill", flux_color)
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 3) %>%
      xml_set_attr("stroke", flux_color)
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 4) %>%
      xml_set_attr("stroke", flux_color)
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 3) %>%
      xml_set_attr("visibility", ifelse(flux_sign < 0, "visible", "hidden"))
    fileXML %>%
      xml_node(path) %>%
      xml_child(search = 4) %>%
      xml_set_attr("visibility", ifelse(flux_sign > 0, "visible", "hidden"))
  }
  rsvg_png(charToRaw(as.character(fileXML)), "out.png")
  HTML(fileXML %>% toString())
}
log10_ceiling <- function(x) {
  sign(x)*10^(ceiling(log10(abs(x))))
}
log10_floor <- function(x) {
  10^(floor(log10(x)))
}
plot2Variables <-
  function(SUBSET="Lactate>=0",
           GLOBAL_CURRENT_CELL_OBJECT,
           var1 = "LDH",
           var2 = "PDH",
           var3 = "density",outliers=T) {
    CELL_OBJECT <-
      GLOBAL_CURRENT_CELL_OBJECT %>% filter(eval(parse(text = SUBSET)))
    CELL_OBJECT <- CELL_OBJECT %>% filter(cycle_model == 9999)
    
    vvar1<-CELL_OBJECT[[var1]]
    vvar2<-CELL_OBJECT[[var2]]
    
    min_x <- min(vvar1)
    max_x <- max(vvar1)
    min_y <- min(vvar2)
    max_y <- max(vvar2)
    if (var1 == "LDH" && var2 == "PDH") {
      min_x <- 0
      max_x <- 12
      min_y <- 0
      max_y <- 1.2
    }
    if(outliers==F){
      min_x <- min(removeOutliers(vvar1))
      max_x <- max(removeOutliers(vvar1))
      min_y <- min(removeOutliers(vvar2))
      max_y <- max(removeOutliers(vvar2))
      #vvar1<-min_max_norm_data(CELL_OBJECT[[var1]],T,F)
      #vvar2<-min_max_norm_data(CELL_OBJECT[[var2]],T,F)
    }
    
    
    val2col <- function(val, min, max) {
      return(floor(255 * (val - min) / (max - min)))
    }
    xy2angle <- function(x, y) {
      angle <- atan2(y, x) * (180 / pi)
      if (angle < 0) {
        angle <- 360 + angle
      }
      return(angle)
    }
    val2hsl <- function(x, y, z) {
      hsl_list <- c()
      minx <- min(-1000)
      maxx <- max(1000)
      miny <- min(-1000)
      maxy <- max(1000)
      minz <- 0
      maxz <- max(CELL_OBJECT$distance)
      for (i in seq(1, length(x), by = 1)) {
        h <- xy2angle(x[i], y[i])
        s <- 1
        l <- 0.5 * (z[i] - minz) / (maxz - minz)
        hsl_list <- c(hsl_list, hsl2col(as.matrix(c(h, s, l))))
      }
      return(hsl_list)
    }
    
    
    ratio_values <- (max_x - min_x) / (max_y - min_y)
    # d <- ggplot(CELL_OBJECT, aes(LDH, PDH))
    # d <- d + stat_density2d(geom = "raster", bins=3000, aes(fill = after_stat(ndensity)), contour = FALSE, interpolate = FALSE) + scale_fill_viridis_c() + scale_y_continuous(limits = c(min_y, max_y), expand = c(0, 0)) + scale_x_continuous(limits = c(min_x, max_x), expand = c(0, 0)) + coord_fixed(ratio = ratio_values / 1.5)
    # d
    jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                    "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    if (var3 == "density") {
      d <- ggplot() + 
        stat_density_2d_filled(
          aes(vvar1, vvar2),
          bins=100, contour_var = "density",
          adjust=1,
          contour = TRUE
        )+ scale_fill_manual(values = rev(jet.colors(100)))+
        #scale_fill_viridis_c()+
        #geom_point(alpha=0.05,colour="white",size=1,shape=20)+
        scale_y_continuous(limits = c(min_y, max_y), expand = c(0, 0)) +
        scale_x_continuous(limits = c(min_x, max_x), expand = c(0, 0)) +
        theme_classic() + theme(panel.border = element_rect(
          colour = "black",
          fill = NA,
          size = 1
        )) + theme(legend.position = "none")+ coord_fixed(ratio = ratio_values / 1.5)+xlab(var1)+ylab(var2)
      
      return(d)
      
    }
    
    
    if (var3 == "pos") {
      CELL_OBJECT$color <-
        val2hsl(CELL_OBJECT$position.x,
                CELL_OBJECT$position.y,
                CELL_OBJECT$distance)
    } else{
      CELL_OBJECT$color <- CELL_OBJECT[[var3]]
    }
    hslmap <- ggplot(CELL_OBJECT, aes(position.x, position.y)) +
      #stat_density2d_filled(contour = TRUE, bins=20, contour_var = "ndensity") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_continuous(expand = c(0, 0)) +
      coord_fixed(ratio = 1) + theme_void() + theme(legend.position = "none")
    
    d <- ggplot(CELL_OBJECT, aes(.data[[var1]], .data[[var2]])) +
      #stat_density2d_filled(contour = TRUE, bins=20, contour_var = "ndensity") +
      scale_y_continuous(limits = c(min_y, max_y), expand = c(0, 0)) +
      scale_x_continuous(limits = c(min_x, max_x), expand = c(0, 0)) +
      theme_classic() + theme(panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1
      ))
    
    if (var1 == "LDH" && var2 == "PDH") {
      hslmap <- hslmap + geom_point(alpha = 1,
                                    colour = CELL_OBJECT$color,
                                    size = 2)
      d <-
        d + geom_point(alpha = 0.1,
                       colour = CELL_OBJECT$color,
                       size = 2) + coord_fixed(ratio = ratio_values / 1.5)
    }
    if (var3 != "pos") {
      hslmap <-
        hslmap + geom_point(alpha = 1, aes(colour = color)) + scale_colour_viridis_c(option =
                                                                                       "plasma")
      d <-
        d + geom_point(alpha = 0.1, aes(colour = color), size = 2) + scale_colour_viridis_c(option =
                                                                                              "plasma")
    }
    
    return(ggarrange(
      d,
      hslmap,
      nrow = 1,
      ncol = 2,
      widths = c(5, 1)
    ))
  }
plot2VariablesAcrossTime <- function(seq,in_path,out_path,v1="LDH",v2="PDH",v3="density",keepoutliers=T) {
  for (i in seq){
    CELL_OBJECT<-GET_CELL_OBJECT(i,in_path)
    result<-plot2Variables("cycle_model == 9999",CELL_OBJECT,var1 = v1,var2=v2,var3=v3,outliers = keepoutliers)
    ggsave(
      filename = paste0(i,".jpg"),
      plot = result,
      device = "jpg",
      path = out_path,
      units = "px",
      bg = "white"
    )
  }
}
plotRegulationMatrix <- function(GLOBAL_CURRENT_CELL_OBJECT) {
  CELL_OBJECT <- GLOBAL_CURRENT_CELL_OBJECT
  CELL_FACTORS <- CELL_OBJECT %>%
    select(starts_with("GeneGamma")) %>%
    mutate_all(mean)
  FACTOR_MATRIX <-
    matrix(
      as.vector(colMeans(CELL_FACTORS)),
      nrow = length(factor_names),
      ncol = length(factor_names),
      byrow = TRUE
    )
  colnames(FACTOR_MATRIX) <- factor_names
  rownames(FACTOR_MATRIX) <- factor_names
  FACTOR_LIST <- as.data.frame(reshape2::melt(FACTOR_MATRIX))
  FACTOR_LIST <- FACTOR_LIST[FACTOR_LIST$value > 0.0,]
  FACTOR_LIST[FACTOR_LIST["value"] > 15.0, "value"] <- 15.0
  midpoint <-
    (1 - min(FACTOR_LIST$value)) / (max(FACTOR_LIST$value) - min(FACTOR_LIST$value))
  colorScale <-
    data.frame(z = c(zmin = 0, midpoint, zmax = 1),
               col = c("red", "white", "blue"))
  fig <-
    plot_ly(
      x = as.character(FACTOR_LIST$Var1),
      y = as.character(FACTOR_LIST$Var2),
      z = FACTOR_LIST$value,
      type = "heatmap",
      colorscale = colorScale
    ) %>% layout(plot_bgcolor = "black")
  fig
}
plotSpheroidAcrossTime <- function(PARAMETER, SEQUENCE) {
  for (current_time in SEQUENCE) {
    current_cell_object <- GET_CELL_OBJECT(current_time)
    p <- plotSpheroidAtTIME(PARAMETER, current_cell_object)
    filename <- paste0(paste0("plot_", current_time), ".jpg")
    ggsave(
      filename = filename,
      plot = p,
      device = "jpg",
      path = "./SHINY_APP/FIG_EXPORT"
    )
  }
}
### UI
ui <- function(request) {
  dashboardPage(
    dashboardHeader(
      title = "PhysiCell Simulations",
      tags$li(
        tags$style(
          ".navbar{height:65px!important; background-color:#00456d!important}"
        ),
        tags$style(".logo{height:65px!important}"),
        tags$style(
          ".main-header{height:65px!important; position:fixed!important; top:0px; width:100%}"
        ),
        tags$style(".irs-grid-text{color:white!important}"),
        tags$style(".irs-grid-pol{color:white!important}"),
        tags$style(".irs-line{color:white!important}"),
        tags$style(".irs-max{color:white!important}"),
        tags$style(".irs-min{color:white!important}"),
        fluidRow(
          column(
            1,
            checkboxInput("keep_on_last_step", "Keep at last step", FALSE)
          ),
          column(
            3,
            numericInput(
              "TimeStepVal",
              "Timestep:",
              min(GET_TIMESTEP()),
              min = min(GET_TIMESTEP()),
              max = max(GET_TIMESTEP())
            )
          ),
          column(
            8,
            sliderInput(
              inputId = "slider_time",
              label = NULL,
              min = min(GET_TIMESTEP()),
              max = max(GET_TIMESTEP()),
              step = 1,
              value = min(GET_TIMESTEP()),
              width = "100%"
            )
          ),
          style = "margin-right:0px;",
          bookmarkButton(),
          actionButton("LoadBookmarks", "Load", icon("refresh")),
          shinyDirButton("dir", "Input directory", "Upload")
        ),
        class = "dropdown",
        style = "height:50px;width:100%; position:absolute; left:7px; top:2px; color:white!important"
      )
    ),
    ## Sidebar content
    dashboardSidebar(
      collapsed = TRUE,
      disable = TRUE,
      sidebarMenu(
        style = "position:fixed;width: 220px;",
        menuItem(
          "Dashboard",
          tabName = "dashboard",
          icon = icon("dashboard")
        ),
        # menuItem("Widgets", tabName = "widgets", icon = icon("th"))
        # sliderInput("slider_time", "Timestep:", min=min(GET_TIMESTEP()), max=max(GET_TIMESTEP()), step=1, value=max(GET_TIMESTEP()))
        div("TIMELINE", style = "text-align:center")
      )
    ),
    ## Body content
    dashboardBody(# First tab content
      tabItem(
        tabName = "dashboard",
        fluidPage(
          tags$style(".content{padding-top:70px;}"),
          tags$style("#svgout svg{height:900px;width:100%;}"),
          column(7,
                 fluidRow(
                   box(
                     selectizeInput(
                       "Input_p1",
                       "Parameters",
                       NULL,
                       multiple = FALSE,
                       selected = "ID"
                     ),
                     plotlyOutput("plot1"),
                     width = 6
                   ),
                   box(
                     selectizeInput(
                       "Input_p2",
                       "Parameters",
                       NULL,
                       multiple = FALSE,
                       selected = "ID"
                     ),
                     plotOutput("plot2"),
                     width = 6
                   )
                 ),
                 fluidRow(
                   box(
                     selectizeInput(
                       "Input_p3",
                       "Parameters",
                       NULL,
                       multiple = FALSE,
                       selected = "ID"
                     ),
                     plotlyOutput("plot3"),
                     width = 6
                   ),
                   box(
                     selectizeInput(
                       "Input_p5",
                       "Parameters",
                       NULL,
                       multiple = FALSE,
                       selected = "ID"
                     ),
                     plotOutput("plot5"),
                     width = 6
                   )
                 )),
          column(5,
                 fluidRow(
                   box(
                     selectizeInput(
                       "Input_p4",
                       "Parameters",
                       NULL,
                       multiple = FALSE,
                       selected = "ID"
                     ),
                     rglwidgetOutput("plot4", height = "900px", width = "677px"),
                     width = 12
                   )
                 )),
          column(7,
                 fluidRow(
                   box(
                     actionButton("UpdateLandscape", "Update", icon("refresh")),
                     fluidRow(
                       column(
                         3,
                         textInput(
                           "selectcondition0",
                           "Condition",
                           value = "Lactate > 0",
                           width = NULL,
                           placeholder = "Write a logical condition"
                         )
                       ),
                       column(
                         3,
                         selectizeInput(
                           "Input_p8",
                           "Parameter 1",
                           NULL,
                           multiple = FALSE,
                           selected = "LDH"
                         )
                       ),
                       column(
                         3,
                         selectizeInput(
                           "Input_p9",
                           "Parameter 2",
                           NULL,
                           multiple = FALSE,
                           selected = "PDH"
                         )
                       ),
                       column(
                         3,
                         textInput(
                           "selectcondition1",
                           "Couleur",
                           value = "density",
                           width = NULL,
                           placeholder = "Write a logical condition"
                         )
                       )
                     ),
                     plotOutput("plot9", height = "900px"),
                     width = 12
                   )
                 )),
          column(5,
                 fluidRow(
                   box(
                     fluidRow(column(
                       3, actionButton("UpdateGraph", "Update", icon("refresh"))
                     ), column(
                       5,
                       numericInput(
                         "obs",
                         "Sampling each () steps:",
                         500,
                         min = 1,
                         max = 1000,
                         width = "200px"
                       )
                     )),
                     plotlyOutput("plot6", height = "900px"),
                     width = 12
                   )
                 )),
          column(12,
                 fluidRow(
                   box(
                     fluidRow(
                       column(3, actionButton(
                         "UpdateGraph2", "Update", icon("refresh")
                       )),
                       column(
                         4,
                         numericInput(
                           "obs2",
                           "Sampling each () steps:",
                           500,
                           min = 1,
                           max = 1000,
                           width = "200px"
                         )
                       ),
                       column(
                         4,
                         numericInput(
                           "levels",
                           "Levels:",
                           4,
                           min = 1,
                           max = 20,
                           width = "200px"
                         )
                       )
                     ),
                     selectizeInput(
                       "Input_p6",
                       "Parameters",
                       NULL,
                       multiple = FALSE,
                       selected = "ID"
                     ),
                     plotlyOutput("plot7", height = "900px"),
                     width = 12
                   )
                 )),
          column(12,
                 fluidRow(
                   box(
                     fluidRow(
                       column(3, actionButton(
                         "UpdateGraph3", "Update", icon("refresh")
                       )),
                       column(
                         4,
                         numericInput(
                           "obs3",
                           "Sampling each () steps:",
                           500,
                           min = 1,
                           max = 1000,
                           width = "200px"
                         )
                       ),
                       column(
                         4,
                         numericInput(
                           "levels2",
                           "Levels:",
                           4,
                           min = 1,
                           max = 20,
                           width = "200px"
                         )
                       )
                     ),
                     selectizeInput("Input_p7", "Parameters", NULL, multiple = FALSE),
                     plotlyOutput("plot8", height = "900px"),
                     width = 12
                   )
                 )),
          column(12,
                 fluidRow(
                   box(
                     actionButton("UpdateImg", "Update", icon("refresh")),
                     textInput(
                       "selectcondition",
                       "Condition",
                       value = "Lactate > 0",
                       width = NULL,
                       placeholder = "Write a logical condition"
                     ),
                     uiOutput("svgout"),
                     width = 12
                   )
                 ))
        )
      ))
  )
}
server <- function(input, output, session) {
  rootpath <- "C:/Users/StevenSummey/Documents/CancerSimulationPhysiCell/SHINY_APP/"
  shinyDirChoose(input, 'dir', roots = c(home = rootpath))
  Simdirectory <- reactive(input$dir)
  
  LoadSim("C:/Users/StevenSummey/Desktop/VB_Share/output/")
  
  # if (TRUE %in% grepl("^GeneGamma", cell_measured_data_name)) {
  #   cell_measured_data_name[grepl("^GeneGamma", cell_measured_data_name)] <-GammaNamesList
  #     #loadGammaNames(cell_measured_data_name[grepl("^GeneGamma", cell_measured_data_name)])
  #   cell_measured_data_name_filtered <<-
  #     cell_measured_data_name[!grepl("^GeneGamma", cell_measured_data_name)]
  # }
  # cell_measured_data_name <<-
  #   c(cell_measured_data_name, "PHENOTYPE", "PHENOTYPE2")
  
  updateInputLocal <- function() {
    updateSelectizeInput(session,
                         "Input_p1",
                         choices = envir_data,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p2",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p3",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p4",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p5",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p6",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p7",
                         choices = envir_data,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p8",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
    updateSelectizeInput(session,
                         "Input_p9",
                         choices = cell_measured_data_name_filtered,
                         server = TRUE)
  }
  
  
  options(rgl.useNULL = TRUE)
  save <- options(rgl.inShiny = TRUE)
  on.exit(options(save))
  updateInputLocal()
  autoInvalidate <- reactiveTimer(5000, session)
  
  
  
  observeEvent(ignoreNULL = TRUE,
               eventExpr = {
                 input$dir
               },
               handlerExpr = {
                 req(is.list(input$dir))
                 home <- normalizePath(rootpath)
                 path <-
                   paste0(file.path(home, paste(
                     unlist(Simdirectory()$path[-1]), collapse = .Platform$file.sep
                   )), "/")
                 print(path)
                 LoadSim(path)
                 # if (TRUE %in% grepl("^GeneGamma", cell_measured_data_name)) {
                 #   cell_measured_data_name[grepl("^GeneGamma", cell_measured_data_name)] <-GammaNamesList
                 #     #loadGammaNames(cell_measured_data_name[grepl("^GeneGamma", cell_measured_data_name)])
                 #   cell_measured_data_name_filtered <<-
                 #     cell_measured_data_name[!grepl("^GeneGamma", cell_measured_data_name)]
                 # }
                 # cell_measured_data_name <<-
                 #   c(cell_measured_data_name, "PHENOTYPE", "PHENOTYPE2")
                 updateSliderInput(session,
                                   "slider_time",
                                   value = 0,
                                   max = max(GET_TIMESTEP()))
                 updateInputLocal()
               })
  
  
  onBookmark(function(state) {
    saved_state <- sub(".*/", "", state$dir)
    files.names <- dir("shiny_bookmarks")
    print(saved_state)
    for (files in files.names) {
      path <- paste0("shiny_bookmarks/", files)
      if (!grepl(saved_state, files, fixed = TRUE) &&
          files != "state.txt") {
        unlink(path, recursive = TRUE)
      }
    }
    write.table(
      saved_state,
      "shiny_bookmarks/state.txt",
      append = FALSE,
      col.names = FALSE,
      row.names = FALSE
    )
  })
  
  observeEvent(input$LoadBookmarks, {
    RDS_path <-
      paste0("shiny_bookmarks/", as.character(read.table("shiny_bookmarks/state.txt")$V1))
    state <- readRDS(paste0(RDS_path, "/input.rds"))
    updateSelectizeInput(
      session,
      "Input_p1",
      selected = state$Input_p1,
      choices = envir_data,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p2",
      selected = state$Input_p2,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p3",
      selected = state$Input_p3,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p4",
      selected = state$Input_p4,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p5",
      selected = state$Input_p5,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p6",
      selected = state$Input_p6,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p7",
      selected = state$Input_p7,
      choices = envir_data,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p8",
      selected = state$Input_p8,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "Input_p9",
      selected = state$Input_p9,
      choices = cell_measured_data_name_filtered,
      server = TRUE
    )
  })
  
  REACTIVE_CELL_OBJECT <- reactive({
    GET_CELL_OBJECT(input$slider_time)
  })
  
  observeEvent({
    autoInvalidate()
    input$keep_on_last_step
    1
  },
  {
    if (input$keep_on_last_step == FALSE) {
      updateSliderInput(session,
                        "slider_time",
                        value = NULL,
                        max = max(GET_TIMESTEP()))
    } else {
      updateSliderInput(
        session,
        "slider_time",
        value = max(GET_TIMESTEP()),
        max = max(GET_TIMESTEP())
      )
    }
  },
  priority = 2)
  
  observeEvent(input$TimeStepVal,
               {
                 if ((as.numeric(input$TimeStepVal) != input$slider_time) &
                     input$TimeStepVal != "" & input$slider_time != "") {
                   updateSliderInput(session,
                                     "slider_time",
                                     value = input$TimeStepVal)
                 } else {
                   if (input$TimeStepVal == "") {
                     updateSliderInput(session, value = 0)
                   }
                 }
               },
               priority = 2)
  
  
  observeEvent(input$slider_time,
               {
                 if ((as.numeric(input$TimeStepVal) != input$slider_time) &
                     input$slider_time != "" & input$TimeStepVal != "") {
                   updateNumericInput(
                     session = session,
                     inputId = "TimeStepVal",
                     value = input$slider_time
                   )
                 }
               },
               priority = 2)
  
  
  
  
  output$plot1 <- renderPlotly({
    EnvironmentMAP(input$Input_p1,
                   input$slider_time,
                   REACTIVE_CELL_OBJECT())
  })
  output$plot2 <- renderPlot({
    plotSpheroidAtTIME(input$Input_p2, REACTIVE_CELL_OBJECT())
  })
  observeEvent(input$UpdateLandscape, {
    output$plot9 <- renderPlot({
      isolate({
        plot2Variables(
          input$selectcondition0,
          REACTIVE_CELL_OBJECT(),
          input$Input_p8,
          input$Input_p9,
          input$selectcondition1
        )
      })
    })
  })
  output$plot3 <- renderPlotly({
    plotCellatTime(input$Input_p3, REACTIVE_CELL_OBJECT())
  })
  output$plot4 <- renderRglwidget({
    # plotSpheroid3DAtTIME(input$Input_p4, REACTIVE_CELL_OBJECT())
  })
  
  output$plot5 <- renderPlot({
    plotSpheroidAtTIME(input$Input_p5, REACTIVE_CELL_OBJECT())
  })
  
  observeEvent(input$UpdateGraph, {
    output$plot6 <- renderPlotly({
      isolate({
        plotCellPopAcrossTime(seq(
          from = 0,
          to = max(GET_TIMESTEP()),
          by = input$obs
        ))
      })
    })
  })
  
  observeEvent(input$UpdateGraph2, {
    output$plot7 <- renderPlotly({
      isolate({
        PLOT_CELL_TIME(input$Input_p6, input$obs2, input$levels)
      })
    })
  })
  
  observeEvent(input$UpdateGraph3, {
    output$plot8 <- renderPlotly({
      isolate({
        PLOT_MEDIUM_TIME(input$Input_p7, input$obs3, input$levels2)
      })
    })
  })
  
  output$svgout <- renderUI({
    input$UpdateImg
    
    isolate({
      plotsvg(input$selectcondition, REACTIVE_CELL_OBJECT())
    })
  })
}
LoadSim <- function(path) {
  #folder <- "/stock/Physicell_data_backup/22-01-19/"
  folder <<- path
  #folder <-"/home/jacquepi/Téléchargements/Physicell sample/"
  #folder <- 'C:/Users/StevenSummey/Desktop/VB_Share/'
  ## GET FIRST XML
  data_param <<- read_xml(paste0(folder, "output00000000.xml"))
  measured_data <<-
    data_param %>%
    xml_find_all("//label") %>%
    xml_text()
  envir_data <<-
    data_param %>%
    xml_find_all("//variable") %>%
    xml_attr(attr = "name")
  if ("H" %in% envir_data == TRUE) {
    envir_data <<- c(envir_data, "pH")
  }
  
  ## RETRIEVE MEDIUM'S AND CELL'S PARAMETERS NAMES
  cell_measured_data_name <<- c()
  for (i in 1:length(measured_data)) {
    j <- 1
    vector_name <- c(".x", ".y", ".z")
    vector_size <-
      (
        data_param %>% xml_find_all("//label") %>% xml_attr(attr =
                                                              "size") %>% as.numeric()
      )[i]
    while (j <= vector_size) {
      if (vector_size > 1) {
        if (measured_data[i] == "position") {
          cell_measured_data_name <<-
            c(cell_measured_data_name,
              paste0(measured_data[i], vector_name[j]))
        } else {
          cell_measured_data_name <<-
            c(cell_measured_data_name,
              paste0(measured_data[i], seq(1:5000)[j]))
        }
      } else {
        cell_measured_data_name <<-
          c(cell_measured_data_name, measured_data[i])
      }
      j <- j + 1
    }
  }
  cell_measured_data_name <<-
    c(cell_measured_data_name, "PHENOTYPE", "PHENOTYPE2")
  if (TRUE %in% grepl("^GeneGamma", cell_measured_data_name)) {
    cell_measured_data_name[grepl("^GeneGamma", cell_measured_data_name)] <-GammaNamesList
    
    cell_measured_data_name_filtered <<-
      cell_measured_data_name[!grepl("^GeneGamma", cell_measured_data_name)]
  }
}
###############################################################################
# basepath <- "/stock/Physicell_data_backup/"
basepath <- folder
date <- "2024-03-24"
# path <- paste0(paste0(basepath, date), "/")
path <- basepath
object <- list()
index <- 0
jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
N_EVAL <- last(GET_TIMESTEP(path))
timestepsequence <-seq(0, last(GET_TIMESTEP(path)), floor(last(GET_TIMESTEP(path)) / (N_EVAL -1)))
name<-"lactate"

LoadSim(path)

for (i in timestepsequence) {
  #p<-EnvironmentMAPGG(name,i,GET_CELL_OBJECT(i,"/stock/Physicell_data_backup/22-06-22/"))+geom_contour(aes(z=oxygen),breaks =seq(3.0,7.5,0.1),color="white")
  #p<-p+ scale_fill_stepsn(name,colours = rev(jet.colors(1000)), n.breaks=50,limits = c(6, 7.5), oob=squish)
  p<-EnvironmentMAPGG(name,i,GET_CELL_OBJECT(i,path))+geom_contour(aes(z=lactate),color="white")
  p<-p+ scale_fill_stepsn(name,colours = rev(jet.colors(1000)), n.breaks=100,limits = c(0, 1), oob=squish,guide="colourbar")
  panel_height = unit(1,"npc") - sum(ggplotGrob(p)[["heights"]][-3]) - unit(1,"line")
  p<-p + guides(fill= guide_colorsteps(barheight=panel_height, even.steps = FALSE,label.theme=element_text(size=6)))+theme(legend.position = c(1.12,0))
  p
  ggsave(
    filename = paste0(
      paste0(paste0(paste0(
        paste0(date, "_"), name
      ), "_"), i),
      ".jpg"
    ),
    plot = p,
    device = "jpg",
    path = folder,
    units = "px",
    bg = "white"
  )
  index <- index + 1
}

for (i in timestepsequence) {
  object[[index + 1]] <- GET_CELL_OBJECT(i, path)
  index <- index + 1
}

for (name in c('PHENOTYPE',"PHENOTYPE2")) {
  plots <- list()
  index <- 0
  for (i in seq(0, (N_EVAL - 1))) {
    #object[[index + 1]] <- GET_CELL_OBJECT(i, path)
    p <- plotSpheroidAtTIME(name, as.data.frame(object[[index + 1]]))
    plots[[index + 1]] <- p + theme(legend.position = "none")
    index <- index + 1
  }
  
  result <- ggarrange(plotlist = plots,
                      nrow = 2,
                      ncol = 3)
  result
  ggsave(
    filename = paste0(
      paste0(paste0(paste0(
        paste0(date, "_"), name
      ), "_"), timestepsequence[2] - timestepsequence[1]),
      ".jpg"
    ),
    plot = result,
    device = "jpg",
    path = folder,
    units = "px",
    bg = "white"
  )
}




# number_of_living_cells=c()
# number_of_dead_cells=c()
# number_of_quiescent_cells=c()
# for( i in seq(from=1,to=400,by=50)){
# current<-GET_CELL_OBJECT(i)
# phase<-current$current_phase
# living<-phase[phase %in% c(106,107,108,109)]
# dead<-phase[phase %in% c(100,101,102,103,104)]
# quiescent<-phase[phase %in% c(105)]
# number_of_living_cells=c(number_of_living_cells,length(living))
# number_of_dead_cells=c(number_of_dead_cells,length(dead))
# number_of_quiescent_cells=c(number_of_quiescent_cells,length(quiescent))
# }
#
# ggplot(data=data.frame(),aes(x=seq(from=1,to=400,by=50),y=number_of_living_cells)) +
#   geom_line(aes(y = number_of_living_cells, colour = "number_of_living_cells")) +
#     geom_line(aes(y = number_of_dead_cells, colour = "number_of_dead_cells")) +
#   geom_line(aes(y = number_of_quiescent_cells, colour = "number_of_quiescent_cells"))
#
#
# library("gridExtra")
basepath <- folder
date <- "2024-03-24"
path <- paste0(paste0(basepath, date), "/")
object <- list()
index <- 0
N_EVAL <- 6
timestepsequence <-
  seq(0, last(GET_TIMESTEP(path)), floor(last(GET_TIMESTEP(path)) / (N_EVAL -
                                                                       1)))
for (i in timestepsequence) {
  object[[index + 1]] <- GET_CELL_OBJECT(i, path)
  index <- index + 1
}

for (name in c('PHENOTYPE',"PHENOTYPE2")) {
  plots <- list()
  index <- 0
  for (i in seq(0, (N_EVAL - 1))) {
    p <- plotSpheroidAtTIME(name, as.data.frame(object[[index + 1]]))
    plots[[index + 1]] <- p + theme(legend.position = "none")
    index <- index + 1
  }
  
  result <- ggarrange(plotlist = plots,
                      nrow = 2,
                      ncol = 3)
  result
  ggsave(
    filename = paste0(
      paste0(paste0(paste0(
        paste0(date, "_"), name
      ), "_"), timestepsequence[2] - timestepsequence[1]),
      ".jpg"
    ),
    plot = result,
    device = "jpg",
    path = folder,
    units = "px",
    bg = "white"
  )
}



library("gridExtra")


library(ggpubr)


meanEnvirOverTime <- function(PARAMETER) {
  mediums <- data.frame()
  index <- 1
  for (i in seq(from = 1,
                to = N_EVAL,
                by = 1)) {
    temp <- GET_MEDIUM_OBJECT(i, path)
    temp$TIME <- i
    mediums <- bind_rows(mediums, temp)
    index <- index + 1
  }
  MEDIUM_COPY = mediums
  MEDIUM_COPY <- MEDIUM_COPY %>% filter(Distance <= max(MEDIUM_COPY$x))
  MEDIUM_COPY$distance_group = cut(
    MEDIUM_COPY$Distance,
    breaks = seq(
      min(MEDIUM_COPY$Distance),
      max(MEDIUM_COPY$Distance) + min(MEDIUM_COPY$Distance),
      min(MEDIUM_COPY$Distance)
    ),
    include.lowest = TRUE,
    labels = FALSE
  ) * min(MEDIUM_COPY$Distance)
  
  
  MEDIUM_COPY <- MEDIUM_COPY %>%
    group_by(TIME) %>%
    group_by(distance_group, add = TRUE) %>%
    dplyr::summarize(Mean = mean(!!sym(PARAMETER) / (1e-15), na.rm = TRUE)) #Divided by volume of voxel
  MEDIUM_COPY$PARAMETER <- PARAMETER
  return(MEDIUM_COPY)
}

g <- meanEnvirOverTime("glucose")
o <- meanEnvirOverTime("oxygen")
p <- meanEnvirOverTime("pH")
l <- meanEnvirOverTime("lactate")

themeNoLegends1 <-
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 18),
  )
themeNoLegends2 <-
  theme(axis.title.y = element_blank(),
        text = element_text(size = 18),)

#ggaxis<-grid_arrange_shared_legend(ggplot(g,aes(x=TIME,y=distance_group,fill=Mean)) +geom_raster(interpolate = T) + scalePARAMETER +scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +themeNoLegends2)
ggg <-
  ggplot(g, aes(x = TIME, y = distance_group, fill = Mean)) + geom_raster(interpolate = F) + scale_fill_viridis_c(option = "plasma", "Glucose") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks =
                                                              pretty_breaks(n = 5)) + themeNoLegends1
ggo <-
  ggplot(o, aes(x = TIME, y = distance_group, fill = Mean)) + geom_raster(interpolate = F) + scale_fill_viridis_c(option = "plasma", "Oxygen ") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks =
                                                              pretty_breaks(n = 5)) + themeNoLegends1
ggp <-
  ggplot(p, aes(x = TIME, y = distance_group, fill = Mean)) + geom_raster(interpolate = F) + scale_fill_viridis_c(option = "plasma", "pH        ") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks =
                                                              pretty_breaks(n = 5)) + themeNoLegends1
ggl <-
  ggplot(l, aes(x = TIME, y = distance_group, fill = Mean)) + geom_raster(interpolate = F) + scale_fill_viridis_c(option = "plasma", "Lactate ") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks =
                                                              pretty_breaks(n = 5)) + themeNoLegends2

ggarrange(
  ggo,
  ggg,
  ggp,
  ggl,
  nrow = 4,
  common.legend = F,
  align = "none",
  heights = c(1, 1, 1, 1.25)
)

############
timepoint = 120 #hours
selected <- cell_measured_data_name_filtered[237:308]
current <- GET_CELL_OBJECT(timepoint,path)
reactions <- current %>% filter(cycle_model == 9999) %>% select(r1:r30)
columns <- colnames(reactions)
reactions <- as.data.frame(melt(reactions))
maxval <- summary(reactions$value)[6]
q3 <- quantile(reactions$value, 0.97)
q1 <- quantile(reactions$value, 0.001)
minval <- summary(reactions$value)[1]


ggplot(reactions,
       aes(
         x = variable,
         y = value,
         group = variable,
         fill = variable
       )) +
  geom_boxplot(outlier.shape = NA) + theme(legend.position = "none") +
  coord_cartesian(ylim = quantile(reactions$value, c(0.01, 0.99))) +
  labs(title = "Reactions",
       subtitle = "Quantiles | TIME: 120 Hours")

############ correlating multiple simulations (look at simlist)
library(GGally)
library(RColorBrewer)
#236:308 gene
 simlist <-
   c(
     "C:/Users/StevenSummey/Desktop/VB_Share/Apr_29_2024_0914/home/vboxuser/PhysiCell/output/", #daily media change ()
     "C:/Users/StevenSummey/Desktop/VB_Share/May_15_2024_0901/home/vboxuser/PhysiCell/output/", #no media change conditions (constant injection control)
     "C:/Users/StevenSummey/Desktop/VB_Share/May_17_2024_2249/home/vboxuser/PhysiCell/output/",
     "C:/Users/StevenSummey/Desktop/VB_Share/May_22_2024_0847/home/vboxuser/PhysiCell/output/",
     "C:/Users/StevenSummey/Desktop/VB_Share/May_28_2024_0903/home/vboxuser/PhysiCell/output/",#changed ATP threshold
     "C:/Users/StevenSummey/Desktop/VB_Share/May_29_2024_1519/home/vboxuser/PhysiCell/output/" #better condition changes
   )
 # "C:/Users/StevenSummey/Desktop/VB_Share/May_28_2024_0903/home/vboxuser/PhysiCell/output/"
 simlist <- unique(simlist)
#across <- GET_CELL_OBJECT(26, "/stock/Physicell_data_backup/22-06-27/")
across <- GET_CELL_OBJECT(24,path)
#across <- GET_ALL_SIM(simlist)
across <- GET_ALL_SIM(path)
#across <- GET_ALL_SIM_SAMPLE(simlist, 5000)
#across <- GET_ALL_SIM_SAMPLE(path,100)
across$metabolites_parameters118 <-
  -log10(across$metabolites_parameters118)
across <- across %>% filter(metabolites_parameters118 <= 8.0)

generateCorrelationMatrix <- function(object) {
  vect1 <- cell_measured_data_name_filtered[c(28:112)]
  vect2 <- cell_measured_data_name_filtered[c(117, 121, 122, 234)]
  finalmatrix <- list()
  index <- 1
  
  for (sim_ in unique(object$sim)) {
    sub_object <- object %>% filter(sim == sim_)
    matrixcorr <- expand.grid(vect1, vect2)
    matrixcorr$correlation <-
      apply(matrixcorr, 1, function(row)
        cor(sub_object[row[1]], sub_object[row[2]], method = "spearman"))

    #matrixcorr<-matrixcorr[matrixcorr$Var1!=matrixcorr$Var2,]
    matrixtrue <- acast(matrixcorr, Var1 ~ Var2, value.var = "correlation")
    colnames(matrixtrue)[colnames(matrixtrue) == 'metabolites_parameters1'] <-
      'GlucoseExtra'
    colnames(matrixtrue)[colnames(matrixtrue) == 'metabolites_parameters5'] <-
      'LactateExtra'
    colnames(matrixtrue)[colnames(matrixtrue) == 'metabolites_parameters6'] <-
      'OxygenExtra'
    rownames(matrixtrue)[rownames(matrixtrue) == 'metabolites_parameters1'] <-
      'GlucoseExtra'
    rownames(matrixtrue)[rownames(matrixtrue) == 'metabolites_parameters5'] <-
      'LactateExtra'
    rownames(matrixtrue)[rownames(matrixtrue) == 'metabolites_parameters6'] <-
      'OxygenExtra'
    colnames(matrixtrue)[colnames(matrixtrue) == 'metabolites_parameters118'] <-
      'pHExtra'
    rownames(matrixtrue)[rownames(matrixtrue) == 'metabolites_parameters118'] <-
      'pHExtra'
    
    #matrixtrue[is.na(matrixtrue)]<-0
    #matrixtrue <-matrixtrue[rowSums(abs(matrixtrue) >= 0.5, na.rm = TRUE) > 0,]
    print(nrow(matrixtrue))
    matrixtrue <-
      matrixtrue[rowSums(is.na(matrixtrue), na.rm = TRUE) == 0,]

    
    finalmatrix[[index]] <- matrixtrue
    index <- index + 1
  }
  return(finalmatrix)
}

finalmatrices<-generateCorrelationMatrix(across)
meanmatrice<-Reduce("+", finalmatrices) / length(finalmatrices)
chosen<-finalmatrices[[3]]
chosen<-par[,-nearZeroVar(par)]
chosen<-(cor(chosen%>%select(where(is.numeric)),method="spearman"))
chosen <-  t(chosen[rowSums(is.na(chosen), na.rm = TRUE) == 0,])
chosen <-  t(chosen[rowSums(is.na(chosen), na.rm = TRUE) == 0,])
chosen <-t(chosen[rowSums(abs(chosen) >= 0.15, na.rm = TRUE) > 1,])
chosen <-t(chosen[rowSums(abs(chosen) >= 0.15, na.rm = TRUE) > 1,])
corrplot(chosen, type="lower", method="circle",order = 'hclust', diag = T,tl.srt = 0.01,tl.cex = 0.6,tl.offset = 0.3, tl.col = "black" )


#####
cond_o2 <- c("anoxia", "hypoxia", "physioxia")
cond_glu <- c("low", "medium", "high")
cond_lac <- c("unfixed", "low", "high")
cond_pH <- c("unfixed", "low", "neutral")
conditions <- expand.grid(cond_o2, cond_glu, cond_lac, cond_pH)
conditions$viability <- "viable"
conditions[conditions$Var1 == "anoxia", "viability"] <- "NO"
length(conditions[conditions$viability == "viable", "viability"]) / length(conditions$viability)

######
matrixlist<-list()
par<-current%>% filter(cycle_model != 9999)
par<-par[,c("r12","r21","Pyruvate","Lactate","Glucose","metabolites_parameters6","metabolites_parameters118","metabolites_parameters1","metabolites_parameters5","PHENOTYPE","ATP","HIF","LDH","PDH","PHENOTYPE2")]
par<-par[,cell_measured_data_name_filtered[c(28:114,116:310)]]
par$metabolites_parameters118<--log10(par$metabolites_parameters118)
colnames(par)[colnames(par) == 'metabolites_parameters1'] <-
  'Glucose Extra'
colnames(par)[colnames(par) == 'metabolites_parameters5'] <-
  'Lactate Extra'
colnames(par)[colnames(par) == 'metabolites_parameters6'] <-
  'Oxygen Extra'
colnames(par)[colnames(par) == 'metabolites_parameters118'] <-
  'pH Extra'

colnames(par)[colnames(par) == 'r12'] <-
  'Flux LDH'
colnames(par)[colnames(par) == 'r21'] <-
  'Flux MCT'
# correlate_all<-cor(par2%>%select(where(is.numeric)),method="spearman")
# correlate_all_f <-t(correlate_all[rowSums(abs(correlate_all) >= 0.8, na.rm = TRUE) > 1,])
# correlate_all_f <-t(correlate_all_f[rowSums(abs(correlate_all_f) >= 0.8, na.rm = TRUE) > 1,])
# corrplot(correlate_all_f, type="lower", method="circle",order = 'hclust', diag = F,tl.srt = 90,tl.offset = 1, tl.col = "black" )

par2<-par[,-nearZeroVar(par)]%>%mutate(ID=row_number())
parmeta<-par2%>%select(ID,PHENOTYPE,PHENOTYPE2)
umpafit<-par2%>%select(where(is.numeric)) %>% column_to_rownames("ID") %>%scale()%>%umap()
umap_df <- umpafit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(parmeta, by="ID")%>%inner_join(par2%>%select(where(is.numeric)), by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2, 
             color = ATP))+
  geom_point()+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "UMAP plot")

#################
dev.off()
par(mfrow=c(1,1))
par(mfrow=c(7,11),mar=c(2,1,2,1))
colnames <- dimnames(par)[[2]]
for (i in 208:280) {
  #hist(current[,i], main=colnames[i], probability=TRUE, col="gray", border="white")
  hist(par[,i],30, main=colnames[i],col="white",xlab=NA, ylab=NA)
}
for (i in 1:65) {
  #hist(current[,i], main=colnames[i], probability=TRUE, col="gray", border="white")
plot(par[,"VEGF->SOD"],par[,i], main=colnames[i],xlab=NA, ylab=NA,pch=20, cex=0.2)
}


#################
kd <- with(current, MASS::kde2d(PDH, LDH, n = 500, lims = c(-0.1,1.2,-0.1,11)))

plot_ly(x = kd$y, y=kd$x, z =-kd$z, contours = list(
  z = list(show = TRUE, start = -0.3, end = 0, size = 0.05, color = 'white'))
  )%>%add_surface(colorscale='Jet')%>%
  layout(scene = list(xaxis = list(title = 'LDH',
  range=c(-0.1,11)
),yaxis = list(
  title = 'PDH',
  range=c(-0.1,1.2)
),zaxis=list(
  title="0    Density    1",
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
),
aspectmode='manual',aspectratio = list(x=1.4, y=1.4, z=0.5),camera = list(projection = list(type = 'orthographic'))))

library(misc3d)

current<-GET_CELL_OBJECT(11, folder)

plot_ly(current,x =~r12 , y=~r18, z =~Lactate, marker = list(color =~Lactate, colorscale = "Viridis", showscale = TRUE, size=1)) %>% add_markers()%>% layout(scene = list(xaxis = list(title = 'r12'),yaxis = list(title = 'r18'),zaxis = list(title = 'r11')))
# fig

################# BARPLOT
current<-GET_CELL_OBJECT(13,folder)
column<-current$PHENOTYPE2
column[column == 100] <- "Dead"#0
column[column == 101] <- "Dead"#0
column[column == 102] <- "Dead"#0
column[column == 103] <- "Dead"#0
column[column == 104] <- "Dead"#0
column[column == 105] <- "Quiescent"#1
column[column == 106] <- "Proliferating"#2
column[column == 107] <- "Proliferating"#2
column[column == 108] <- "Proliferating"#2
column[column == 109] <- "Proliferating"#2

column<-column [! column =="DEAD"]
column<-factor(c(column, "WARBURG MODERN","WARBURG ORIGINAL","REVERSE WARBURG"))

result<-ggplot(data.frame(column), aes(y=column,fill=column)) +
  geom_bar(aes(x = (..count..)/sum(..count..))) + 
  scale_x_continuous(labels=percent, limits = c(0,1))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(
  filename = "WarburgPH67.jpg",
  plot = result,
  device = "jpg",
  path = folder,
  units = "px",
  bg = "white",
  width = 1050,
  height = 300
)
#===============================================

  
  getCellPopAcrossTime <- function(sample,path=folder,PARAMETER) {
    living_cells <- c()
    quiescent_cells <- c()
    dead_cells <- c()
    all_cells <- c()
    for (value in 1:length(sample)) {
      all_cells_at_step <- GET_CELL_OBJECT(sample[value],path = path)$current_phase
      living_cell_at_step <-
        length(all_cells_at_step[(all_cells_at_step %in% c(106, 107, 108, 109))])
      living_cells <- c(living_cells, living_cell_at_step)
      
      quiescent_cells_at_step <-length(all_cells_at_step[(all_cells_at_step %in% c(105))])
      quiescent_cells<- c(quiescent_cells, quiescent_cells_at_step)
      
      all_cells <- c(all_cells, length(all_cells_at_step))
      dead_cells <-
        c(dead_cells, length(all_cells_at_step) - living_cell_at_step - quiescent_cells_at_step)
    }
    return(list("allcells"=all_cells,"livingcells"=living_cells,"quiescentcells"=quiescent_cells,"deadcells"=dead_cells))
  }

get_MEDIUM_TIME <- function(PARAMETER, SAMPLE, RADIUS) {
  across <-
    GET_ENV_ACROSS_TIME(PARAMETER, seq(
      from = 0,
      to = N_EVAL,
      by = SAMPLE
    ), RADIUS)

  return(across)
}

phenobourgeon<-getCellPopAcrossTime(seq(0,N_EVAL,1),path,"PHENOT")
folder<-path
mediumbourgeonO2<-get_MEDIUM_TIME("oxygen",1,1)%>%filter(time!=0)
mediumbourgeonglu<-get_MEDIUM_TIME("glucose",1,1)%>%filter(time!=0)
mediumbourgeonlac<-get_MEDIUM_TIME("lactate",1,1)%>%filter(time!=0)
mediumbourgeonpH<-get_MEDIUM_TIME("pH",1,1)%>%filter(time!=0)
tmpdf<-as.data.frame(phenobourgeon)
tmpdf$time<- 1:nrow(tmpdf)

df<-merge(tmpdf, mediumbourgeonO2, by = "time")                                 # Merge data frames by columns 


fig <- plot_ly(df)
fig <-
  fig %>% add_trace(
    x = ~ time,
    y = ~ livingcells,
    name = "proliferating cells",
    mode = "lines",
    type = "scatter"
  )%>% add_trace(
    x = ~ time,
    y = ~ quiescentcells,
    name = "quiescent cells",
    mode = "lines",
    type = "scatter"
  )%>% add_trace(
    x = ~ time,
    y = ~ livingcells+quiescentcells,
    name = "living cells",
    mode = "lines",
    line = list(color = 'blue', width = 2),
    type = "scatter"
  )%>% add_trace(
    x = ~ time,
    y = ~ deadcells,
    name = "dead cells",
    mode = "lines",
    type = "scatter"
  )%>% add_trace(
    x = ~ time,
    y = ~ allcells,
    name = "all cells",
    mode = "lines",
    type = "scatter"
  )



ay <- list(
  tickfont = list(color = "orange"),
  overlaying = "y",
  side = "right",
  title = "<b>oxygen level</b> "
)


fig <-
  fig %>% add_trace(
    x = ~ time,
    y = ~ mean*1e15,
    name = "oxygen level",
    yaxis = "y2",
    mode = "lines",
    line = list(color = 'orange', width = 2),
    type = "scatter"
  )


# Set figure title, x and y-axes titles

fig <- fig %>% layout(
  title = "Double Y Axis Example",
  yaxis2 = ay
) %>%
  
  layout(

    xaxis = list(
      zerolinecolor = '#ffff',
      
      zerolinewidth = 2,
      
      gridcolor = 'ffff'
    ),
    
    yaxis = list(
      zerolinecolor = '#ffff',
      
      zerolinewidth = 2,
      
      gridcolor = 'ffff'
    )
    
  )


fig%>%layout(xaxis=list(tickfont=list(size=22)),yaxis=list(tickfont=list(size=22)),yaxis2=list(tickfont=list(size=22)))



#ggplot(mediumbourgeonlac, aes(x = time, y = as.numeric(group), fill = mean*1e15)) + geom_raster(interpolate = F) + scale_fill_viridis_c(option = "plasma", "Lactate") +scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0), breaks =pretty_breaks(n = 5)) + themeNoLegends1
oo <- ggplot(mediumbourgeonO2, aes(time))+
  geom_ribbon(aes(ymin = q1*1e15, ymax = q3*1e15), fill="#009aab",alpha=0.2) +
  geom_line(aes(y = median*1e15))+themeNoLegends1
oo

gg <- ggplot(mediumbourgeonglu, aes(time))+
  geom_ribbon(aes(ymin = q1*1e15, ymax = q3*1e15), fill="#ffb300",alpha=0.2) +
  geom_line(aes(y = median*1e15))+themeNoLegends1

pp <- ggplot(mediumbourgeonpH, aes(time))+
  geom_ribbon(aes(ymin = q1*1e15, ymax = q3*1e15), fill="#2cb000",alpha=0.2) +
  geom_line(aes(y = median*1e15))+themeNoLegends1

ll <- ggplot(mediumbourgeonlac, aes(time))+
  geom_ribbon(aes(ymin = q1*1e15, ymax = q3*1e15), fill="#3600b3",alpha=0.2) +
  geom_line(aes(y = median*1e15))+themeNoLegends2

ggarrange(
  oo,
  gg,
  pp,
  ll,
  nrow = 4,
  common.legend = F,
  align = "v",
  heights = c(1, 1, 1, 1.25)
)

