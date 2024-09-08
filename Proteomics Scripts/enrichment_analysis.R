library(ggplot2)
library(forcats)
library(tidyverse)
library(xml2)
library(rsvg)
library(svglite)

#plotting = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/EnrichmentAnalysis_PoC_CSV_filtered.csv")
#plotting = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/slimGO_BP.csv")
#plotting = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/slimGO_MF.csv")
plotting = read.csv("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/slimGO_CC.csv")



##GO SLIM STUFF -> append url prefix to idGO
prefix = "https://amigo.geneontology.org/amigo/term/"
temp = data.frame(Link=paste0(prefix,plotting$idGO))
plotting = cbind(plotting, temp)

graph = ggplot(data = plotting, aes(x = fct_rev(fct_reorder(Term, Percentages)), y = Percentages, color = Category)) + geom_col() + labs(x = NULL)
graph

plot <- ggplot(plotting, aes(x=fct_rev(fct_reorder(annotation,logP)),y = logP, fill=Comp, color = Category))
plot <- plot + geom_bar(stat = "identity", position = 'dodge')
plot <- plot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot

#slim GO BP
numUnique <- length(unique(plotting$Term))
plotSplit <- split(plotting, f=plotting$Term)

plotSplitA <- do.call(rbind.data.frame, plotSplit[c(seq(1:(numUnique/2)))])
plotSplitB <- do.call(rbind.data.frame, plotSplit[-c(seq(1:(numUnique/2)))])

plot <- ggplot(plotSplitA, aes(x=fct_rev(fct_reorder(Term,Percentages)),y = Percentages, fill=Sample))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', width = 0.6) + theme_dark() 
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("slimGO Terms") +
  labs(title = "Percent of Genes Binned to Annotated slimGO Terms, by Sample",
       subtitle = "Biological Processes")
plot

ggsave( tf1 <- tempfile(fileext=".svg"),plot)
links <- with(plotSplitA, setNames(Link,Term))
xml <-read_xml(tf1)
xml %>%
  xml_find_all(xpath="//d1:text") %>%
  keep(xml_text(.) %in% names(links)) %>%
  xml_add_parent("a", "xlink:href" = links[xml_text(.)], target="_blank")
write_xml(xml, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/BP_slimGO_1.svg")

write_xml(xml, tf2 <- tempfile(fileext = ".svg"))
rsvg_pdf(tf2, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/BP_slimGO_1.pdf")

plot <- ggplot(plotSplitB, aes(x=fct_rev(fct_reorder(Term,Percentages)),y = Percentages, fill=Sample))
plot <- plot + geom_bar(stat = "identity", position =  'dodge', width = 0.6) + theme_dark()
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("slimGO Terms") +
  labs(title = "Percent of Genes Binned to Annotated slimGO Terms, by Sample",
       subtitle = "Biological Processes")
plot

ggsave( tf1 <- tempfile(fileext=".svg"),plot)
links <- with(plotSplitB, setNames(Link,Term))
xml <-read_xml(tf1)
xml %>%
  xml_find_all(xpath="//d1:text") %>%
  keep(xml_text(.) %in% names(links)) %>%
  xml_add_parent("a", "xlink:href" = links[xml_text(.)], target="_blank")
write_xml(xml, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/BP_slimGO_2.svg")

write_xml(xml, tf2 <- tempfile(fileext = ".svg"))
rsvg_pdf(tf2, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/BP_slimGO_2.pdf")

#slim GO MF
numUnique <- length(unique(plotting$Term))

plot <- ggplot(plotting, aes(x=fct_rev(fct_reorder(Term,Percentages)),y = Percentages, fill=Sample))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', width = 0.6) + theme_dark() 
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("slimGO Terms") +
  labs(title = "Percent of Genes Binned to Annotated slimGO Terms, by Sample",
       subtitle = "Molecular Function")
plot

ggsave( tf1 <- tempfile(fileext=".svg"),plot)
links <- with(plotting, setNames(Link,Term))
xml <-read_xml(tf1)
xml %>%
  xml_find_all(xpath="//d1:text") %>%
  keep(xml_text(.) %in% names(links)) %>%
  xml_add_parent("a", "xlink:href" = links[xml_text(.)], target="_blank")
write_xml(xml, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/MF_slimGO.svg")

write_xml(xml, tf2 <- tempfile(fileext = ".svg"))
rsvg_pdf(tf2, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/MF_slimGO.pdf")

#slim GO CC
numUnique <- length(unique(plotting$Term))

plot <- ggplot(plotting, aes(x=fct_rev(fct_reorder(Term,Percentages)),y = Percentages, fill=Sample))
plot <- plot + geom_bar(stat = "identity", position = 'dodge', width = 0.6) + theme_dark() 
plot <- plot + theme(axis.title.x = element_blank(),
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("slimGO Terms") +
  labs(title = "Percent of Genes Binned to Annotated slimGO Terms, by Sample",
       subtitle = "Cellular Component")
plot

ggsave( tf1 <- tempfile(fileext=".svg"),plot)
links <- with(plotting, setNames(Link,Term))
xml <-read_xml(tf1)
xml %>%
  xml_find_all(xpath="//d1:text") %>%
  keep(xml_text(.) %in% names(links)) %>%
  xml_add_parent("a", "xlink:href" = links[xml_text(.)], target="_blank")
write_xml(xml, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/CC_slimGO.svg")

write_xml(xml, tf2 <- tempfile(fileext = ".svg"))
rsvg_pdf(tf2, "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/CC_slimGO.pdf")
