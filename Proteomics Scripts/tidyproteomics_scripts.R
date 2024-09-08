library(dplyr)
library(tidyproteomics)
#import to tidyproteomics
data_prot <- "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_protein_groups_filtered.xlsx" %>%
  import("ProteomeDiscoverer","proteins")# %>%
  # reassign('Pre', 'a7') %>%
  # reassign('EV', 'b7') %>%
  # reassign('Pep', 'c7')
data_prot <- data_prot %>%
  reassign(sample == 'a7', .replace = 'Pre') %>%
  reassign(sample == 'b7', .replace = 'EV') %>%
  reassign(sample == 'c7', .replace = 'Pep')

data_prot <- data_prot %>%
  # plot some simple summary stats
  plot_counts(destination = "png") %>%
  plot_quantrank(destination = "png") %>%
  plot_venn(destination = "png") %>%
  plot_euler(destination = "png")

data_prot <- data_prot %>%
  # save a table of simple summary stats
  summary("sample", destination = "save") %>%
  # save a report on contamination
  summary(contamination = "CRAP", destination = "save") %>%
  # remove contamination
  subset(!description %like% "^CRAP")

suppressWarnings(
data_prot <- data_prot %>%
  # normalize via several methods, best method will be automatically selected
  normalize(c("median","linear","limma","randomforest")) %>%
  # impute with a minimum value (this is a knock-out)
  impute(.function = base::min, method = 'row')) %>% #impute(base::mean)
 # plot visualizations comparing normalization methods
 # plot_normalization(destination = "png") %>%
 #   plot_variation_cv(destination = "png") %>%
 #   plot_variation_pca(destination = "png") %>%
 #   plot_dynamic_range(destination = "png") %>%
 #   # plot visualizations of unbiased clustering
  
  plot_heatmap(destination = "png") %>%
  plot_pca(destination = "png")

data_prot <- data_prot %>%
  # calculate the expression between experiment: ko and control: wt
  expression(EV/Pre) %>%
  # plot the expression analysis
  plot_volcano(EV/Pre, destination = "png", significance_column = "p_value") %>% 
  plot_proportion(EV/Pre, destination = "png")

data_prot <- data_prot %>%
  # calculate the enrichment of the GO term(s) using the results
  # from the expression analysis
  enrichment(EV/Pre, .term = "biological_process") %>%
  enrichment(EV/Pre, .term = "cellular_component") %>%
  enrichment(EV/Pre, .term = "molecular_function") %>%
  # plot the enrichment analysis
  plot_enrichment(EV/Pre, .term = "biological_process", destination = "png") %>%
  plot_enrichment(EV/Pre, .term = "cellular_component", destination = "png") %>%
  plot_enrichment(EV/Pre, .term = "molecular_function", destination = "png") %>%
  plot_enrichment(EV/Pre, .term = "biological_process", destination = "save") %>%
  plot_enrichment(EV/Pre, .term = "cellular_component", destination = "save") %>%
  plot_enrichment(EV/Pre, .term = "molecular_function", destination = "save")









data_prot %>% summary(by='sample')
data_prot %>% plot_quantrank()
data_prot %>% plot_pca()


data_prot %>%
  expression(Pre/EV) %>%
  enrichment(Pre/EV, .term = 'biological_process') %>%
  enrichment(Pre/EV, .term = 'molecular_function')