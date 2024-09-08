library(tidyverse)
library(protti)
library(xlsx)

#must be csv 
#pd_prot_data <- read_protti(protPath_Pre <- "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_2024-0130_02_Pre human percolator.csv")
pd_prot_data <- read_protti("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_2024-0130_ human percolator consensus score.csv")
# Select relevant columns
pd_prot_selected <- pd_prot_data %>%
  select(
    master_protein_accessions,
    description,
    sum_pep_score,
    starts_with("number_peptides"),# select all columns that start with "number_peptides"
    starts_with("found_in"),# select all columns that start with "impute"
    starts_with("abundance_"),# select all columns that start with "impute"
    starts_with("mascot"),# select all columns that start with "impute"
  )

# Filter data frame
pd_prot_data_filtered <- pd_prot_selected %>%
  #filter(contaminant == FALSE) %>% # remove annotated contaminants
  filter(number_peptides > 1) # select proteins with more than one identified peptide

# Convert into long format
pd_prot_long <- pd_prot_data_filtered %>%
  pivot_longer(
    cols = starts_with("abundance"),
    names_to = "sample",
    names_prefix = "abundance_",
    values_to = "abundance_raw"
  )
# Convert into long format
pd_prot_long2 <- pd_prot_data_filtered %>%
  pivot_longer(
    cols = starts_with("mascot"),
    names_to = "sample2",
    values_to = "mascot"
  )
# Convert into long format
pd_prot_long3 <- pd_prot_data_filtered %>%
  pivot_longer(
    cols = starts_with("found_in"),
    names_to = "sample3",
    values_to = "found_in"
  )

colRemoveNames = c("found_in_sample_ev","found_in_sample_pre","found_in_sample_peptide",
  "mascot_ev","mascot_pre","mascot_peptide")
for(i in 1:length(colRemoveNames)){
  longInd <- which(colnames(pd_prot_long) == colRemoveNames[i])
  pd_prot_long <- pd_prot_long[,-longInd]
}

protein_data_long <- cbind(pd_prot_long,mascot=pd_prot_long2$mascot,found_in=pd_prot_long3$found_in)
# Make an annotation data frame and merge it with your data frame to obtain conditions
# We are annotating sample 1-3 as controls and samples 4-6 as treated conditions

sample <- c( # make sure that the names are the same name as in your report
  "pre",
  "ev",
  "peptide"
)
condition <- c(
  "control",
  "treated",
  "treated"
)

annotation <- data.frame(sample,condition)

# Combine your long data frame with the annotation
protein_groups_annotated <- protein_data_long %>% 
  left_join(y = annotation, by = "sample")

#export csv
write.csv(protein_groups_annotated,
          "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_protein_groups_annotated.csv",
          row.names=FALSE)
##
#export xlsx
write.xlsx(protein_groups_annotated,
          "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_protein_groups_annotated.xlsx",
          row.names=FALSE)

library(tidyproteomics)
#import to tidyproteomics
data_proteins <- "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_protein_groups_annotated.xlsx" %>%
  import("ProteomeDiscoverer","proteins")









write.xlsx(pd_prot_data_filtered,
           "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_protein_groups_filtered.xlsx")

library(tidyproteomics)
#import to tidyproteomics
data_proteins <- "C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_protein_groups_filtered_peptides.xlsx" %>%
  import("ProteomeDiscoverer","peptides")
