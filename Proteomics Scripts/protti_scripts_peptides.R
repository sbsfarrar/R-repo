install.packages("protti", dependencies = TRUE)
# install.packages("devtools")
devtools::install_github("jpquast/protti", dependencies = TRUE)

library(protti)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)

pd_pep_data <- read_protti("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_2024-0130_combined peptides_ForAnalysis.csv")

# Select relevant columns
pd_pep_selected <- pd_pep_data %>%
  select(
    pep_sequence,
    modifications,
    number_proteins,
    contaminant,
    master_protein_accessions,
    starts_with("abundance_"), # select all columns that start with "abundances_grouped"
    quan_info
  )

# Filter data frame
pd_pep_filtered <- pd_pep_selected %>%
  filter(contaminant == FALSE) %>% # remove annotated contaminants
  filter(number_proteins == 1) %>% # select proteotypic peptides
  filter(quan_info != "No Quan Values")# remove peptides that have no quantification values

# Convert into long format
pd_pep_long <- pd_pep_filtered %>%
  pivot_longer(
    cols = starts_with("abundance"),
    names_to = "sample",
    names_prefix = "abundance_",
    values_to = "abundance_raw"
  ) %>%
  # combine peptide sequence and modifications to make a precursor column
  mutate(precursor = paste(pep_sequence, modifications)) 


#####test to filter treatment group?
 # index = which(pd_pep_long$sample == "pep")
 # pd_pep_long <- pd_pep_long[-index,]

# Make an annotation data frame and merge it with your data frame to obtain conditions
# We are annotating sample 1-3 as controls and samples 4-6 as treated conditions

sample <- c( # make sure that the names are the same name as in your report
  "pre",
  "ev",
  "pep"
)
condition <- c(
  "control",
  "A",
  "B"
)

annotation <- data.frame(sample,condition)

# Combine your long data frame with the annotation
pd_pep_long_annotated <- pd_pep_long %>% 
  left_join(y = annotation, by = "sample")
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#QUALITY CONTROL (can do later)
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#SINGLE DOSE TREATMENT 
data_normalised <- pd_pep_long_annotated %>%
  mutate(intensity_log2 = log2(abundance_raw)) %>%
  normalise(sample = sample, 
            intensity_log2 = intensity_log2,
            method = "median")

 data_filtered <- data_normalised %>%
   filter_cv(grouping = master_protein_accessions, 
             condition = condition, 
             log2_intensity = intensity_log2, 
             cv_limit = 0.40,
             min_conditions = 1)

#Fetching database information and assigning peptide types
#data_filtered <- data_normalised
uniprot_ids <- unique(data_filtered$master_protein_accessions)

uniprot <-
  fetch_uniprot(
    uniprot_ids = uniprot_ids,
    columns = c(
      "protein_name",
      "gene_names",
      "go_f",
      "xref_string",
      "cc_interaction",
      "ft_act_site",
      "ft_binding",
      "xref_pdb",
      "length",
      "sequence"
    )
  ) %>% 
  rename(master_protein_accessions = accession)

data_filtered_uniprot <- data_filtered %>%
  left_join(y = uniprot,
            by = "master_protein_accessions") %>%
  find_peptide(protein_sequence = sequence,
               peptide_sequence = pep_sequence) %>%
  assign_peptide_type(aa_before = aa_before,
                      last_aa = last_aa, 
                      aa_after = aa_after) %>%
  calculate_sequence_coverage(protein_sequence = sequence,
                              peptides = pep_sequence)

qc_sequence_coverage(
  data = data_filtered_uniprot,
  protein_identifier = master_protein_accessions,
  coverage = coverage
)

#DATA ANALYSIS 

diff_abundance_data <- data_filtered_uniprot %>%
  assign_missingness(
    sample = sample,
    condition = condition,
    grouping = precursor,
    intensity = normalised_intensity_log2,
    ref_condition = "control",
    completeness_MAR = 0.7,
    completeness_MNAR = 0.25,
    retain_columns = c(master_protein_accessions, 
                       go_f, 
                       xref_string, 
                       start, 
                       end, 
                       length, 
                       coverage)
  ) %>%
  # impute(
  #   sample = sample,
  #   grouping = precursor,
  #   intensity_log2 = normalised_intensity_log2,
  #   condition = condition,
  #   comparison = comparison,
  #   missingness = missingness,
  #   noise = NULL,
  #   method = "ludovic",
  #   skip_log2_transform_error = FALSE,
  #   retain_columns = c(master_protein_accessions,
  #                      go_f,
  #                      xref_string,
  #                      start,
  #                      end,
  #                      length,
  #                      coverage)
  # ) %>%
  calculate_diff_abundance(
    sample = sample,
    condition = condition,
    grouping = precursor,
    intensity_log2 = normalised_intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "t-test",
    retain_columns = c(master_protein_accessions, 
                       go_f, 
                       xref_string, 
                       start, 
                       end, 
                       length, 
                       coverage)
  ) 

