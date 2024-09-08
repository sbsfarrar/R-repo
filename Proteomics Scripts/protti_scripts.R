install.packages("protti", dependencies = TRUE)
# install.packages("devtools")
devtools::install_github("jpquast/protti", dependencies = TRUE)

library(protti)
library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)

pd_prot_data <- read_protti("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_2024-0130_combined.csv")
pd_prot_data <- read_protti("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_2024-0130_combined peptides_ForAnalysis.csv")
#B <- read_protti("C:/Users/StevenSummey/OneDrive - Aruna Bio/Characterization/Proteomics/ArunaBio_2024-0130_02_Pre human percolator.csv")
# Select relevant columns
pd_prot_selected <- pd_prot_data %>%
  select(
    accession,
    description,
    gene_symbol,
    sum_pep_score,
    biological_process,
    molecular_function,
    cellular_component,
    kegg_pathways,
    reactome_pathways,
    coverage_percent,
    sequence,
    starts_with("number_peptides"),# select all columns that start with "number_peptides"
    starts_with("found_in"),# select all columns that start with "impute"
    starts_with("abundance_"),# select all columns that start with "impute"
    starts_with("score_mascot"),# select all columns that start with "impute"
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
    cols = starts_with("score_mascot"),
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
# Convert into long format
pd_prot_long4 <- pd_prot_data_filtered %>%
  pivot_longer(
    cols = starts_with("number_peptides_by"),
    names_to = "sample4",
    values_to = "number_peptides_by"
  )

colRemoveNames = c("found_in_sample_f1_sample_a1",
                   "found_in_sample_f2_sample_b1",
                   "found_in_sample_f3_sample_c1",
                   "score_mascot_f1_sample_a1",
                   "score_mascot_f2_sample_b1",
                   "score_mascot_f3_sample_c1",
                   "number_peptides_by_search_engine_f1_sample_a1",
                   "number_peptides_by_search_engine_f2_sample_b1",
                   "number_peptides_by_search_engine_f3_sample_c1"
                   )
for(i in 1:length(colRemoveNames)){
  longInd <- which(colnames(pd_prot_long) == colRemoveNames[i])
  pd_prot_long <- pd_prot_long[,-longInd]
}

protein_data_long <- cbind(pd_prot_long,
                           mascot=pd_prot_long2$mascot,
                           found_in=pd_prot_long3$found_in,
                           number_peptides_by=pd_prot_long4$number_peptides_by)
# Make an annotation data frame and merge it with your data frame to obtain conditions
# We are annotating sample 1-3 as controls and samples 4-6 as treated conditions

sample <- c( # make sure that the names are the same name as in your report
  "f1_sample_a1",
  "f2_sample_b1",
  "f3_sample_c1"
)
condition <- c(
  "Pre",
  "Ev",
  "Pep"
)

annotation <- data.frame(sample,condition)

# Combine your long data frame with the annotation
protein_groups_annotated <- protein_data_long %>% 
  left_join(y = annotation, by = "sample")
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#QUALITY CONTROL
input <- protein_groups_annotated
  # # as the data is log2 transformed, we need to transform it back before calculating the CVs
  # mutate(raw_intensity = 2^peptide_intensity_missing) 

qc_cvs(data = input,
       grouping = gene_symbol,
       condition = condition,
       intensity = abundance_raw, 
       plot = FALSE)
       #plot_style = "violin")

qc_ids(data = input,
       sample = sample,
       grouping = accession,
       intensity = abundance_raw,
       condition = condition, 
       plot = FALSE)


qc_ids(data = input,
       sample = sample,
       grouping = accession,
       intensity = abundance_raw,
       condition = condition, 
       title = "Protein identifications per sample",
       plot = TRUE)

qc_intensity_distribution(data = input,
                          sample = sample,
                          grouping = accession,
                          intensity = abundance_raw,
                          plot_style = "boxplot")


qc_sequence_coverage(data = input, 
                     protein_identifier = accession, 
                     coverage = coverage_percent)


qc_sample_correlation(data = input,
                      sample = sample, 
                      grouping = accession, 
                      intensity_log2 = abundance_raw, 
                      condition = condition, 
                      interactive = FALSE)

qc_pca(
  data = input,
  sample = sample,
  grouping = accession,
  intensity = abundance_raw,
  condition = condition,
  digestion = NULL,
  plot_style = "scree"
)

qc_pca(
  data = input,
  sample = sample,
  grouping = accession,
  intensity = abundance_raw,
  condition = condition,
  components = c("PC1", "PC2"), 
  plot_style = "pca"
)

# # Plot ranked peptide intensities
# qc_ranked_intensities(
#   data = input,
#   sample = sample,
#   grouping = accession,
#   intensity_log2 = abundance_raw,
#   plot = TRUE,
# )
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# SINGLE DOSE TREATMENT DATA ANALYSIS

data_normalised <- input %>%
  #filter(eg_is_decoy == FALSE) %>%
  mutate(intensity_log2 = log2(abundance_raw)) %>%
  normalise(sample = sample, 
            intensity_log2 = intensity_log2,
            method = "median")

# data_filtered <- data_normalised %>%
#   filter_cv(grouping = accession, 
#             condition = condition, 
#             log2_intensity = intensity_log2, 
#             cv_limit = 0.25,
#             min_conditions = 1)
uniprot_ids <- unique(data_normalised$accession)

uniprot <-
  fetch_uniprot(
    uniprot_ids = uniprot_ids,
    columns = c(
      "protein_name",
      "gene_names",
      "go_f",
      "xref_string", #was, xref_string
      "cc_interaction",
      "ft_act_site",
      "ft_binding",
      "xref_pdb",
      "length",
      "sequence"
    )
  ) %>% 
  #rename(pg_protein_accessions = accession)

data_filtered_uniprot <- data_normalised %>%
  left_join(y = uniprot,
            by = "accession") %>%
  # find_peptide(protein_sequence = sequence)%>%
  #              peptide_sequence = pep_stripped_sequence) 
  #              
  # assign_peptide_type(aa_before = aa_before,
  #                     last_aa = last_aa, 
  #                     aa_after = aa_after) %>%
  # calculate_sequence_coverage(protein_sequence = sequence)
  #                             peptides = pep_stripped_sequence)

diff_abundance_data <- data_filtered_uniprot %>%
  assign_missingness(
    sample = sample,
    condition = condition,
    grouping = accession,
    intensity = normalised_intensity_log2,
    ref_condition = "Pre",
    completeness_MAR = 0.7,
    completeness_MNAR = 0.25,
    retain_columns = c(accession, 
                       go_f, 
                       xref_string, 
                       start, 
                       end, 
                       length, 
                       coverage_percent)
  ) %>%
  calculate_diff_abundance(
    sample = sample,
    condition = condition,
    grouping = accession,
    intensity_log2 = normalised_intensity_log2,
    missingness = missingness,
    comparison = comparison,
    method = "moderated_t-test",
    retain_columns = c(accession, 
                       go_f, 
                       xref_string, 
                       start, 
                       end, 
                       length, 
                       coverage_percent)
  ) 








volcano_plot(
  data = diff_abundance_data,
  grouping = eg_precursor_id,
  log2FC = diff,
  significance = pval,
  method = "target",
  target_column = pg_protein_accessions,
  target = "P62942",
  x_axis_label = "log2(fold change) Rapamycin treated vs. untreated",
  significance_cutoff = c(0.05, "adj_pval") 
)