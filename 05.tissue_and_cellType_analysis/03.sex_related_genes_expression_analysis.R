# Stage-resolved version: only early/late spermatids are merged.
# Complete integrated script: tissue non-biased gains + original scRNA workflow.
# Run from the same project working directory as the original script.
# Eligibility: >=1 species with gonadal bias; non-biased counts as informative.
# Ancestry: unique majority, >=3 supporting species, >50% informative species.
# A single informative species is retained but cannot establish conservation.
# All genes remain in the QC heatmap/donut; no-bias-only genes are Other tissue.
# Low expression and missing data are not non-biased states.
# Majority-inferred ancestry is a heuristic, not phylogenetic reconstruction.

# Complete downstream analysis: majority-inferred tissue turnover + merged unresolved class
# Required existing objects:
#   df_summary: gene, species, E_T, E_O, category
#   my_genes: geneName plus Hab/Cin/Mbr/Dre ID columns

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

# New outputs are written into a dedicated folder.
output_dir <- "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------- settings ----------------
species_order <- c("Ham", "Hem", "Hha","Cha", "Dda",  "Mba")
species_sc_order <- c("Hab", "Cin", "Mbr", "Dre")
species_sc_labels <- c(Hab = "Ha", Cin = "Ch", Mbr = "Mb", Dre = "Dr")
min_biased_tissue_species <- 1
# Non-biased is informative; low_expr and missing_data are not.
min_valid_celltype_species <- 3L
min_ancestral_support_species <- 3

# Tissue enrichment cutoff on the log2 scale.
# 1 corresponds to two-fold enrichment.
gonad_enrichment_cutoff <- 1
min_gonad_TPM <- 1
sc_low_expr_quantile <- 0.50
sc_top1_top2_cutoff <- 1.5

sc_files <- c(
  Dre = "../seahorse_sc/mid_data/combined/pseudobulk/zbf_rmBias_pseudobulk.v2.tsv",
  Hab = "../seahorse_sc/mid_data/combined/pseudobulk/Hab_rmBias_pseudobulk.v2.tsv",
  Mbr = "../seahorse_sc/mid_data/combined/pseudobulk/Mbr_rmBias_pseudobulk.v2.tsv",
  Cin = "../seahorse_sc/mid_data/combined/pseudobulk/Cin_rmBias_pseudobulk.v3.tsv"
)

# ID-conversion files used in the original working scRNA workflow.
syngnathid_idconvert_file <-
  "../seahorse_sc/ref/ortho_zbfens99_Ha_Cin_Mbr/Ha_Cint_Mbr_zbfens99_singleCopy.convert3"

dre_idconvert_file <-
  "../ref/zbf/Lawson_zebrafish_v4_3_2geneinfo.txt"

# ---------------- tissue conserved/turnover ----------------

library(tidyverse)


my_genes <- readxl::read_xlsx(path = "mid_data/DNA_pipe/sex_related_genes/sex_related_genes6.xlsx", sheet = 1)
head(my_genes)
# my_genes <- my_genes %>% filter(if_singleCopy == 1)


my_genes <- my_genes %>%
  mutate(
    # For multi-copy genes, retain only the first copy
    across(
      c(Cin_ID, Mbr_ID, Hab_ID, Dre_ID, Dda_ID, Hhe_ID, Her_ID),
      ~ str_trim(str_split_fixed(as.character(.x), ",", 2)[, 1])
    ),
    
    # Remove species-specific suffixes
    Cin_ID = str_remove(Cin_ID, "_Cint"),
    Mbr_ID = str_remove(Mbr_ID, "_HJ"),
    Hab_ID = str_remove(Hab_ID, "_Ha"),
    Dre_ID = str_remove(Dre_ID, "__?zbfe99"),
    Dda_ID = str_remove(Dda_ID, "_Ddac"),
    Her_ID = str_remove(Her_ID, "_He"),
    
    # Remove transcript suffixes where necessary
    Hab_ID = str_split(Hab_ID, "\\.") %>% map_chr(1),
    Hhe_ID = str_split(Hhe_ID, "\\.") %>% map_chr(1)
  )

head(my_genes)
combined_df <- my_genes


# Cin
Cin_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/papa.Fo.tau.tsv")
head(Cin_df)
tmp_data <- Cin_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df) %>% str_remove("gene-") %>% str_split("_") %>% map_chr(1)
head(new_df)
table(combined_df$Cin_ID %in% new_df$geneID)
combined_df$Cin_ID[!combined_df$Cin_ID %in% new_df$geneID]
combined_df2 <- combined_df %>% left_join(y = new_df, by = c("Cin_ID" = "geneID"))
head(combined_df2)

# Dda
Dda_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/BJ.Fo.tau.tsv")
head(Dda_df)
tmp_data <- Dda_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df) %>% str_remove("gene-") %>% str_split("_") %>% map_chr(1)
head(new_df)
table(combined_df$Dda_ID %in% new_df$geneID)
combined_df$Dda_ID[!combined_df$Dda_ID %in% new_df$geneID]
combined_df3 <- combined_df2 %>% left_join(y = new_df, by = c("Dda_ID" = "geneID"))
head(combined_df3)


# Hhe
Hhe_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/260529/HF.Fo.tau.tsv")
head(Hhe_df)
tmp_data <- Hhe_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df) 
head(new_df)
table(combined_df$Hhe_ID %in% new_df$geneID)
combined_df$Hhe_ID[!combined_df$Hhe_ID %in% new_df$geneID]
combined_df4 <- combined_df3 %>% left_join(y = new_df, by = c("Hhe_ID" = "geneID"))

# Mbr
Mbr_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/HJ.Fo.tau.tsv")
head(Mbr_df)
tmp_data <- Mbr_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Mbr_ID %in% new_df$geneID)
combined_df5 <- combined_df4 %>% left_join(y = new_df, by = c("Mbr_ID" = "geneID"))

# Hab_adult
Hab_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/Ha_Hav.Fo.tau.tsv")
head(Hab_df)
tmp_data <- Hab_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Hab_ID %in% new_df$geneID)
combined_df6 <- combined_df5 %>% left_join(y = new_df, by = c("Hab_ID" = "geneID"))

# Hab, 3 month
Hab_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/Ha_Ha3.Fo.tau.tsv")
head(Hab_df)
tmp_data <- Hab_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Hab_ID %in% new_df$geneID)
combined_df6 <- combined_df6 %>% left_join(y = new_df, by = c("Hab_ID" = "geneID"))

# Hab, 2 month
Hab_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/Ha_Ha2.Fo.tau.tsv")
head(Hab_df)
tmp_data <- Hab_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Hab_ID %in% new_df$geneID)
combined_df6 <- combined_df6 %>% left_join(y = new_df, by = c("Hab_ID" = "geneID"))

# Hab, 1 month
Hab_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/Ha_Ha1.Fo.tau.tsv")
head(Hab_df)
tmp_data <- Hab_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Hab_ID %in% new_df$geneID)
combined_df6 <- combined_df6 %>% left_join(y = new_df, by = c("Hab_ID" = "geneID"))
colnames(combined_df6)



# Her, adult
Her_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/He_Hev.Fo.tau.tsv")
head(Her_df)
tmp_data <- Her_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Her_ID %in% new_df$geneID)
combined_df7 <- combined_df6 %>% left_join(y = new_df, by = c("Her_ID" = "geneID"))

# Her, 3 month
Her_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/He_He3.Fo.tau.tsv")
head(Her_df)
tmp_data <- Her_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Her_ID %in% new_df$geneID)
combined_df7 <- combined_df7 %>% left_join(y = new_df, by = c("Her_ID" = "geneID"))

# Her, 2 month
Her_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/He_He2.Fo.tau.tsv")
head(Her_df)
tmp_data <- Her_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df7$Her_ID %in% new_df$geneID)
combined_df7 <- combined_df7 %>% left_join(y = new_df, by = c("Her_ID" = "geneID"))

# Her, 1 month
Her_df <- read.table("mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/He_He1.Fo.tau.tsv")
head(Her_df)
tmp_data <- Her_df
tmp_data2 <- tmp_data %>% select(-c(contains("tau")))
head(tmp_data2)
new_df <- tmp_data2
new_df$geneID <- rownames(new_df)
head(new_df)
table(combined_df$Her_ID %in% new_df$geneID)
combined_df7 <- combined_df7 %>% left_join(y = new_df, by = c("Her_ID" = "geneID"))

colnames(combined_df7)
Ha_mean_df <- data.frame("Ham_Fa"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Fa" , value = T)),na.rm = T),
                         "Ham_Fb"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Fb" , value = T)),na.rm = T),
                         "Ham_Fk"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Fk" , value = T)),na.rm = T),
                         "Ham_Fo"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Fo" , value = T)),na.rm = T),
                         "Ham_Mp"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Mp" , value = T, ignore.case = T)),na.rm = T),
                         "Ham_Mb"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Mb" , value = T)),na.rm = T),
                         "Ham_Mk"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Mk" , value = T)),na.rm = T),
                         "Ham_Ms"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Ms" , value = T)),na.rm = T),
                         "Ham_Mt"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "Ha._Mt" , value = T)), na.rm = T)
)
rownames(Ha_mean_df) <- combined_df7$geneName
He_mean_df <- data.frame("Hem_Fa"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Fa" , value = T)),na.rm = T),
                         "Hem_Fb"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Fb" , value = T)),na.rm = T),
                         "Hem_Fk"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Fk" , value = T)),na.rm = T),
                         "Hem_Fo"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Fo" , value = T)),na.rm = T),
                         "Hem_Mp"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Mp" , value = T, ignore.case = T)),na.rm = T),
                         "Hem_Mb"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Mb" , value = T)),na.rm = T),
                         "Hem_Mk"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Mk" , value = T)),na.rm = T),
                         "Hem_Ms"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Ms" , value = T)),na.rm = T),
                         "Hem_Mt"=rowMeans(combined_df7 %>% 
                                             select(grep(x = colnames(combined_df7),pattern = "He._Mt" , value = T)), na.rm = T)
)

head(He_mean_df)

rownames(He_mean_df) <- combined_df7$geneName
combined_df9 <- cbind(Ha_mean_df, He_mean_df,
                      combined_df7 %>% select(-grep(x = colnames(combined_df7), pattern= "H[ae]._")) %>% select(-(1:9)))
combined_df9 %>% colnames()
combined_df9 <- combined_df9 %>% select(-grep(pattern = "_Mp|_MP|_Ms|_Fa",x = colnames(combined_df9),value = T))


# calculate T/S, O/S
df_long <- combined_df9 %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "TPM"
  ) %>%
  mutate(
    species = str_extract(sample, "^[^_]+"),
    
    tissue_group = case_when(
      str_detect(sample, "_Mt$") ~ "testis",
      str_detect(sample, "_Fo$") ~ "ovary",
      str_detect(sample, "_Mb$|_Mk$") ~ "male_somatic",
      str_detect(sample, "_Fb$|_Fk$") ~ "female_somatic",
      TRUE ~ NA_character_
    )
  )

df_summary <- df_long %>%
  filter(!is.na(tissue_group)) %>%
  group_by(gene, species, tissue_group) %>%
  summarise(
    TPM = if (all(is.na(TPM))) {
      NA_real_
    } else {
      max(TPM, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = tissue_group,
    values_from = TPM
  ) %>%
  rename(
    T   = testis,
    O   = ovary,
    S_M = male_somatic,
    S_F = female_somatic
  )


e <- 0.1

df_summary <- df_summary %>%
  mutate(
    E_T = log2((T + e) / (S_M + e)),
    E_O = log2((O + e) / (S_F + e))
  )


# Recalculate the tissue-expression category here so that the
# classification, turnover analysis and heatmap borders all use
# the same enrichment cutoff defined above.
df_summary <- df_summary %>%
  mutate(
    category = case_when(
      is.na(T) | is.na(O) |
        is.na(S_M) | is.na(S_F) ~
        "missing_data",
      
      pmax(T, O) <= min_gonad_TPM ~
        "low_expr",
      
      E_T > gonad_enrichment_cutoff &
        T > min_gonad_TPM &
        E_O > gonad_enrichment_cutoff &
        O > min_gonad_TPM ~
        "gonad_biased",
      
      E_T > gonad_enrichment_cutoff &
        T > min_gonad_TPM ~
        "testis_biased",
      
      E_O > gonad_enrichment_cutoff &
        O > min_gonad_TPM ~
        "ovary_biased",
      
      TRUE ~ "other"
    )
  )

df_state <- df_summary %>%
  filter(species %in% species_order) %>%
  mutate(
    # Four informative states: testis, ovary, both gonads, non-biased.
    # Low-expression and missing states do not establish absence of bias.
    tissue_state = case_when(
      category == "testis_biased" ~ "testis",
      category == "ovary_biased" ~ "ovary",
      category == "gonad_biased" ~ "dual",
      category == "other" ~ "non_bias",
      TRUE ~ NA_character_
    )
  )

if (anyDuplicated(df_state[c("gene", "species")]) > 0) {
  stop("Duplicated gene/species rows in df_state.")
}

df_gene_consensus <- df_state %>%
  filter(!is.na(tissue_state)) %>%
  group_by(gene) %>%
  summarise(
    n_species_used = n_distinct(species),
    n_biased_species = sum(tissue_state != "non_bias"),
    n_non_bias = sum(tissue_state == "non_bias"),
    biased_species = paste(species[tissue_state != "non_bias"], collapse = "; "),
    non_biased_species = paste(species[tissue_state == "non_bias"], collapse = "; "),
    tissue_states_observed = paste(sort(unique(tissue_state)), collapse = "; "),
    n_tissue_states = n_distinct(tissue_state),
    n_testis = sum(tissue_state == "testis"),
    n_ovary = sum(tissue_state == "ovary"),
    n_dual = sum(tissue_state == "dual"),
    state_counts = list(c(
      testis = sum(tissue_state == "testis"),
      ovary = sum(tissue_state == "ovary"),
      dual = sum(tissue_state == "dual"),
      non_bias = sum(tissue_state == "non_bias")
    )),
    observed_state_set = list(sort(unique(tissue_state))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    eligible_for_analysis = n_biased_species >= min_biased_tissue_species,
    maximum_state_count = max(state_counts),
    n_states_with_maximum = sum(state_counts == maximum_state_count),
    
    # An ancestral state is accepted only when it:
    #   1) occurs in at least three species;
    #   2) occurs in more than half of the valid species; and
    #   3) is the unique most frequent state.
    ancestral_state_is_supported =
      maximum_state_count >= min_ancestral_support_species &
      maximum_state_count > n_species_used / 2 &
      n_states_with_maximum == 1,
    
    inferred_ancestral_state = if (
      ancestral_state_is_supported
    ) {
      names(state_counts)[which.max(state_counts)]
    } else {
      NA_character_
    },
    
    derived_states = list(
      if (is.na(inferred_ancestral_state)) {
        character(0)
      } else {
        setdiff(observed_state_set, inferred_ancestral_state)
      }
    ),
    derived_states_text = if (
      length(derived_states) == 0
    ) {
      "none"
    } else {
      paste(derived_states, collapse = "; ")
    },
    n_derived_states = length(derived_states),
    
    final_category = case_when(
      !eligible_for_analysis ~ "no_bias_observed",
      n_species_used < 2 ~ "insufficient_data",
      n_species_used >= 4 & n_tissue_states == 1 & n_testis > 0 ~ "conserved_testis_biased",
      n_species_used >= 4 & n_tissue_states == 1 & n_ovary > 0 ~ "conserved_ovary_biased",
      n_species_used >= 4 & n_tissue_states == 1 & n_dual > 0 ~ "conserved_gonad_biased",
      n_tissue_states == 1 ~ "insufficient_data",

      # Unsupported ancestry remains unresolved. Retain multiple gains
      # separately when the supported ancestral state is non-biased.
      !ancestral_state_is_supported ~ "unresolved_turnover",
      inferred_ancestral_state == "non_bias" & n_derived_states > 1 ~
        "non_bias_to_multiple_turnover",
      n_derived_states > 1 ~ "unresolved_turnover",
      inferred_ancestral_state == "non_bias" & derived_states_text == "testis" ~ "non_bias_to_testis_turnover",
      inferred_ancestral_state == "non_bias" & derived_states_text == "ovary" ~ "non_bias_to_ovary_turnover",
      inferred_ancestral_state == "non_bias" & derived_states_text == "dual" ~ "non_bias_to_dual_turnover",
      inferred_ancestral_state == "testis" & derived_states_text == "non_bias" ~ "testis_to_non_bias_turnover",
      inferred_ancestral_state == "ovary" & derived_states_text == "non_bias" ~ "ovary_to_non_bias_turnover",
      inferred_ancestral_state == "dual" & derived_states_text == "non_bias" ~ "dual_to_non_bias_turnover",
      inferred_ancestral_state == "testis" & derived_states_text == "ovary" ~ "testis_to_ovary_turnover",
      inferred_ancestral_state == "testis" & derived_states_text == "dual" ~ "testis_to_dual_turnover",
      inferred_ancestral_state == "ovary" & derived_states_text == "testis" ~ "ovary_to_testis_turnover",
      inferred_ancestral_state == "ovary" & derived_states_text == "dual" ~ "ovary_to_dual_turnover",
      inferred_ancestral_state == "dual" & derived_states_text == "testis" ~ "dual_to_testis_turnover",
      inferred_ancestral_state == "dual" & derived_states_text == "ovary" ~ "dual_to_ovary_turnover",
      TRUE ~ "insufficient_data"
    )
  ) %>%
  ungroup() %>%
  select(-state_counts, -observed_state_set, -derived_states)

# One authoritative gene-level tissue classification used by the heatmap,
# donut and subtype bar plot. This prevents later summary figures from
# reclassifying Other/Insufficient genes differently from the heatmap.
tissue_gene_display <- df_summary %>%
  group_by(gene) %>%
  summarise(
    n_other_tissue_species = sum(category == "other", na.rm = TRUE),
    n_low_expr_tissue_species = sum(category == "low_expr", na.rm = TRUE),
    n_missing_tissue_species = sum(category == "missing_data", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(df_gene_consensus, by = "gene") %>%
  mutate(
    final_category = replace_na(final_category, "insufficient_data"),
    n_biased_species = replace_na(n_biased_species, 0L),
    eligible_for_analysis = replace_na(eligible_for_analysis, FALSE),
    tissue_status = case_when(
      final_category %in% c(
        "conserved_testis_biased",
        "conserved_ovary_biased",
        "conserved_gonad_biased"
      ) ~ "Tissue conserved",
      str_detect(final_category, "_turnover$") ~ "Tissue turnover",
      final_category == "no_bias_observed" ~ "Other tissue",
      TRUE ~ "Insufficient tissue data"
    ),
    tissue_subtype = case_when(
      final_category == "conserved_testis_biased" ~ "Conserved testis",
      final_category == "conserved_ovary_biased" ~ "Conserved ovary",
      final_category == "conserved_gonad_biased" ~ "Conserved both gonads",
      final_category == "testis_to_ovary_turnover" ~ "Testis to ovary",
      final_category == "testis_to_dual_turnover" ~ "Testis to both gonads",
      final_category == "ovary_to_testis_turnover" ~ "Ovary to testis",
      final_category == "ovary_to_dual_turnover" ~ "Ovary to both gonads",
      final_category == "dual_to_testis_turnover" ~ "Both gonads to testis",
      final_category == "dual_to_ovary_turnover" ~ "Both gonads to ovary",
      final_category == "non_bias_to_testis_turnover" ~ "Non-biased to testis",
      final_category == "non_bias_to_ovary_turnover" ~ "Non-biased to ovary",
      final_category == "non_bias_to_dual_turnover" ~ "Non-biased to both gonads",
      final_category == "testis_to_non_bias_turnover" ~ "Testis to non-biased",
      final_category == "ovary_to_non_bias_turnover" ~ "Ovary to non-biased",
      final_category == "dual_to_non_bias_turnover" ~ "Both gonads to non-biased",
      final_category == "non_bias_to_multiple_turnover" ~ "Non-biased to multiple biased states",
      final_category == "unresolved_turnover" ~ "Unresolved turnover",
      final_category == "no_bias_observed" ~ "Other tissue",
      TRUE ~ "Insufficient tissue data"
    ),
    plot_category = case_when(
      final_category == "conserved_ovary_biased" ~ "ovary-biased",
      final_category == "conserved_gonad_biased" ~ "gonad-biased",
      final_category == "conserved_testis_biased" ~ "testis-biased",
      final_category == "testis_to_ovary_turnover" ~ "testis-to-ovary turnover",
      final_category == "testis_to_dual_turnover" ~ "testis-to-dual turnover",
      final_category == "ovary_to_testis_turnover" ~ "ovary-to-testis turnover",
      final_category == "ovary_to_dual_turnover" ~ "ovary-to-dual turnover",
      final_category == "dual_to_testis_turnover" ~ "dual-to-testis turnover",
      final_category == "dual_to_ovary_turnover" ~ "dual-to-ovary turnover",
      final_category == "non_bias_to_testis_turnover" ~ "Non-biased to testis",
      final_category == "non_bias_to_ovary_turnover" ~ "Non-biased to ovary",
      final_category == "non_bias_to_dual_turnover" ~ "Non-biased to both gonads",
      final_category == "testis_to_non_bias_turnover" ~ "Testis to non-biased",
      final_category == "ovary_to_non_bias_turnover" ~ "Ovary to non-biased",
      final_category == "dual_to_non_bias_turnover" ~ "Both gonads to non-biased",
      final_category == "non_bias_to_multiple_turnover" ~ "Non-biased to multiple biased states",
      final_category == "unresolved_turnover" ~ "unresolved turnover",
      final_category == "no_bias_observed" ~ "Other tissue",
      TRUE ~ "Insufficient tissue data"
    )
  )


# ============================================================
# SECOND PASS: rescue ONLY genes originally classified as
# "Insufficient tissue data". All other genes are untouched.
# The first pass above retains the ORIGINAL gonad-only low_expr rule.
# ============================================================
tissue_first_pass_qc <- tissue_gene_display %>%
  select(gene, original_final_category = final_category,
         original_tissue_status = tissue_status)

rescue_gene_ids <- tissue_first_pass_qc %>%
  filter(original_tissue_status == "Insufficient tissue data") %>%
  pull(gene)

# Change low_expr -> other (= non_bias) only when somatic TPM > 1,
# and ONLY in genes originally classified as insufficient.
# Missing values are not rescued.
rescue_detail <- df_summary %>%
  filter(gene %in% rescue_gene_ids,
         category == "low_expr",
         !is.na(S_M), !is.na(S_F),
         pmax(S_M, S_F) > min_gonad_TPM) %>%
  transmute(gene, species, T, O, S_M, S_F,
            original_state = "low_expr", revised_state = "non_bias")

df_summary <- df_summary %>%
  mutate(category = if_else(
    gene %in% rescue_gene_ids & category == "low_expr" &
      !is.na(S_M) & !is.na(S_F) &
      pmax(S_M, S_F) > min_gonad_TPM,
    "other", category
  ))

# Recompute tissue classifications using the ORIGINAL rules.
# Because only originally insufficient genes changed their per-species
# state, all other genes retain precisely their original results.
df_state <- df_summary %>%
  filter(species %in% species_order) %>%
  mutate(
    # Four informative states: testis, ovary, both gonads, non-biased.
    # Low-expression and missing states do not establish absence of bias.
    tissue_state = case_when(
      category == "testis_biased" ~ "testis",
      category == "ovary_biased" ~ "ovary",
      category == "gonad_biased" ~ "dual",
      category == "other" ~ "non_bias",
      TRUE ~ NA_character_
    )
  )

if (anyDuplicated(df_state[c("gene", "species")]) > 0) {
  stop("Duplicated gene/species rows in df_state.")
}

df_gene_consensus <- df_state %>%
  filter(!is.na(tissue_state)) %>%
  group_by(gene) %>%
  summarise(
    n_species_used = n_distinct(species),
    n_biased_species = sum(tissue_state != "non_bias"),
    n_non_bias = sum(tissue_state == "non_bias"),
    biased_species = paste(species[tissue_state != "non_bias"], collapse = "; "),
    non_biased_species = paste(species[tissue_state == "non_bias"], collapse = "; "),
    tissue_states_observed = paste(sort(unique(tissue_state)), collapse = "; "),
    n_tissue_states = n_distinct(tissue_state),
    n_testis = sum(tissue_state == "testis"),
    n_ovary = sum(tissue_state == "ovary"),
    n_dual = sum(tissue_state == "dual"),
    state_counts = list(c(
      testis = sum(tissue_state == "testis"),
      ovary = sum(tissue_state == "ovary"),
      dual = sum(tissue_state == "dual"),
      non_bias = sum(tissue_state == "non_bias")
    )),
    observed_state_set = list(sort(unique(tissue_state))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    eligible_for_analysis = n_biased_species >= min_biased_tissue_species,
    maximum_state_count = max(state_counts),
    n_states_with_maximum = sum(state_counts == maximum_state_count),
    
    # An ancestral state is accepted only when it:
    #   1) occurs in at least three species;
    #   2) occurs in more than half of the valid species; and
    #   3) is the unique most frequent state.
    ancestral_state_is_supported =
      maximum_state_count >= min_ancestral_support_species &
      maximum_state_count > n_species_used / 2 &
      n_states_with_maximum == 1,
    
    inferred_ancestral_state = if (
      ancestral_state_is_supported
    ) {
      names(state_counts)[which.max(state_counts)]
    } else {
      NA_character_
    },
    
    derived_states = list(
      if (is.na(inferred_ancestral_state)) {
        character(0)
      } else {
        setdiff(observed_state_set, inferred_ancestral_state)
      }
    ),
    derived_states_text = if (
      length(derived_states) == 0
    ) {
      "none"
    } else {
      paste(derived_states, collapse = "; ")
    },
    n_derived_states = length(derived_states),
    
    final_category = case_when(
      !eligible_for_analysis ~ "no_bias_observed",
      n_species_used < 2 ~ "insufficient_data",
      n_species_used >= 4 & n_tissue_states == 1 & n_testis > 0 ~ "conserved_testis_biased",
      n_species_used >= 4 & n_tissue_states == 1 & n_ovary > 0 ~ "conserved_ovary_biased",
      n_species_used >= 4 & n_tissue_states == 1 & n_dual > 0 ~ "conserved_gonad_biased",
      n_tissue_states == 1 ~ "insufficient_data",
      # Unsupported ancestry remains unresolved. Retain multiple gains
      # separately when the supported ancestral state is non-biased.
      !ancestral_state_is_supported ~ "unresolved_turnover",
      inferred_ancestral_state == "non_bias" & n_derived_states > 1 ~
        "non_bias_to_multiple_turnover",
      n_derived_states > 1 ~ "unresolved_turnover",
      inferred_ancestral_state == "non_bias" & derived_states_text == "testis" ~ "non_bias_to_testis_turnover",
      inferred_ancestral_state == "non_bias" & derived_states_text == "ovary" ~ "non_bias_to_ovary_turnover",
      inferred_ancestral_state == "non_bias" & derived_states_text == "dual" ~ "non_bias_to_dual_turnover",
      inferred_ancestral_state == "testis" & derived_states_text == "non_bias" ~ "testis_to_non_bias_turnover",
      inferred_ancestral_state == "ovary" & derived_states_text == "non_bias" ~ "ovary_to_non_bias_turnover",
      inferred_ancestral_state == "dual" & derived_states_text == "non_bias" ~ "dual_to_non_bias_turnover",
      inferred_ancestral_state == "testis" & derived_states_text == "ovary" ~ "testis_to_ovary_turnover",
      inferred_ancestral_state == "testis" & derived_states_text == "dual" ~ "testis_to_dual_turnover",
      inferred_ancestral_state == "ovary" & derived_states_text == "testis" ~ "ovary_to_testis_turnover",
      inferred_ancestral_state == "ovary" & derived_states_text == "dual" ~ "ovary_to_dual_turnover",
      inferred_ancestral_state == "dual" & derived_states_text == "testis" ~ "dual_to_testis_turnover",
      inferred_ancestral_state == "dual" & derived_states_text == "ovary" ~ "dual_to_ovary_turnover",
      TRUE ~ "insufficient_data"
    )
  ) %>%
  ungroup() %>%
  select(-state_counts, -observed_state_set, -derived_states)

# One authoritative gene-level tissue classification used by the heatmap,
# donut and subtype bar plot. This prevents later summary figures from
# reclassifying Other/Insufficient genes differently from the heatmap.
tissue_gene_display <- df_summary %>%
  group_by(gene) %>%
  summarise(
    n_other_tissue_species = sum(category == "other", na.rm = TRUE),
    n_low_expr_tissue_species = sum(category == "low_expr", na.rm = TRUE),
    n_missing_tissue_species = sum(category == "missing_data", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(df_gene_consensus, by = "gene") %>%
  mutate(
    final_category = replace_na(final_category, "insufficient_data"),
    n_biased_species = replace_na(n_biased_species, 0L),
    eligible_for_analysis = replace_na(eligible_for_analysis, FALSE),
    tissue_status = case_when(
      final_category %in% c(
        "conserved_testis_biased",
        "conserved_ovary_biased",
        "conserved_gonad_biased"
      ) ~ "Tissue conserved",
      str_detect(final_category, "_turnover$") ~ "Tissue turnover",
      final_category == "no_bias_observed" ~ "Other tissue",
      TRUE ~ "Insufficient tissue data"
    ),
    tissue_subtype = case_when(
      final_category == "conserved_testis_biased" ~ "Conserved testis",
      final_category == "conserved_ovary_biased" ~ "Conserved ovary",
      final_category == "conserved_gonad_biased" ~ "Conserved both gonads",
      final_category == "testis_to_ovary_turnover" ~ "Testis to ovary",
      final_category == "testis_to_dual_turnover" ~ "Testis to both gonads",
      final_category == "ovary_to_testis_turnover" ~ "Ovary to testis",
      final_category == "ovary_to_dual_turnover" ~ "Ovary to both gonads",
      final_category == "dual_to_testis_turnover" ~ "Both gonads to testis",
      final_category == "dual_to_ovary_turnover" ~ "Both gonads to ovary",
      final_category == "non_bias_to_testis_turnover" ~ "Non-biased to testis",
      final_category == "non_bias_to_ovary_turnover" ~ "Non-biased to ovary",
      final_category == "non_bias_to_dual_turnover" ~ "Non-biased to both gonads",
      final_category == "testis_to_non_bias_turnover" ~ "Testis to non-biased",
      final_category == "ovary_to_non_bias_turnover" ~ "Ovary to non-biased",
      final_category == "dual_to_non_bias_turnover" ~ "Both gonads to non-biased",
      final_category == "non_bias_to_multiple_turnover" ~ "Non-biased to multiple biased states",
      final_category == "unresolved_turnover" ~ "Unresolved turnover",
      final_category == "no_bias_observed" ~ "Other tissue",
      TRUE ~ "Insufficient tissue data"
    ),
    plot_category = case_when(
      final_category == "conserved_ovary_biased" ~ "ovary-biased",
      final_category == "conserved_gonad_biased" ~ "gonad-biased",
      final_category == "conserved_testis_biased" ~ "testis-biased",
      final_category == "testis_to_ovary_turnover" ~ "testis-to-ovary turnover",
      final_category == "testis_to_dual_turnover" ~ "testis-to-dual turnover",
      final_category == "ovary_to_testis_turnover" ~ "ovary-to-testis turnover",
      final_category == "ovary_to_dual_turnover" ~ "ovary-to-dual turnover",
      final_category == "dual_to_testis_turnover" ~ "dual-to-testis turnover",
      final_category == "dual_to_ovary_turnover" ~ "dual-to-ovary turnover",
      final_category == "non_bias_to_testis_turnover" ~ "Non-biased to testis",
      final_category == "non_bias_to_ovary_turnover" ~ "Non-biased to ovary",
      final_category == "non_bias_to_dual_turnover" ~ "Non-biased to both gonads",
      final_category == "testis_to_non_bias_turnover" ~ "Testis to non-biased",
      final_category == "ovary_to_non_bias_turnover" ~ "Ovary to non-biased",
      final_category == "dual_to_non_bias_turnover" ~ "Both gonads to non-biased",
      final_category == "non_bias_to_multiple_turnover" ~ "Non-biased to multiple biased states",
      final_category == "unresolved_turnover" ~ "unresolved turnover",
      final_category == "no_bias_observed" ~ "Other tissue",
      TRUE ~ "Insufficient tissue data"
    )
  )

# Auditable before/after QC; verify all non-insufficient genes unchanged.
tissue_second_pass_qc <- tissue_gene_display %>%
  select(gene, revised_final_category = final_category,
         revised_tissue_status = tissue_status)
tissue_rescue_comparison <- tissue_first_pass_qc %>%
  full_join(tissue_second_pass_qc, by = "gene") %>%
  mutate(status_changed = original_tissue_status != revised_tissue_status |
           original_final_category != revised_final_category)

unexpected_changes <- tissue_rescue_comparison %>%
  filter(original_tissue_status != "Insufficient tissue data", status_changed)
if (nrow(unexpected_changes) > 0) {
  stop("Unexpected tissue classification changes outside original Insufficient tissue data genes")
}
readr::write_tsv(rescue_detail,
                 file.path(output_dir, "tissue_insufficient_rescued_species.tsv"))
readr::write_tsv(tissue_rescue_comparison,
                 file.path(output_dir, "tissue_insufficient_before_after.tsv"))
print(tissue_rescue_comparison %>%
        count(original_tissue_status, revised_tissue_status, name = "n_genes"))
# ============================================================

plot_category_levels <- c(
  "gonad-biased", "testis-biased", "ovary-biased",
  "testis-to-ovary turnover", "testis-to-dual turnover",
  "ovary-to-testis turnover", "ovary-to-dual turnover",
  "dual-to-testis turnover", "dual-to-ovary turnover",
  "Non-biased to testis",
  "Non-biased to ovary",
  "Non-biased to both gonads",
  "Testis to non-biased",
  "Ovary to non-biased",
  "Both gonads to non-biased",
  "Non-biased to multiple biased states",
  "unresolved turnover", "Other tissue", "Insufficient tissue data"
)

# All-gene QC retains genes with no informative species.
readr::write_tsv(tissue_gene_display, file.path(output_dir, "all_genes_tissue_QC.tsv"))
genes_with_any_bias <- tissue_gene_display %>%
  filter(eligible_for_analysis)
non_bias_gain_genes <- genes_with_any_bias %>%
  filter(inferred_ancestral_state == "non_bias")
non_bias_gain_details <- df_state %>%
  semi_join(non_bias_gain_genes, by = "gene") %>%
  left_join(non_bias_gain_genes %>%
              select(gene, inferred_ancestral_state, final_category), by = "gene") %>%
  arrange(gene, match(species, species_order))
readr::write_tsv(genes_with_any_bias,
                 file.path(output_dir, "tissue_genes_with_any_bias.tsv"))
readr::write_tsv(non_bias_gain_genes,
                 file.path(output_dir, "non_bias_ancestral_gain_genes.tsv"))
readr::write_tsv(non_bias_gain_details,
                 file.path(output_dir, "non_bias_ancestral_gain_details.tsv"))
print(genes_with_any_bias %>% count(final_category, name = "n_genes"))

# ---------------- read scRNA pseudobulk ----------------
missing_files <- sc_files[!file.exists(sc_files)]
if (length(missing_files) > 0) {
  stop("Missing scRNA pseudobulk file(s): ", paste(missing_files, collapse = ", "))
}

make_long_pb <- function(path, species_name) {
  read.table(path, sep = "\t", check.names = FALSE) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    pivot_longer(-gene_id, names_to = "celltype", values_to = "expr") %>%
    mutate(species = species_name, gene_id = as.character(gene_id), expr = as.numeric(expr))
}

pb_long <- map2_dfr(unname(sc_files), names(sc_files), make_long_pb)
# Correct the Hab gsdf feature ID
pb_long <- pb_long %>%
  mutate(
    gene_id = case_when(
      species == "Hab" &
        gene_id == "gsdf" ~ "HaiMa04881",
      TRUE ~ gene_id
    )
  )
# ---------------- original working ID conversion ----------------
# Hab/Cin/Mbr pseudobulk row names use LLgeneSymbol when an ortholog
# conversion is available. Dre pseudobulk row names also use
# LLgeneSymbol, obtained from the zebrafish gene-info table.

if (!file.exists(syngnathid_idconvert_file)) {
  stop("Missing ID-conversion file: ", syngnathid_idconvert_file)
}

if (!file.exists(dre_idconvert_file)) {
  stop("Missing Dre ID-conversion file: ", dre_idconvert_file)
}

IDconvert <- read.table(
  syngnathid_idconvert_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_idconvert_cols <- c("Mbr", "LLgeneSymbol")
if (!all(required_idconvert_cols %in% colnames(IDconvert))) {
  stop(
    "The syngnathid ID-conversion table must contain: ",
    paste(required_idconvert_cols, collapse = ", ")
  )
}

# Avoid column-name collisions if the upstream object was already used
# in a previous run of the analysis.
my_genes_sc <- my_genes %>%
  select(
    -any_of(
      c(
        "Hab_cID", "Cin_cID", "Mbr_cID", "Dre_cID",
        "LLgeneSymbol", "LLgeneSymbol.x", "LLgeneSymbol.y",
        "Dre_LLgeneSymbol"
      )
    )
  ) %>%
  left_join(
    IDconvert %>%
      select(Mbr, LLgeneSymbol) %>%
      distinct(),
    by = c("Mbr_ID" = "Mbr")
  ) %>%
  mutate(
    # This reproduces the original ifelse() logic:
    # use LLgeneSymbol when present; otherwise retain the species ID.
    Cin_cID = coalesce(
      as.character(LLgeneSymbol),
      as.character(Cin_ID)
    ),
    Mbr_cID = coalesce(
      as.character(LLgeneSymbol),
      as.character(Mbr_ID)
    ),
    Hab_cID = coalesce(
      as.character(LLgeneSymbol),
      as.character(Hab_ID)
    )
  )

Dre_IDconvert <- read.delim2(
  dre_idconvert_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  select(
    Ens99geneIDversion,
    LLgeneSymbol
  ) %>%
  mutate(
    Ens99geneIDversion = str_remove(
      as.character(Ens99geneIDversion),
      "\\.[0-9]+$"
    )
  ) %>%
  rename(
    Dre_LLgeneSymbol = LLgeneSymbol
  ) %>%
  distinct()

my_genes_sc <- my_genes_sc %>%
  mutate(
    Dre_join_ID = as.character(Dre_ID) %>%
      str_remove("_{1,2}zbfe99$") %>%
      str_remove("\\.[0-9]+$")
  ) %>%
  left_join(
    Dre_IDconvert,
    by = c("Dre_join_ID" = "Ens99geneIDversion")
  ) %>%
  mutate(
    Dre_cID = as.character(Dre_LLgeneSymbol)
  )

gene_map <- my_genes_sc %>%
  select(
    geneName,
    Hab_cID,
    Cin_cID,
    Mbr_cID,
    Dre_cID
  ) %>%
  pivot_longer(
    cols = c(Hab_cID, Cin_cID, Mbr_cID, Dre_cID),
    names_to = "species_col",
    values_to = "gene_id"
  ) %>%
  mutate(
    species = recode(
      species_col,
      Hab_cID = "Hab",
      Cin_cID = "Cin",
      Mbr_cID = "Mbr",
      Dre_cID = "Dre"
    ),
    gene_id = str_trim(as.character(gene_id))
  ) %>%
  select(geneName, species, gene_id) %>%
  filter(!is.na(geneName), geneName != "", !is.na(gene_id), gene_id != "") %>%
  distinct()

gene_map <- gene_map %>%
  mutate(
    gene_id = case_when(
      species == "Mbr" &
        geneName == "gsdf1" &
        gene_id == "g4136" ~ "gsdf",
      TRUE ~ gene_id
    )
  )
standardize_celltype <- function(x) {
  # Retain spermatogonia and spermatocyte; merge only early/late spermatids.
  x <- str_to_lower(str_trim(as.character(x)))
  case_when(
    x %in% c("early.spermatid", "late.spermatid") ~ "spermatid",
    TRUE ~ x
  )
}

pb_mapped <- pb_long %>%
  mutate(
    gene_id = str_trim(as.character(gene_id))
  ) %>%
  inner_join(gene_map, by = c("species", "gene_id")) %>%
  mutate(celltype = standardize_celltype(celltype)) %>%
  group_by(geneName, species, celltype) %>%
  summarise(
    expr = if (all(is.na(expr))) NA_real_ else max(expr, na.rm = TRUE),
    .groups = "drop"
  )

# Fail early if any scRNA species has no mapped target genes.
mapping_check <- pb_mapped %>%
  distinct(species, geneName) %>%
  count(species, name = "n_mapped_target_genes") %>%
  complete(
    species = species_sc_order,
    fill = list(n_mapped_target_genes = 0)
  )

print(mapping_check)

if (any(mapping_check$n_mapped_target_genes == 0)) {
  stop(
    "No mapped target genes for: ",
    paste(
      mapping_check$species[mapping_check$n_mapped_target_genes == 0],
      collapse = ", "
    )
  )
}

# ------------------------------------------------------------
# Species-specific low-expression threshold
#
# For each species:
#   1. For every gene in the complete pseudobulk dataset,
#      calculate its maximum expression across all cell types.
#   2. Use the median of these gene-wise maximum expression
#      values as the low-expression threshold.
# ------------------------------------------------------------

threshold_df <- pb_long %>%
  group_by(species, gene_id) %>%
  summarise(
    gene_max = if (all(is.na(expr))) {
      NA_real_
    } else {
      max(expr, na.rm = TRUE)
    },
    .groups = "drop"
  ) %>%
  group_by(species) %>%
  summarise(
    threshold = quantile(
      gene_max,
      probs = sc_low_expr_quantile,   # 0.50 = median
      na.rm = TRUE,
      names = FALSE
    ),
    n_genes_for_threshold = sum(!is.na(gene_max)),
    .groups = "drop"
  )

print(threshold_df)





get_major_celltypes <- function(celltype, expr) {
  keep <- !is.na(celltype) & !is.na(expr)
  celltype <- celltype[keep]
  expr <- expr[keep]
  if (length(expr) == 0) {
    return(tibble(major_celltype1 = NA_character_, major_celltype2 = NA_character_,
                  top1_expr = NA_real_, top2_expr = NA_real_, ratio_top1_top2 = NA_real_))
  }
  ord <- order(expr, decreasing = TRUE)
  celltype <- celltype[ord]
  expr <- expr[ord]
  top2_expr <- if (length(expr) >= 2) expr[2] else NA_real_
  ratio <- if (is.na(top2_expr)) Inf else if (top2_expr <= 0 && expr[1] > 0) Inf else if (top2_expr <= 0) NA_real_ else expr[1] / top2_expr
  tibble(
    major_celltype1 = celltype[1],
    major_celltype2 = if (length(celltype) >= 2) celltype[2] else NA_character_,
    top1_expr = expr[1], top2_expr = top2_expr, ratio_top1_top2 = ratio
  )
}

sc_ranked <- pb_mapped %>%
  group_by(geneName, species) %>%
  group_modify(~ get_major_celltypes(.x$celltype, .x$expr)) %>%
  ungroup() %>%
  left_join(threshold_df, by = "species") %>%
  mutate(
    top1_is_valid = !is.na(top1_expr) & !is.na(threshold) & top1_expr >= threshold,
    top2_is_valid = !is.na(top2_expr) & !is.na(threshold) & top2_expr >= threshold,
    top1_top2_are_close = !is.na(ratio_top1_top2) & ratio_top1_top2 <= sc_top1_top2_cutoff,
    celltype_class = case_when(
      is.na(top1_expr) | is.na(threshold) ~ "missing",
      !top1_is_valid ~ "low_expr",
      top2_is_valid & top1_top2_are_close & !is.na(major_celltype2) & major_celltype1 != major_celltype2 ~ "two_celltypes",
      TRUE ~ "one_celltype"
    ),
    major_celltype2 = if_else(celltype_class == "two_celltypes", major_celltype2, NA_character_)
  )

all_genes <- unique(as.character(df_summary$gene))
sc_states <- expand_grid(geneName = all_genes, species = species_sc_order) %>%
  left_join(sc_ranked, by = c("geneName", "species")) %>%
  mutate(
    celltype_class = replace_na(celltype_class, "missing"),
    plot_celltype1 = case_when(
      celltype_class == "missing" ~ "missing",
      celltype_class == "low_expr" ~ "low_expr",
      celltype_class %in% c("one_celltype", "two_celltypes") & !is.na(major_celltype1) ~ major_celltype1,
      TRUE ~ "missing"
    ),
    plot_celltype2 = if_else(celltype_class == "two_celltypes", major_celltype2, NA_character_)
  )

# ---------------- cell-type conserved/turnover ----------------
# Same criteria as the updated genome-wide analysis:
# >=3 of 4 species with a valid predominant cell type; "other" is informative.
# Conserved: at least one major cell type shared by all valid species.
# Turnover: no major cell type shared by all valid species.
# Low-expression and missing species are excluded from the comparison.
valid_celltype_states_parsimony <- sc_states %>%
  filter(celltype_class %in% c("one_celltype", "two_celltypes")) %>%
  rowwise() %>%
  mutate(celltype_set = list(sort(unique(na.omit(c(major_celltype1, major_celltype2)))))) %>%
  ungroup() %>%
  filter(lengths(celltype_set) > 0)

get_parsimonious_celltype_solution <- function(celltype_sets) {
  celltype_sets <- celltype_sets[lengths(celltype_sets) > 0]
  if (!length(celltype_sets)) {
    return(tibble(minimum_distinct_celltypes = 0L,
                  n_optimal_solutions = 0L,
                  parsimonious_celltypes_text = "none",
                  globally_shared_celltypes_text = "none",
                  ambiguous_parsimonious_solution = FALSE))
  }
  candidate_grid <- do.call(expand.grid, c(celltype_sets, stringsAsFactors = FALSE))
  solution_sets <- lapply(seq_len(nrow(candidate_grid)), function(i)
    sort(unique(as.character(unlist(candidate_grid[i, , drop = FALSE])))))
  minimum_distinct <- min(lengths(solution_sets))
  optimal_text <- sort(unique(vapply(
    solution_sets[lengths(solution_sets) == minimum_distinct],
    paste, collapse = " ↔ ", FUN.VALUE = character(1))))
  globally_shared <- Reduce(intersect, celltype_sets)
  tibble(minimum_distinct_celltypes = as.integer(minimum_distinct),
         n_optimal_solutions = length(optimal_text),
         parsimonious_celltypes_text = paste(optimal_text, collapse = " / "),
         globally_shared_celltypes_text = if (length(globally_shared))
           paste(sort(globally_shared), collapse = "; ") else "none",
         ambiguous_parsimonious_solution = length(optimal_text) > 1)
}

celltype_gene_status_parsimony <- valid_celltype_states_parsimony %>%
  group_by(geneName) %>%
  group_modify(~ get_parsimonious_celltype_solution(.x$celltype_set)) %>%
  ungroup() %>%
  left_join(valid_celltype_states_parsimony %>%
              group_by(geneName) %>%
              summarise(n_valid_celltype_species = n_distinct(species),
                        celltype_states_observed = paste(sort(unique(
                          map_chr(celltype_set, ~paste(sort(.x), collapse = "+")))), collapse = "; "),
                        .groups = "drop"), by = "geneName") %>%
  mutate(n_globally_shared_celltypes = if_else(
    globally_shared_celltypes_text == "none", 0L,
    str_count(globally_shared_celltypes_text, ";") + 1L),
    celltype_status = case_when(
      n_valid_celltype_species < min_valid_celltype_species ~ "Insufficient cell-type data",
      minimum_distinct_celltypes == 1 ~ "Conserved cell type",
      minimum_distinct_celltypes > 1 ~ "Cell-type turnover",
      TRUE ~ "Insufficient cell-type data")) %>%
  rename(gene = geneName)

celltype_gene_status <- tibble(gene = all_genes) %>%
  left_join(celltype_gene_status_parsimony, by = "gene") %>%
  mutate(n_valid_celltype_species = replace_na(n_valid_celltype_species, 0L),
         minimum_distinct_celltypes = replace_na(minimum_distinct_celltypes, 0L),
         n_optimal_solutions = replace_na(n_optimal_solutions, 0L),
         parsimonious_celltypes_text = replace_na(parsimonious_celltypes_text, "none"),
         globally_shared_celltypes_text = replace_na(globally_shared_celltypes_text, "none"),
         n_globally_shared_celltypes = replace_na(n_globally_shared_celltypes, 0L),
         celltype_states_observed = replace_na(celltype_states_observed, "none"),
         ambiguous_parsimonious_solution = replace_na(ambiguous_parsimonious_solution, FALSE),
         celltype_status = replace_na(celltype_status, "Insufficient cell-type data"))

readr::write_tsv(celltype_gene_status,
                 file.path(output_dir, "celltype_classification_min3.tsv"))

# ---------------- order genes ----------------
celltype_status_levels <- c(
  "Conserved cell type", "Cell-type turnover", "Insufficient cell-type data"
)

plot_df <- df_summary %>%
  left_join(tissue_gene_display, by = "gene") %>%
  left_join(celltype_gene_status, by = "gene") %>%
  mutate(
    celltype_status = replace_na(celltype_status, "Insufficient cell-type data")
  )

gene_order_df <- plot_df %>%
  group_by(gene, final_category, plot_category, celltype_status) %>%
  summarise(
    mean_ET = if (all(is.na(E_T))) NA_real_ else mean(E_T, na.rm = TRUE),
    mean_EO = if (all(is.na(E_O))) NA_real_ else mean(E_O, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    plot_category = factor(plot_category, levels = plot_category_levels),
    celltype_status = factor(celltype_status, levels = celltype_status_levels)
  ) %>%
  arrange(plot_category, celltype_status, mean_EO, mean_ET)

gene_order <- gene_order_df$gene

# Counts requested for each tissue category.
celltype_count_by_tissue <- gene_order_df %>%
  count(plot_category, celltype_status, name = "n_genes") %>%
  complete(plot_category, celltype_status, fill = list(n_genes = 0))
print(celltype_count_by_tissue)

# ---------------- tissue matrix ----------------
# 仅修改绘图数据，不改变原始表达值和分类结果
plot_wide <- plot_df %>%
  mutate(
    E_T = if_else(category == "low_expr", NA_real_, E_T),
    E_O = if_else(category == "low_expr", NA_real_, E_O),
    species = factor(species, levels = species_order)
  ) %>%
  select(gene, species, E_O, E_T) %>%
  pivot_longer(
    c(E_O, E_T),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric, E_O = "O", E_T = "T"),
    colname = paste(species, metric, sep = "_")
  ) %>%
  select(gene, colname, value) %>%
  pivot_wider(names_from = colname, values_from = value)
mat <- plot_wide %>% column_to_rownames("gene") %>% as.matrix()
expected_cols <- as.vector(rbind(paste0(species_order, "_O"), paste0(species_order, "_T")))
for (nm in setdiff(expected_cols, colnames(mat))) {
  missing_col <- matrix(
    NA_real_, nrow = nrow(mat), ncol = 1,
    dimnames = list(rownames(mat), nm)
  )
  mat <- cbind(mat, missing_col)
}
mat <- mat[gene_order, expected_cols, drop = FALSE]

add_gap_columns <- function(x, species_order) {
  out <- list()
  for (i in seq_along(species_order)) {
    sp <- species_order[i]
    out[[length(out) + 1]] <- x[, c(paste0(sp, "_O"), paste0(sp, "_T")), drop = FALSE]
    if (i < length(species_order)) {
      gap <- matrix(NA_real_, nrow(x), 1, dimnames = list(rownames(x), paste0("gap_", i)))
      out[[length(out) + 1]] <- gap
    }
  }
  do.call(cbind, out)
}
mat_gap <- add_gap_columns(mat, species_order)

column_labels <- colnames(mat_gap)
column_labels[grepl("_O$", column_labels)] <- "O"
column_labels[grepl("_T$", column_labels)] <- "T"
column_labels[grepl("^gap_", column_labels)] <- ""
species_block_with_gap <- unlist(map(seq_along(species_order), ~ c(species_order[.x], species_order[.x], if (.x < length(species_order)) "")))

row_info <- gene_order_df %>% column_to_rownames("gene")
row_info <- row_info[rownames(mat_gap), , drop = FALSE]
row_split <- factor(row_info$plot_category, levels = plot_category_levels)

# Two row-annotation columns shown on the far left.
# Tissue retains the detailed categories used to split/order the heatmap;
# Cell type shows conservation/turnover status independently.
tissue_annotation_col <- c(
  "ovary-biased" = "#0072B2",
  "gonad-biased" = "#7B3294",
  "testis-biased" = "#D73027",
  "testis-to-ovary turnover" = "#56B4E9",
  "testis-to-dual turnover" = "#F4A261",
  "ovary-to-testis turnover" = "#E76F51",
  "ovary-to-dual turnover" = "#5E60CE",
  "dual-to-testis turnover" = "#E69F00",
  "dual-to-ovary turnover" = "#4895EF",
  "Non-biased to testis" = "#009E73",
  "Non-biased to ovary" = "#66A61E",
  "Non-biased to both gonads" = "#1B9E77",
  "Testis to non-biased" = "#A6761D",
  "Ovary to non-biased" = "#D95F02",
  "Both gonads to non-biased" = "#7570B3",
  "Non-biased to multiple biased states" = "#006D2C",
  "unresolved turnover" = "#6A3D9A",
  "Other tissue" = "#9C755F",
  "Insufficient tissue data" = "#F0F0F0"
)

stopifnot(all(plot_category_levels %in% names(tissue_annotation_col)))

celltype_status_col <- c(
  "Conserved cell type" = "#009E73",
  "Cell-type turnover" = "#D55E00",
  "Insufficient cell-type data" = "#E5E5E5"
)

# Use only nonempty factor levels, in the exact heatmap slice order.
category_labels <- levels(droplevels(row_split))
category_label_gp <- gpar(fontsize = 9)
category_label_width <- ComplexHeatmap::max_text_width(
  category_labels, gp = category_label_gp
) + unit(5, "mm")
heatmap_row_gap <- unit(4, "mm")
heatmap_body_height <- nrow(mat_gap) * unit(3, "mm") +
  max(0, length(category_labels) - 1) * heatmap_row_gap

left_row_annotation <- rowAnnotation(
  Category = anno_block(
    labels = category_labels,
    labels_rot = 0,
    labels_gp = category_label_gp,
    labels_just = "right",
    labels_offset = unit(1, "npc") - unit(2, "mm"),
    gp = gpar(fill = NA, col = NA),
    width = category_label_width
  ),
  Tissue = as.character(row_info$plot_category),
  `Cell type` = as.character(row_info$celltype_status),
  col = list(
    Tissue = tissue_annotation_col,
    `Cell type` = celltype_status_col
  ),
  na_col = "white",
  simple_anno_size = unit(3.2, "mm"),
  gap = unit(1, "mm"),
  annotation_name_side = "bottom",
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 8),
  show_annotation_name = c(Category = FALSE, Tissue = TRUE, `Cell type` = TRUE)
)

# ---------------- scRNA matrices ----------------
sc_wide <- sc_states %>%
  filter(geneName %in% gene_order, species %in% species_sc_order) %>%
  select(geneName, species, plot_celltype1, plot_celltype2) %>%
  distinct() %>%
  mutate(species = factor(species, levels = species_sc_order)) %>%
  pivot_wider(names_from = species, values_from = c(plot_celltype1, plot_celltype2), names_sep = "__")

expected_sc_cols <- c(paste0("plot_celltype1__", species_sc_order), paste0("plot_celltype2__", species_sc_order))
for (nm in setdiff(expected_sc_cols, colnames(sc_wide))) sc_wide[[nm]] <- NA_character_
sc_wide <- sc_wide %>% column_to_rownames("geneName")
missing_genes <- setdiff(gene_order, rownames(sc_wide))
if (length(missing_genes) > 0) {
  extra <- matrix(NA_character_, length(missing_genes), ncol(sc_wide), dimnames = list(missing_genes, colnames(sc_wide)))
  sc_wide <- rbind(as.matrix(sc_wide), extra)
}
sc_wide <- sc_wide[gene_order, , drop = FALSE]
sc_mat1 <- as.matrix(sc_wide[, paste0("plot_celltype1__", species_sc_order), drop = FALSE])
sc_mat2 <- as.matrix(sc_wide[, paste0("plot_celltype2__", species_sc_order), drop = FALSE])
colnames(sc_mat1) <- unname(species_sc_labels[species_sc_order])
colnames(sc_mat2) <- unname(species_sc_labels[species_sc_order])
sc_mat1[is.na(sc_mat1) | sc_mat1 == ""] <- "missing"
sc_mat2[sc_mat2 == ""] <- NA_character_
stopifnot(identical(rownames(mat_gap), rownames(sc_mat1)))

# ---------------- colors and heatmaps ----------------
col_fun <- colorRamp2(c(-7, 0, 1, 7), c("#2166AC", "#F7F7F7", "#FDD49E", "#B2182B"))
celltype_col <- c(
  spermatogonia = "#0067a5", spermatocyte = "#f38400",
  spermatid = "#008856",
  sertoli = "#E68FAC", leydig = "#B3446C", immune = "#654522",
  erythroid = "#E9BD1B", other = "#BDBDBD", low_expr = "#E5E5E5", missing = "#FFFFFF"
)
get_celltype_color <- function(x) if (is.na(x) || !x %in% names(celltype_col)) "#FFFFFF" else unname(celltype_col[x])

species_topBar <- HeatmapAnnotation(
  species = anno_text(species_block_with_gap, gp = gpar(fontsize = 10, fontface = "bold")),
  annotation_height = unit(6, "mm")
)

# Build a matrix indicating gonadal enrichment
# using the updated tissue classifications

border_df <- df_summary %>%
  mutate(
    O_border = category %in% c(
      "ovary_biased", "gonad_biased"
    ),
    T_border = category %in% c(
      "testis_biased", "gonad_biased"
    )
  ) %>%
  select(gene, species, O_border, T_border) %>%
  pivot_longer(
    cols = c(O_border, T_border),
    names_to = "tissue",
    values_to = "border"
  ) %>%
  mutate(
    tissue = if_else(
      tissue == "O_border", "O", "T"
    ),
    colname = paste(species, tissue, sep = "_")
  ) %>%
  select(gene, colname, border) %>%
  pivot_wider(
    names_from = colname,
    values_from = border
  )

border_mat <- border_df %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Match the heatmap matrix, including gap columns
border_mat_gap <- matrix(
  FALSE,
  nrow = nrow(mat_gap),
  ncol = ncol(mat_gap),
  dimnames = dimnames(mat_gap)
)

common_genes <- intersect(
  rownames(mat_gap),
  rownames(border_mat)
)

common_cols <- intersect(
  colnames(mat_gap),
  colnames(border_mat)
)

border_mat_gap[
  common_genes,
  common_cols
] <- border_mat[
  common_genes,
  common_cols,
  drop = FALSE
]

stopifnot(
  identical(dimnames(border_mat_gap), dimnames(mat_gap))
)

# 与热图行列完全对应的低表达标记矩阵
low_expr_mat_gap <- matrix(
  FALSE,
  nrow = nrow(mat_gap),
  ncol = ncol(mat_gap),
  dimnames = dimnames(mat_gap)
)

low_expr_rows <- df_summary %>%
  filter(category == "low_expr") %>%
  distinct(gene, species)

for (k in seq_len(nrow(low_expr_rows))) {
  g <- as.character(low_expr_rows$gene[k])
  sp <- as.character(low_expr_rows$species[k])
  cols <- intersect(
    paste0(sp, c("_O", "_T")),
    colnames(low_expr_mat_gap)
  )
  
  if (g %in% rownames(low_expr_mat_gap) && length(cols) > 0) {
    low_expr_mat_gap[g, cols] <- TRUE
  }
}


tissue_ht <- Heatmap(
  mat_gap, name = "log2(gonad/soma)", col = col_fun, na_col = "white",
  width = ncol(mat_gap) * unit(6, "mm"), height = heatmap_body_height,
  cluster_rows = FALSE, cluster_columns = FALSE, row_split = row_split, row_gap = heatmap_row_gap,
  cluster_row_slices = FALSE, row_title = NULL,
  # Gene names are displayed once, on the right of the scRNA heatmap.
  show_row_names = FALSE,
  show_column_names = TRUE, column_labels = column_labels, column_names_gp = gpar(fontsize = 10),
  rect_gp = gpar(col = "grey90", lwd = 0.5), top_annotation = species_topBar,
  column_title = "Tissue RNA-seq: O = ovary/female soma; T = testis/male soma",
  cell_fun = function(j, i, x, y, width, height, fill) {
    
    # 低表达：覆盖为灰色
    if (isTRUE(low_expr_mat_gap[i, j])) {
      grid.rect(
        x, y, width, height,
        gp = gpar(
          fill = "#E5E5E5",
          col = "grey90",
          lwd = 0.5
        )
      )
      
    } else if (isTRUE(border_mat_gap[i, j])) {
      # 有效的性腺富集格子：保留黑色边框
      grid.rect(
        x, y, width, height,
        gp = gpar(
          fill = NA,
          col = "#333333",
          lwd = 1
        )
      )
    }
  }
)

sc_ht <- Heatmap(
  sc_mat1, name = "Major cell type", col = celltype_col, na_col = "white",
  rect_gp = gpar(type = "none"), width = ncol(sc_mat1) * unit(7, "mm"),
  height = nrow(sc_mat1) * unit(3, "mm"), cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = row_split, row_gap = heatmap_row_gap,
  cluster_row_slices = FALSE, row_title = NULL,
  show_row_names = TRUE, row_names_side = "right",
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  show_column_names = TRUE, column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title = "Single-cell RNA-seq",
  cell_fun = function(j, i, x, y, width, height, fill) {
    c1 <- as.character(sc_mat1[i, j]); c2 <- as.character(sc_mat2[i, j])
    f1 <- get_celltype_color(c1); f2 <- get_celltype_color(c2)
    if (is.na(c2) || c2 == "" || c2 == c1) {
      grid.rect(x, y, width, height, gp = gpar(fill = f1, col = "grey88", lwd = 0.5))
    } else {
      grid.rect(x - width / 4, y, width / 2, height, gp = gpar(fill = f1, col = NA))
      grid.rect(x + width / 4, y, width / 2, height, gp = gpar(fill = f2, col = NA))
      grid.rect(x, y, width, height, gp = gpar(fill = NA, col = "grey88", lwd = 0.5))
    }
  }
)

combined_ht <- left_row_annotation + tissue_ht + sc_ht
draw(
  combined_ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
  ht_gap = unit(c(2, 5), "mm")
)

dir.create("plots/DNA_pipe/sex_related_genes", recursive = TRUE, showWarnings = FALSE)
pdf_height <- max(10, nrow(mat_gap) * 3 / 25.4 +
                    max(0, length(category_labels) - 1) * 4 / 25.4 + 3)
pdf("plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/sex_related_genes.tissue_celltype.combined.synced.v14_nonbias_layout.cutoff1.pdf", width = 18, height = pdf_height)
draw(
  combined_ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
  ht_gap = unit(c(2, 5), "mm")
)
dev.off()

# Useful output tables
write.table(
  tissue_gene_display,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_gene_display.synced.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  df_gene_consensus,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_turnover_majority_unresolved_merged.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  celltype_gene_status,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/celltype_conserved_turnover.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  celltype_count_by_tissue,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/celltype_counts_by_tissue_category.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ================================================================
# FINAL FIGURE: TWO-RING DONUT + TISSUE-SUBTYPE BAR PLOT
# ================================================================

# Reuse the same authoritative tissue classification as the heatmap.
panel_gene_df <- tissue_gene_display %>%
  left_join(celltype_gene_status, by = "gene") %>%
  mutate(
    celltype_detail = case_when(
      celltype_status == "Conserved cell type" ~ "Cell-type conserved",
      celltype_status == "Cell-type turnover" ~ "Cell-type turnover",
      TRUE ~ "Insufficient cell-type data"
    ),
    tissue_status = factor(
      tissue_status,
      levels = c(
        "Tissue conserved",
        "Tissue turnover",
        "Other tissue",
        "Insufficient tissue data"
      )
    ),
    celltype_detail = factor(
      celltype_detail,
      levels = c(
        "Cell-type conserved",
        "Cell-type turnover",
        "Insufficient cell-type data"
      )
    )
  ) %>%
  distinct(gene, .keep_all = TRUE)

panel_gene_df <- panel_gene_df %>%
  mutate(
    celltype_subtype = case_when(
      celltype_detail == "Cell-type conserved" ~ paste0(
        "Conserved: ",
        str_to_title(str_replace_all(globally_shared_celltypes_text, "; ", " + "))
      ),
      celltype_detail == "Cell-type turnover" &
        ambiguous_parsimonious_solution ~
        "Ambiguous parsimonious turnover",
      celltype_detail == "Cell-type turnover" &
        minimum_distinct_celltypes == 2 ~ paste0(
          "Turnover: ", str_to_title(parsimonious_celltypes_text)
        ),
      celltype_detail == "Cell-type turnover" &
        minimum_distinct_celltypes >= 3 ~ paste0(
          "Complex turnover: ", str_to_title(parsimonious_celltypes_text)
        ),
      TRUE ~ "Insufficient cell-type data"
    )
  )

# ---------------- two-ring donut ----------------

panel_celltype_levels <- factor(
  c(
    "Cell-type conserved",
    "Cell-type turnover",
    "Insufficient cell-type data"
  ),
  levels = c(
    "Cell-type conserved",
    "Cell-type turnover",
    "Insufficient cell-type data"
  )
)

panel_donut_grid <- crossing(
  tissue_status = factor(
    c(
      "Tissue conserved",
      "Tissue turnover",
      "Other tissue",
      "Insufficient tissue data"
    ),
    levels = c(
      "Tissue conserved",
      "Tissue turnover",
      "Other tissue",
      "Insufficient tissue data"
    )
  ),
  celltype_detail = panel_celltype_levels
)

panel_donut_counts <- panel_gene_df %>%
  count(tissue_status, celltype_detail, name = "n_genes", .drop = FALSE) %>%
  right_join(
    panel_donut_grid,
    by = c("tissue_status", "celltype_detail")
  ) %>%
  mutate(n_genes = replace_na(n_genes, 0L)) %>%
  arrange(tissue_status, celltype_detail)

panel_inner_ring <- panel_donut_counts %>%
  arrange(tissue_status, celltype_detail) %>%
  mutate(
    ymax = cumsum(n_genes),
    ymin = lag(ymax, default = 0),
    label_position = (ymin + ymax) / 2,
    legend_key = as.character(celltype_detail)
  )

panel_outer_ring <- panel_donut_counts %>%
  group_by(tissue_status) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  arrange(tissue_status) %>%
  mutate(
    ymax = cumsum(n_genes),
    ymin = lag(ymax, default = 0),
    label_position = (ymin + ymax) / 2,
    legend_key = as.character(tissue_status)
  )

# Fisher's exact test excludes both insufficient categories.
panel_counts_for_test <- panel_donut_counts %>%
  filter(
    tissue_status != "Other tissue",
    tissue_status != "Insufficient tissue data",
    celltype_detail != "Insufficient cell-type data"
  ) %>%
  droplevels()

panel_fisher_matrix <- xtabs(
  n_genes ~ tissue_status + celltype_detail,
  data = panel_counts_for_test
)

panel_fisher_result <- fisher.test(panel_fisher_matrix)

panel_fisher_summary <- tibble(
  test = "Fisher's exact test",
  comparison = paste(
    "Cell-type conserved/turnover proportions between",
    "tissue-conserved and tissue-turnover genes"
  ),
  n_genes_plotted = nrow(panel_gene_df),
  n_genes_tested = sum(panel_counts_for_test$n_genes),
  odds_ratio = unname(panel_fisher_result$estimate),
  p_value = panel_fisher_result$p.value,
  significance = case_when(
    panel_fisher_result$p.value < 0.001 ~ "***",
    panel_fisher_result$p.value < 0.01 ~ "**",
    panel_fisher_result$p.value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
)

# Separate exact binomial tests at the tissue and cell-type levels.
# H0 for each test: conserved and turnover genes occur in a 1:1 ratio.
n_tissue_conserved <- sum(
  panel_gene_df$tissue_status == "Tissue conserved"
)
n_tissue_turnover <- sum(
  panel_gene_df$tissue_status == "Tissue turnover"
)

tissue_binomial_result <- binom.test(
  x = n_tissue_conserved,
  n = n_tissue_conserved + n_tissue_turnover,
  p = 0.5,
  alternative = "two.sided"
)

n_celltype_conserved <- sum(
  panel_gene_df$celltype_detail == "Cell-type conserved"
)
n_celltype_turnover <- sum(
  panel_gene_df$celltype_detail == "Cell-type turnover"
)

celltype_binomial_result <- binom.test(
  x = n_celltype_conserved,
  n = n_celltype_conserved + n_celltype_turnover,
  p = 0.5,
  alternative = "two.sided"
)

get_significance_label <- function(p) {
  case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

panel_binomial_summary <- bind_rows(
  tibble(
    level = "Tissue level",
    test = "Exact binomial test",
    null_hypothesis = "Conserved:turnover = 1:1",
    n_conserved = n_tissue_conserved,
    n_turnover = n_tissue_turnover,
    n_tested = n_tissue_conserved + n_tissue_turnover,
    observed_proportion_conserved = unname(tissue_binomial_result$estimate),
    expected_proportion_conserved = 0.5,
    confidence_interval_low = tissue_binomial_result$conf.int[1],
    confidence_interval_high = tissue_binomial_result$conf.int[2],
    p_value = tissue_binomial_result$p.value,
    significance = get_significance_label(tissue_binomial_result$p.value)
  ),
  tibble(
    level = "Cell-type level",
    test = "Exact binomial test",
    null_hypothesis = "Conserved:turnover = 1:1",
    n_conserved = n_celltype_conserved,
    n_turnover = n_celltype_turnover,
    n_tested = n_celltype_conserved + n_celltype_turnover,
    observed_proportion_conserved = unname(celltype_binomial_result$estimate),
    expected_proportion_conserved = 0.5,
    confidence_interval_low = celltype_binomial_result$conf.int[1],
    confidence_interval_high = celltype_binomial_result$conf.int[2],
    p_value = celltype_binomial_result$p.value,
    significance = get_significance_label(celltype_binomial_result$p.value)
  )
)

format_test_p <- function(p) {
  if (p < 0.001) {
    "< 0.001"
  } else {
    paste0("= ", formatC(p, digits = 3, format = "f"))
  }
}

panel_p_label <- paste0(
  "Tissue binomial P ", format_test_p(tissue_binomial_result$p.value),
  "; Cell-type binomial P ", format_test_p(celltype_binomial_result$p.value),
  "\nTissue-cell-type association: Fisher's P ",
  format_test_p(panel_fisher_result$p.value),
  "; n = ", sum(panel_counts_for_test$n_genes)
)

print(panel_binomial_summary)
print(panel_fisher_summary)

panel_donut_colors <- c(
  "Tissue conserved" = "#4C78A8",
  "Tissue turnover" = "#E76F51",
  "Other tissue" = "#9C755F",
  "Insufficient tissue data" = "#D9D9D9",
  "Cell-type conserved" = "#72B7B2",
  "Cell-type turnover" = "#F2CF5B",
  "Insufficient cell-type data" = "#BDBDBD"
)

two_ring_donut_plot <- ggplot() +
  geom_rect(
    data = panel_inner_ring,
    aes(
      ymin = ymin, ymax = ymax,
      xmin = 1.02, xmax = 1.72,
      fill = legend_key
    ),
    color = "white", linewidth = 0.8
  ) +
  geom_rect(
    data = panel_outer_ring,
    aes(
      ymin = ymin, ymax = ymax,
      xmin = 1.80, xmax = 2.50,
      fill = legend_key
    ),
    color = "white", linewidth = 0.8
  ) +
  geom_text(
    data = panel_inner_ring %>% filter(n_genes > 0),
    aes(x = 1.37, y = label_position, label = n_genes),
    size = 3.5, color = "#222222"
  ) +
  geom_text(
    data = panel_outer_ring %>% filter(n_genes > 0),
    aes(x = 2.15, y = label_position, label = n_genes),
    size = 3.8, fontface = "bold", color = "#222222"
  ) +
  annotate(
    "text", x = 0, y = 0,
    label = paste0("Genes\n(n = ", nrow(panel_gene_df), ")"),
    size = 4.0, fontface = "bold", lineheight = 1.0
  ) +
  coord_polar(theta = "y", start = 0) +
  xlim(c(0, 2.56)) +
  scale_fill_manual(
    values = panel_donut_colors,
    breaks = c(
      "Tissue conserved", "Tissue turnover", "Other tissue",
      "Insufficient tissue data",
      "Cell-type conserved", "Cell-type turnover",
      "Insufficient cell-type data"
    ),
    name = NULL
  ) +
  labs(
    title = "Tissue and cell-type expression turnover",
    subtitle = panel_p_label
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE))

# ---------------- tissue-subtype horizontal bar plot ----------------

panel_subtype_order <- c(
  "Conserved both gonads",
  "Conserved testis",
  "Conserved ovary",
  "Testis to ovary",
  "Testis to both gonads",
  "Ovary to testis",
  "Ovary to both gonads",
  "Both gonads to testis",
  "Both gonads to ovary",
  "Non-biased to testis",
  "Non-biased to ovary",
  "Non-biased to both gonads",
  "Testis to non-biased",
  "Ovary to non-biased",
  "Both gonads to non-biased",
  "Non-biased to multiple biased states",
  "Unresolved turnover",
  "Other tissue",
  "Insufficient tissue data"
)

stopifnot(all(panel_gene_df$tissue_subtype %in% panel_subtype_order))

panel_subtype_counts <- panel_gene_df %>%
  count(tissue_status, tissue_subtype, name = "n_genes") %>%
  mutate(
    tissue_subtype = factor(
      tissue_subtype,
      levels = rev(panel_subtype_order)
    )
  )

tissue_subtype_bar_plot <- ggplot(
  panel_subtype_counts,
  aes(x = n_genes, y = tissue_subtype, fill = tissue_status)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = n_genes),
    hjust = -0.20, size = 3.4, color = "#222222"
  ) +
  facet_grid(
    tissue_status ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_fill_manual(
    values = c(
      "Tissue conserved" = "#4C78A8",
      "Tissue turnover" = "#E76F51",
      "Other tissue" = "#9C755F",
      "Insufficient tissue data" = "#D9D9D9"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.18)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    title = "Detailed tissue-expression categories",
    x = "Number of genes",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
    axis.text.y = element_text(size = 8.5),
    panel.spacing.y = unit(6, "mm"),
    plot.margin = margin(12, 12, 12, 18)
  )

# ---------------- cell-type-subtype horizontal bar plot ----------------

celltype_subtype_counts <- panel_gene_df %>%
  filter(
    celltype_detail %in% c("Cell-type conserved", "Cell-type turnover")
  ) %>%
  group_by(celltype_detail, celltype_subtype) %>%
  summarise(
    n_genes = n_distinct(gene),
    genes = paste(sort(unique(gene)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    # Turnover groups are small enough to label directly with gene names.
    # Conserved groups retain compact count-only labels.
    bar_label = case_when(
      celltype_detail == "Cell-type turnover" ~ paste0(
        "n = ", n_genes, "\n", str_wrap(genes, width = 35)
      ),
      TRUE ~ paste0("n = ", n_genes)
    )
  ) %>%
  arrange(celltype_detail, n_genes, celltype_subtype) %>%
  mutate(
    # Prefix with the facet name so identical labels can be ordered
    # independently in the two facets.
    celltype_subtype_axis = paste(
      as.character(celltype_detail),
      celltype_subtype,
      sep = "|||"
    ),
    celltype_subtype_axis = factor(
      celltype_subtype_axis,
      levels = unique(celltype_subtype_axis)
    )
  )

celltype_subtype_bar_plot <- ggplot(
  celltype_subtype_counts,
  aes(x = n_genes, y = celltype_subtype_axis, fill = celltype_detail)
) +
  geom_col(width = 0.68, color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = bar_label),
    hjust = -0.05,
    size = 3.0,
    color = "#222222",
    lineheight = 0.95
  ) +
  facet_grid(
    celltype_detail ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(
    labels = function(x) sub("^.*\\|\\|\\|", "", x)
  ) +
  scale_fill_manual(
    values = c(
      "Cell-type conserved" = "#72B7B2",
      "Cell-type turnover" = "#F2CF5B"
    ),
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.75)),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Detailed cell-type categories",
    x = "Number of genes",
    y = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
    axis.text.y = element_text(size = 8.2),
    panel.spacing.y = unit(4, "mm"),
    plot.margin = margin(10, 110, 10, 18)
  )

# ---------------- combine and export ----------------

draw_donut_bar_figure <- function() {
  grid.newpage()
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        widths = unit(c(1.05, 1.00), "null"),
        heights = unit(c(1.20, 0.80), "null")
      )
    )
  )
  print(
    two_ring_donut_plot,
    newpage = FALSE,
    vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1)
  )
  print(
    tissue_subtype_bar_plot,
    newpage = FALSE,
    vp = viewport(layout.pos.row = 1, layout.pos.col = 2)
  )
  print(
    celltype_subtype_bar_plot,
    newpage = FALSE,
    vp = viewport(layout.pos.row = 2, layout.pos.col = 2)
  )
}

draw_donut_bar_figure()

pdf(
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.v14_nonbias_layout.pdf",
  width = 15.5, height = 10.0
)
draw_donut_bar_figure()
dev.off()



panel_donut_proportions <- panel_donut_counts %>%
  group_by(tissue_status) %>%
  mutate(
    tissue_total = sum(n_genes),
    proportion_within_tissue = if_else(
      tissue_total > 0,
      n_genes / tissue_total,
      NA_real_
    )
  ) %>%
  ungroup()

write.table(
  panel_donut_proportions,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.donut_counts.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  panel_subtype_counts,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.tissue_subtype_counts.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  celltype_subtype_counts,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.celltype_subtype_counts.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  panel_fisher_summary,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.fisher.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  panel_binomial_summary,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.binomial.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  panel_gene_df,
  "plots/DNA_pipe/sex_related_genes/non_bias_gain/separate_germ_subtypes/tissue_celltype_parsimonious.turnover_gene_labels.final.genes.v14_nonbias_layout.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
