# Updated: Insufficient-only somatic-expression rescue (100-gene v4 matched).
# Ancestry retains >=3 supporters, >50% support and a unique majority.
# One biased species qualifies for inclusion; one informative species alone
# is retained as insufficient_data, never called conserved.
# Existing species-level enrichment/low-expression thresholds are preserved in the first pass.
# Heatmap includes eligible genes; complete tables retain all orthogroups.
# Majority-based ancestral inference is not phylogenetic reconstruction.


# ================================================================
# Genome-wide gonadal expression conservation and turnover
# Six-species 1:1 orthologs
#
# Hab: average TPM across Hav, Ha3, Ha2, Ha1
# Her: average TPM across Hev, He3, He2, He1
#
# E_T = log2((Testis + 0.1)/(Male soma + 0.1))
# E_O = log2((Ovary  + 0.1)/(Female soma + 0.1))
#
# Enrichment cutoff = 1
# Low expression cutoff = 1 TPM
# Inclusion: at least one biased species; non-biased is informative
#
# No quantile normalization
# ================================================================


# ================================================================
# 1. Libraries and settings
# ================================================================

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(grid)

species_order <- c(
  "Hab", "Her", "Hhe", "Cin", "Dda", "Mbr"
)

gonad_enrichment_cutoff <- 1
low_expression_cutoff <- 1
pseudocount <- 0.1

min_biased_tissue_species <- 1
min_ancestral_support_species <- 3

outdir <- "plots/DNA_pipe/sex_related_genes/all_singlecopy_genes/all_ortholog_gonadal_bias/non_bias_gain"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ================================================================
# 2. Read six-species single-copy ortholog table
# ================================================================

ortho_file <- paste0(
  "mid_data/RNA_counts/tissue_specific_analysis/",
  "02.analysis_241013/01.orthoFinder_result/",
  "OrthoTable_SingleCopy.table"
)

if (!file.exists(ortho_file)) {
  stop("Ortholog file not found: ", ortho_file)
}

ortho_table <- read.table(
  ortho_file
)

if (ncol(ortho_table) != 7) {
  stop(
    "Expected 7 columns in ortholog table, found ",
    ncol(ortho_table)
  )
}

names(ortho_table) <- c(
  "Orthogroup",
  "Cin",
  "Dda",
  "Hhe",
  "Mbr",
  "Hab",
  "Her"
)

# Same gene ID conversions as the original code.

ortho_table <- ortho_table %>%
  mutate(
    across(everything(), as.character),
    
    Cin = str_replace_all(
      Cin, "Cint", "Cin"
    ),
    
    Dda = str_replace_all(
      Dda, "Ddac", "Dda"
    ),
    
    Hhe = str_remove(
      Hhe, "\\.t[0-9]+$"
    ),
    
    Hab = str_remove(
      Hab, "\\.t[0-9]+$"
    )
  )

if (anyDuplicated(ortho_table$Orthogroup)) {
  stop("Duplicated Orthogroup IDs.")
}

if (anyNA(ortho_table)) {
  stop("Missing gene IDs in ortholog table.")
}

cat(
  "\nTotal 1:1 orthologs:",
  nrow(ortho_table),
  "\n"
)

ortholog_long <- ortho_table %>%
  pivot_longer(
    cols = all_of(species_order),
    names_to = "species",
    values_to = "gene_id"
  ) %>%
  rename(gene = Orthogroup)


# ================================================================
# 3. Expression files
#
# Hab and Her each have four developmental stages.
# ================================================================

base_dir <- paste0(
  "mid_data/RNA_counts/tissue_specific_analysis/",
  "02.analysis_241013/01.tis_spec_genes/raw_tau"
)

tpm_files <- list(
  
  Cin = file.path(
    base_dir,
    "250325/papa.Fo.tau.tsv"
  ),
  
  Dda = file.path(
    base_dir,
    "250325/BJ.Fo.tau.tsv"
  ),
  
  Hhe = file.path(
    base_dir,
    "260529/HF.Fo.tau.tsv"
  ),
  
  Mbr = file.path(
    base_dir,
    "250325/HJ.Fo.tau.tsv"
  ),
  
  Hab = c(
    file.path(base_dir, "250325/Ha_Hav.Fo.tau.tsv"),
    file.path(base_dir, "Ha_Ha3.Fo.tau.tsv"),
    file.path(base_dir, "Ha_Ha2.Fo.tau.tsv"),
    file.path(base_dir, "Ha_Ha1.Fo.tau.tsv")
  ),
  
  Her = c(
    file.path(base_dir, "250325/He_Hev.Fo.tau.tsv"),
    file.path(base_dir, "He_He3.Fo.tau.tsv"),
    file.path(base_dir, "He_He2.Fo.tau.tsv"),
    file.path(base_dir, "He_He1.Fo.tau.tsv")
  )
  
)

all_input_files <- unlist(
  tpm_files,
  use.names = FALSE
)

missing_files <- all_input_files[
  !file.exists(all_input_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing expression files:\n",
    paste(missing_files, collapse = "\n")
  )
}


# ================================================================
# 4. Function: read one TPM file
#
# Species-specific gene ID conversions are preserved.
#
# Exclude:
# _Mp / _MP / _Ms / _Fa
#
# Retain:
# Mt / Fo / Mb / Mk / Fb / Fk
# ================================================================

read_species_expression <- function(
    path,
    species,
    stage
) {
  
  dat <- read.table(
    path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  
  dat <- dat %>%
    as.data.frame() %>%
    select(
      -contains("tau")
    ) %>%
    select(
      -matches("_Mp|_MP|_Ms|_Fa")
    ) %>%
    rownames_to_column("gene_id")
  
  # ------------------------------------------------------------
  # Gene ID conversion
  # ------------------------------------------------------------
  
  if (species == "Cin") {
    
    dat$gene_id <- dat$gene_id %>%
      str_remove("^gene-") %>%
      str_split_fixed("_", 2) %>%
      .[, 1] %>%
      paste0("_Cin")
    
  } else if (species == "Dda") {
    
    dat$gene_id <- dat$gene_id %>%
      str_remove("^gene-") %>%
      str_split_fixed("_", 2) %>%
      .[, 1] %>%
      paste0("_Dda")
    
  } else if (species == "Mbr") {
    
    dat$gene_id <- paste0(
      dat$gene_id, "_HJ"
    )
    
  } else if (species == "Her") {
    
    dat$gene_id <- paste0(
      dat$gene_id, "_He"
    )
    
  }
  
  # Hhe and Hab retain original gene IDs.
  
  # ------------------------------------------------------------
  # Convert sample columns into common tissue labels
  # ------------------------------------------------------------
  
  dat %>%
    pivot_longer(
      cols = -gene_id,
      names_to = "sample",
      values_to = "TPM"
    ) %>%
    mutate(
      species = species,
      stage = stage,
      TPM = as.numeric(TPM),
      
      tissue = case_when(
        str_detect(sample, "_Mt$") ~ "T",
        str_detect(sample, "_Fo$") ~ "O",
        str_detect(sample, "_Mb$") ~ "Mb",
        str_detect(sample, "_Mk$") ~ "Mk",
        str_detect(sample, "_Fb$") ~ "Fb",
        str_detect(sample, "_Fk$") ~ "Fk",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(tissue)
    ) %>%
    select(
      gene_id,
      species,
      stage,
      tissue,
      TPM
    )
  
}


# ================================================================
# 5. Read all developmental stages
# ================================================================

expression_raw <- purrr::imap_dfr(
  
  tpm_files,
  
  function(files, species) {
    
    stage_names <- if (species == "Hab") {
      
      c("Hav", "Ha3", "Ha2", "Ha1")
      
    } else if (species == "Her") {
      
      c("Hev", "He3", "He2", "He1")
      
    } else {
      
      species
      
    }
    
    purrr::map2_dfr(
      
      files,
      stage_names,
      
      function(path, stage) {
        
        read_species_expression(
          path = path,
          species = species,
          stage = stage
        )
        
      }
      
    )
    
  }
  
)


# ================================================================
# 6. Average TPM across developmental stages
#
# Hab:
# mean(Hav, Ha3, Ha2, Ha1)
#
# Her:
# mean(Hev, He3, He2, He1)
#
# Other species:
# one dataset, unchanged.
#
# Missing values:
# mean of available stages, NA if all stages missing.
# ================================================================

# Ensure one expression value per gene/stage/tissue.

# ================================================================
#  Remove unm records and average TPM
# ================================================================

# Remove unm records only if they are not in the
# single-copy ortholog table.

stopifnot(
  !("unm_Cin" %in% ortho_table$Cin),
  !("unm_Dda" %in% ortho_table$Dda)
)

expression_raw <- expression_raw %>%
  filter(
    !(species == "Cin" & gene_id == "unm_Cin"),
    !(species == "Dda" & gene_id == "unm_Dda")
  )

# Check duplicates again.

duplicate_check <- expression_raw %>%
  count(
    species,
    stage,
    gene_id,
    tissue
  ) %>%
  filter(n > 1)

print(duplicate_check)

stopifnot(nrow(duplicate_check) == 0)


# ================================================================
# Average TPM across developmental stages
# ================================================================

expression_long <- expression_raw %>%
  group_by(
    gene_id,
    species,
    tissue
  ) %>%
  summarise(
    TPM = if (all(is.na(TPM))) {
      NA_real_
    } else {
      mean(TPM, na.rm = TRUE)
    },
    
    n_stages = sum(!is.na(TPM)),
    
    .groups = "drop"
  )
# ================================================================
# 7. Check developmental-stage coverage
# ================================================================

stage_summary <- expression_long %>%
  
  filter(
    species %in% c("Hab", "Her")
  ) %>%
  
  group_by(
    species,
    tissue
  ) %>%
  
  summarise(
    
    n_genes = n(),
    
    n_four_stages = sum(
      n_stages == 4
    ),
    
    n_partial_stages = sum(
      n_stages > 0 &
        n_stages < 4
    ),
    
    n_no_expression_data = sum(
      n_stages == 0
    ),
    
    .groups = "drop"
    
  )

cat("\nHab/Her developmental-stage coverage:\n")
print(stage_summary, n = Inf)

# Check mean TPM of selected genes, if needed.

print(
  expression_long %>%
    filter(
      species %in% c("Hab", "Her"),
      tissue %in% c("T", "O")
    ) %>%
    head(12)
)


# ================================================================
# 8. Map expression to orthogroups
# ================================================================

ortholog_expression <- ortholog_long %>%
  left_join(
    expression_long %>%
      select(
        gene_id,
        species,
        tissue,
        TPM
      ),
    by = c("species", "gene_id")
  )

mapping_check <- ortholog_expression %>%
  
  group_by(species) %>%
  
  summarise(
    
    n_orthologs = n_distinct(gene),
    
    n_mapped = n_distinct(
      gene[!is.na(TPM)]
    ),
    
    .groups = "drop"
    
  ) %>%
  
  mutate(
    mapping_rate = n_mapped / n_orthologs
  )

cat("\nOrtholog expression mapping:\n")
print(mapping_check)

if (any(mapping_check$n_mapped == 0)) {
  stop(
    "At least one species has zero mapped orthologs."
  )
}


# ================================================================
# 9. Construct tissue expression table
# ================================================================

df_tissue <- ortholog_expression %>%
  select(
    gene,
    species,
    tissue,
    TPM
  ) %>%
  pivot_wider(
    names_from = tissue,
    values_from = TPM
  )

required_tissues <- c(
  "T", "O", "Mb", "Mk", "Fb", "Fk"
)

for (nm in required_tissues) {
  
  if (!nm %in% names(df_tissue)) {
    df_tissue[[nm]] <- NA_real_
  }
  
}

safe_max <- function(...) {
  
  x <- c(...)
  
  if (all(is.na(x))) {
    
    NA_real_
    
  } else {
    
    max(x, na.rm = TRUE)
    
  }
  
}

df_summary <- df_tissue %>%
  
  rowwise() %>%
  
  mutate(
    S_M = safe_max(Mb, Mk),
    S_F = safe_max(Fb, Fk)
  ) %>%
  
  ungroup() %>%
  
  select(
    gene,
    species,
    S_F,
    S_M,
    O,
    T
  )

stopifnot(
  nrow(df_summary) ==
    nrow(ortho_table) * length(species_order)
)


# ================================================================
# 10. Calculate gonadal enrichment
#
# E_T = log2((T + 0.1)/(S_M + 0.1))
# E_O = log2((O + 0.1)/(S_F + 0.1))
# ================================================================

df_summary <- df_summary %>%
  
  mutate(
    
    E_T = log2(
      (T + pseudocount) /
        (S_M + pseudocount)
    ),
    
    E_O = log2(
      (O + pseudocount) /
        (S_F + pseudocount)
    )
    
  )


# ================================================================
# 11. Species-level tissue classification
#
# Cutoff = 1
# ================================================================

df_summary <- df_summary %>%
  mutate(
    category = case_when(
      is.na(T) | is.na(O) |
        is.na(S_M) | is.na(S_F) ~
        "missing_data",
      
      pmax(T, O) <= low_expression_cutoff ~
        "low_expr",
      
      E_T > gonad_enrichment_cutoff &
        T > low_expression_cutoff &
        E_O > gonad_enrichment_cutoff &
        O > low_expression_cutoff ~
        "gonad_biased",
      
      E_T > gonad_enrichment_cutoff &
        T > low_expression_cutoff ~
        "testis_biased",
      
      E_O > gonad_enrichment_cutoff &
        O > low_expression_cutoff ~
        "ovary_biased",
      
      TRUE ~ "other"
    )
  )
cat("\nSpecies-level tissue categories:\n")

print(
  table(
    df_summary$species,
    df_summary$category
  )
)


# ================================================================
# 12. Informative tissue states
# ================================================================

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


# ================================================================
# 14. Complete classification for ALL orthogroups
# ================================================================

all_gene_classification <- tibble(
  gene = ortho_table$Orthogroup
) %>%
  
  left_join(
    df_gene_consensus,
    by = "gene"
  ) %>%
  
  mutate(
    
    n_species_used = replace_na(
      n_species_used,
      0L
    ),
    
    n_biased_species = replace_na(n_biased_species, 0L),
    n_non_bias = replace_na(n_non_bias, 0L),
    eligible_for_analysis = replace_na(eligible_for_analysis, FALSE),
    final_category = replace_na(
      final_category,
      "insufficient_data"
    )
    
  )

stopifnot(
  nrow(all_gene_classification) == nrow(ortho_table),
  !anyDuplicated(all_gene_classification$gene)
)


# ================================================================
# SECOND PASS: Insufficient-only rescue, matching the 100-gene v4 script
# The original species-level classification and first-pass gene
# classification above are kept unchanged.
# ================================================================
first_pass_categories <- all_gene_classification %>%
  transmute(gene, original_final_category = final_category)

rescue_gene_ids <- first_pass_categories %>%
  filter(original_final_category == "insufficient_data") %>%
  pull(gene)

# Only originally insufficient genes can be rescued. A low_expr
# species becomes non_bias when at least one sex-matched somatic
# reference exceeds 1 TPM. Require both references to be available,
# as in the 100-gene v4 analysis. Never rescue missing_data.
rescue_detail <- df_summary %>%
  filter(gene %in% rescue_gene_ids,
         category == "low_expr",
         !is.na(S_M), !is.na(S_F),
         pmax(S_M, S_F) > low_expression_cutoff) %>%
  transmute(gene, species, T, O, S_M, S_F,
            original_state = "low_expr", revised_state = "non_bias")

df_summary <- df_summary %>%
  mutate(category = if_else(
    gene %in% rescue_gene_ids & category == "low_expr" &
      !is.na(S_M) & !is.na(S_F) &
      pmax(S_M, S_F) > low_expression_cutoff,
    "other", category
  ))

# Re-run the same consensus rules on the rescued species states.
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


# ================================================================
# 14. Complete classification for ALL orthogroups
# ================================================================

all_gene_classification <- tibble(
  gene = ortho_table$Orthogroup
) %>%
  
  left_join(
    df_gene_consensus,
    by = "gene"
  ) %>%
  
  mutate(
    
    n_species_used = replace_na(
      n_species_used,
      0L
    ),
    
    n_biased_species = replace_na(n_biased_species, 0L),
    n_non_bias = replace_na(n_non_bias, 0L),
    eligible_for_analysis = replace_na(eligible_for_analysis, FALSE),
    final_category = replace_na(
      final_category,
      "insufficient_data"
    )
    
  )

stopifnot(
  nrow(all_gene_classification) == nrow(ortho_table),
  !anyDuplicated(all_gene_classification$gene)
)



# Audit: genes not initially insufficient must never change.
rescue_comparison <- first_pass_categories %>%
  left_join(
    all_gene_classification %>%
      transmute(gene, revised_final_category = final_category),
    by = "gene"
  ) %>%
  mutate(category_changed = original_final_category != revised_final_category)

unexpected_changes <- rescue_comparison %>%
  filter(original_final_category != "insufficient_data", category_changed)
if (nrow(unexpected_changes) > 0) {
  stop("Unexpected classification changes outside originally insufficient genes")
}

write_tsv(rescue_detail,
          file.path(outdir, "tissue_insufficient_rescued_species.tsv"))
write_tsv(rescue_comparison,
          file.path(outdir, "tissue_insufficient_before_after.tsv"))
cat("\nInsufficient-only rescue transition counts:\n")
print(rescue_comparison %>%
        count(original_final_category, revised_final_category, name = "n_genes"),
      n = Inf)

# Dedicated tables for the new inclusion rule and non-biased ancestry.
genes_with_any_bias <- all_gene_classification %>% filter(eligible_for_analysis)
non_bias_gain_genes <- genes_with_any_bias %>%
  filter(inferred_ancestral_state == "non_bias")
non_bias_gain_details <- df_state %>%
  semi_join(non_bias_gain_genes, by = "gene") %>%
  left_join(non_bias_gain_genes %>% select(gene, final_category), by = "gene")
write_tsv(genes_with_any_bias, file.path(outdir, "genes_with_any_bias.tsv"))
write_tsv(non_bias_gain_genes, file.path(outdir, "non_bias_ancestral_gain_genes.tsv"))
write_tsv(non_bias_gain_details, file.path(outdir, "non_bias_ancestral_gain_details.tsv"))

# ================================================================
# 15. Category statistics
# ================================================================

category_levels <- c(
  "conserved_gonad_biased",
  "conserved_testis_biased",
  "conserved_ovary_biased",
  "testis_to_ovary_turnover",
  "testis_to_dual_turnover",
  "ovary_to_testis_turnover",
  "ovary_to_dual_turnover",
  "dual_to_testis_turnover",
  "dual_to_ovary_turnover",
  "non_bias_to_testis_turnover",
  "non_bias_to_ovary_turnover",
  "non_bias_to_dual_turnover",
  "testis_to_non_bias_turnover",
  "ovary_to_non_bias_turnover",
  "dual_to_non_bias_turnover",
  "non_bias_to_multiple_turnover",
  "unresolved_turnover",
  "no_bias_observed",
  "insufficient_data"
)

category_counts <- all_gene_classification %>%
  
  mutate(
    final_category = factor(
      final_category,
      levels = category_levels
    )
  ) %>%
  
  count(
    final_category,
    name = "n_genes",
    .drop = FALSE
  ) %>%
  
  mutate(
    percentage = 100 * n_genes / sum(n_genes)
  )

cat("\n========================================\n")
cat("Genome-wide category counts\n")
cat("========================================\n")

print(
  category_counts,
  n = Inf
)


# ================================================================
# 16. Select genes for heatmap
#
# At least one biased species AND sufficient tissue data
# Insufficient genes remain in the full classification table.
# ================================================================

genes_for_heatmap <- all_gene_classification %>%
  filter(
    eligible_for_analysis,
    final_category != "insufficient_data"
  ) %>%
  pull(gene)

cat(
  "\nGenes included in heatmap:",
  length(genes_for_heatmap),
  "\n"
)

if (length(genes_for_heatmap) == 0) {
  stop(
    "No genes meet the heatmap inclusion criterion."
  )
}


# ================================================================
# 17. Heatmap category labels
# ================================================================

category_label_map <- c(
  conserved_gonad_biased = "Conserved both gonads",
  conserved_testis_biased = "Conserved testis",
  conserved_ovary_biased = "Conserved ovary",
  testis_to_ovary_turnover = "Testis to ovary",
  testis_to_dual_turnover = "Testis to both gonads",
  ovary_to_testis_turnover = "Ovary to testis",
  ovary_to_dual_turnover = "Ovary to both gonads",
  dual_to_testis_turnover = "Both gonads to testis",
  dual_to_ovary_turnover = "Both gonads to ovary",
  non_bias_to_testis_turnover = "Non-biased to testis",
  non_bias_to_ovary_turnover = "Non-biased to ovary",
  non_bias_to_dual_turnover = "Non-biased to both gonads",
  testis_to_non_bias_turnover = "Testis to non-biased",
  ovary_to_non_bias_turnover = "Ovary to non-biased",
  dual_to_non_bias_turnover = "Both gonads to non-biased",
  non_bias_to_multiple_turnover = "Non-biased to multiple biased states",
  unresolved_turnover = "Unresolved turnover",
  no_bias_observed = "Other tissue",
  insufficient_data = "Insufficient tissue data"
)

category_info <- all_gene_classification %>%
  
  filter(
    gene %in% genes_for_heatmap
  ) %>%
  
  mutate(
    plot_category = unname(
      category_label_map[final_category]
    )
  )

plot_category_levels <- unname(
  category_label_map
)

category_info <- category_info %>%
  
  mutate(
    plot_category = factor(
      plot_category,
      levels = plot_category_levels
    )
  )

category_counts_heatmap <- category_info %>%
  
  count(
    plot_category,
    name = "n_genes",
    .drop = FALSE
  )

print(
  category_counts_heatmap,
  n = Inf
)


# ================================================================
# 18. Build heatmap matrix
#
# O = ovary enrichment
# T = testis enrichment
# ================================================================

heatmap_long <- df_summary %>%
  
  filter(
    gene %in% genes_for_heatmap
  ) %>%
  
  select(
    gene,
    species,
    E_O,
    E_T
  ) %>%
  
  pivot_longer(
    cols = c(E_O, E_T),
    names_to = "metric",
    values_to = "value"
  ) %>%
  
  mutate(
    
    metric = recode(
      metric,
      E_O = "O",
      E_T = "T"
    ),
    
    colname = paste(
      species,
      metric,
      sep = "_"
    )
    
  ) %>%
  
  select(
    gene,
    colname,
    value
  )

heatmap_mat <- heatmap_long %>%
  
  pivot_wider(
    names_from = colname,
    values_from = value
  ) %>%
  
  column_to_rownames("gene") %>%
  
  as.matrix()

expected_cols <- as.vector(
  rbind(
    paste0(species_order, "_O"),
    paste0(species_order, "_T")
  )
)

for (nm in setdiff(expected_cols, colnames(heatmap_mat))) {
  
  missing_col <- matrix(
    NA_real_,
    nrow = nrow(heatmap_mat),
    ncol = 1,
    dimnames = list(
      rownames(heatmap_mat),
      nm
    )
  )
  
  heatmap_mat <- cbind(
    heatmap_mat,
    missing_col
  )
  
}

heatmap_mat <- heatmap_mat[
  ,
  expected_cols,
  drop = FALSE
]


# ================================================================
# 19. Match categories to matrix rows
# ================================================================

category_info <- category_info[
  match(
    rownames(heatmap_mat),
    category_info$gene
  ),
]

stopifnot(
  identical(
    rownames(heatmap_mat),
    category_info$gene
  )
)

category_label_df <- category_counts_heatmap %>%
  
  mutate(
    category_label = paste0(
      as.character(plot_category),
      " (n = ",
      n_genes,
      ")"
    )
  )

category_label_lookup <- setNames(
  category_label_df$category_label,
  as.character(category_label_df$plot_category)
)

category_label_levels <- unname(
  category_label_lookup[plot_category_levels]
)

category_label_levels <- category_label_levels[
  !is.na(category_label_levels)
]

category_info <- category_info %>%
  
  mutate(
    plot_category_label = factor(
      unname(
        category_label_lookup[
          as.character(plot_category)
        ]
      ),
      levels = category_label_levels
    )
  )

row_split <- category_info$plot_category_label


# ================================================================
# 20. Heatmap color
#
# Blue: low expression enrichment
# Near-white: cutoff = 1
# Red: high expression enrichment
# ================================================================

heatmap_mat_plot <- pmax(
  pmin(heatmap_mat, 7),
  -7
)

col_fun <- colorRamp2(
  
  c(
    -7,
    gonad_enrichment_cutoff,
    7
  ),
  
  c(
    "#0067a5",
    "#f3f3f4",
    "#be0032"
  )
  
)


# ================================================================
# 21. Species annotation
# ================================================================

column_labels <- rep(
  c("O", "T"),
  length(species_order)
)

species_top <- rep(
  species_order,
  each = 2
)

top_annotation <- HeatmapAnnotation(
  
  Species = anno_text(
    species_top,
    
    gp = gpar(
      fontsize = 10,
      fontface = "bold"
    )
  ),
  
  annotation_height = unit(6, "mm"),
  
  show_annotation_name = FALSE
  
)


# ================================================================
# 22. Construct heatmap
#
# Cluster genes within each category.
# No dendrogram displayed.
# ================================================================

# Plotting masks; original expression values remain available for clustering.
tissue_low_mask <- matrix(FALSE, nrow(heatmap_mat), ncol(heatmap_mat),
                          dimnames = dimnames(heatmap_mat))
tissue_missing_mask <- is.na(heatmap_mat)
mask_rows <- df_summary %>% filter(gene %in% rownames(heatmap_mat))
for (k in seq_len(nrow(mask_rows))) {
  g <- as.character(mask_rows$gene[k])
  sp <- as.character(mask_rows$species[k])
  cols <- intersect(paste0(sp, c("_O", "_T")), colnames(heatmap_mat))
  if (mask_rows$category[k] == "low_expr") tissue_low_mask[g, cols] <- TRUE
  if (mask_rows$category[k] == "missing_data") tissue_missing_mask[g, cols] <- TRUE
}

ht <- Heatmap(
  
  heatmap_mat_plot,
  
  name = "log2(gonad/soma)",
  
  col = col_fun,
  
  na_col = "white",
  layer_fun = function(j, i, x, y, width, height, fill) {
    low <- tissue_low_mask[cbind(i, j)]
    missing <- tissue_missing_mask[cbind(i, j)]
    if (any(low)) grid.rect(x[low], y[low], width, height,
                            gp = gpar(fill = "#E5E5E5", col = NA))
    if (any(missing)) grid.rect(x[missing], y[missing], width, height,
                                gp = gpar(fill = "white", col = NA))
  },
  
  row_split = row_split,
  
  cluster_rows = F,
  
  cluster_row_slices = FALSE,
  
  show_row_dend = FALSE,
  
  row_title_side = "left",
  row_title_rot = 0,
  
  row_title_gp = gpar(
    fontsize = 8,
    fontface = "bold"
  ),
  
  row_gap = unit(2, "mm"),
  
  show_row_names = FALSE,
  
  cluster_columns = FALSE,
  
  show_column_dend = FALSE,
  
  show_column_names = TRUE,
  
  column_labels = column_labels,
  
  column_names_gp = gpar(
    fontsize = 9
  ),
  
  top_annotation = top_annotation,
  
  rect_gp = gpar(
    col = NA
  ),
  
  column_title =
    "Genome-wide gonadal expression of 1:1 orthologs",
  
  column_title_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  
  use_raster = TRUE,
  
  raster_quality = 3,
  
  heatmap_legend_param = list(
    
    title = "log2(gonad/soma)",
    
    at = c(-7, 0, 1, 4, 7),
    
    legend_height = unit(
      35,
      "mm"
    )
    
  )
  
)

draw(
  ht,
  heatmap_legend_side = "right"
)


# ================================================================
# 23. Save heatmap PDF
# ================================================================

pdf(
  file.path(
    outdir,
    "all_1to1_orthologs.gonadal_bias.cutoff1.v2.pdf"
  ),
  
  width = 13,
  height = 16
)

draw(
  ht,
  heatmap_legend_side = "right"
)

dev.off()


# ================================================================
# 24. Export results
# ================================================================

write.table(
  
  df_summary,
  
  file.path(
    outdir,
    "all_orthologs.expression_summary.cutoff1.v2.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

write.table(
  
  all_gene_classification,
  
  file.path(
    outdir,
    "all_orthologs.evolutionary_classification.cutoff1.v2.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

write.table(
  
  category_counts,
  
  file.path(
    outdir,
    "all_orthologs.category_counts.cutoff1.v2.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

write.table(
  
  category_info,
  
  file.path(
    outdir,
    "all_orthologs.heatmap_categories.cutoff1.v2.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

write.table(
  
  ortho_table,
  
  file.path(
    outdir,
    "all_orthologs.gene_ID_mapping.v2.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

write.table(
  
  stage_summary,
  
  file.path(
    outdir,
    "Hab_Her.developmental_stage_coverage.v2.tsv"
  ),
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)


# ================================================================
# 25. Export gene lists by evolutionary category
# ================================================================

for (cat_name in category_levels) {
  
  gene_list <- all_gene_classification %>%
    filter(
      final_category == cat_name
    ) %>%
    pull(gene)
  
  writeLines(
    gene_list,
    file.path(
      outdir,
      paste0(cat_name, ".genes.txt")
    )
  )
  
}


# ================================================================
# 26. Final summary
# ================================================================

cat("\n========================================\n")
cat("FINAL GENOME-WIDE ANALYSIS SUMMARY\n")
cat("========================================\n")

cat(
  "Total 1:1 orthologs:",
  nrow(ortho_table),
  "\n"
)

cat(
  "Genes with >=1 biased species:",
  length(genes_for_heatmap),
  "\n"
)

cat(
  "Gonadal enrichment cutoff:",
  gonad_enrichment_cutoff,
  "\n"
)

cat(
  "Hab/Her expression:",
  "mean TPM across four developmental stages\n"
)

print(
  category_counts,
  n = Inf
)

cat(
  "\nOutput directory:",
  outdir,
  "\n"
)

cat("========================================\n")




# ================================================================
# 27. Export gene Symbols by cross-species category for GO
#
# Mapping:
# Orthogroup -> Mbr gene ID -> Symbol
#
# Required objects:
#   all_gene_classification
#   ortho_table
#   category_levels
#   eligible_for_analysis
# ================================================================

library(tidyverse)

go_outdir <- file.path(outdir, "GO")

dir.create(
  go_outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------
# 1. Read Mbr gene Symbol annotation
# ------------------------------------------------

loc_data <- read.table(
  "ref/eggnog/03.HJ/HJ_withEggnogPannzer2SymbolInfo.tsv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  na.strings = c("", "NA", "<NA>", "-")
)

stopifnot(
  all(c("geneID", "Symbol") %in% colnames(loc_data))
)

symbol_map <- loc_data %>%
  transmute(
    Mbr_geneID = str_trim(as.character(geneID)),
    Symbol = str_trim(as.character(Symbol))
  ) %>%
  mutate(
    Symbol = if_else(
      is.na(Symbol) | Symbol %in% c("", "-", "NA", "<NA>"),
      NA_character_,
      Symbol
    )
  ) %>%
  filter(
    !is.na(Mbr_geneID),
    nzchar(Mbr_geneID)
  ) %>%
  distinct()

# Check whether one gene has multiple different Symbols.
symbol_conflicts <- symbol_map %>%
  filter(!is.na(Symbol)) %>%
  group_by(Mbr_geneID) %>%
  summarise(
    n_symbols = n_distinct(Symbol),
    Symbols = paste(sort(unique(Symbol)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_symbols > 1)

if (nrow(symbol_conflicts) > 0) {
  
  write_tsv(
    symbol_conflicts,
    file.path(go_outdir, "ERROR_conflicting_Symbols.tsv")
  )
  
  stop(
    "Some Mbr genes have multiple Symbols. ",
    "Check ERROR_conflicting_Symbols.tsv."
  )
}

# If both missing and non-missing records exist, keep the annotation.
symbol_map <- symbol_map %>%
  arrange(Mbr_geneID, is.na(Symbol)) %>%
  distinct(Mbr_geneID, .keep_all = TRUE)

# ------------------------------------------------
# 2. Join cross-species categories to Mbr gene IDs
# ------------------------------------------------

stopifnot(
  !anyDuplicated(all_gene_classification$gene),
  !anyDuplicated(ortho_table$Orthogroup),
  all(
    all_gene_classification$gene %in%
      ortho_table$Orthogroup
  ),
  all(
    all_gene_classification$final_category %in%
      category_levels
  )
)

go_mapping <- all_gene_classification %>%
  select(
    gene,
    final_category,
    n_species_used,
    n_biased_species,
    eligible_for_analysis
  ) %>%
  left_join(
    ortho_table %>%
      transmute(
        gene = Orthogroup,
        Mbr_geneID = str_remove(
          str_trim(as.character(Mbr)),
          "_HJ$"
        )
      ),
    by = "gene"
  ) %>%
  left_join(
    symbol_map,
    by = "Mbr_geneID"
  )

stopifnot(
  nrow(go_mapping) == nrow(all_gene_classification),
  !anyNA(go_mapping$Mbr_geneID),
  all(nzchar(go_mapping$Mbr_geneID))
)

if (all(is.na(go_mapping$Symbol))) {
  stop("No Symbols mapped. Check Mbr gene IDs against loc_data$geneID.")
}

# Full mapping table.
write_tsv(
  go_mapping,
  file.path(
    go_outdir,
    "all_categories_Orthogroup_Mbr_Symbol_mapping.tsv"
  )
)

# Genes without a Symbol.
write_tsv(
  go_mapping %>%
    filter(is.na(Symbol)),
  file.path(go_outdir, "genes_without_Symbol.tsv")
)

# ------------------------------------------------
# 3. Helper functions
# ------------------------------------------------

get_symbols <- function(x) {
  sort(unique(x[!is.na(x) & nzchar(x)]))
}

write_gene_list <- function(ids, path) {
  writeLines(
    get_symbols(ids),
    con = path,
    useBytes = TRUE
  )
}

# ------------------------------------------------
# 4. Export one Symbol list per category
# ------------------------------------------------

for (cat_name in category_levels) {
  
  cat_data <- go_mapping %>%
    filter(final_category == cat_name)
  
  # Plain list for copying into a GO website:
  # one Symbol per line, no header, no duplicates.
  write_gene_list(
    cat_data$Symbol,
    file.path(
      go_outdir,
      paste0(cat_name, ".Symbols.txt")
    )
  )
  
  # Detailed mapping for checking.
  write_tsv(
    cat_data,
    file.path(
      go_outdir,
      paste0(cat_name, ".mapping.tsv")
    )
  )
}

# ------------------------------------------------
# 5. Export background lists
# ------------------------------------------------

# Background A:
# All genes eligible for conservation/turnover classification.
write_gene_list(
  go_mapping %>%
    filter(
      eligible_for_analysis,
      final_category != "insufficient_data"
    ) %>%
    pull(Symbol),
  file.path(
    go_outdir,
    "BACKGROUND_informative_genes.Symbols.txt"
  )
)

# Background B:
# All single-copy Orthogroups in the original analysis.
write_gene_list(
  go_mapping$Symbol,
  file.path(
    go_outdir,
    "BACKGROUND_all_orthologs.Symbols.txt"
  )
)

# ------------------------------------------------
# 6. Export mapping statistics
# ------------------------------------------------

go_counts <- go_mapping %>%
  mutate(
    final_category = factor(
      final_category,
      levels = category_levels
    )
  ) %>%
  group_by(
    final_category,
    .drop = FALSE
  ) %>%
  summarise(
    n_orthogroups = n(),
    n_with_Symbol = sum(!is.na(Symbol)),
    n_without_Symbol = sum(is.na(Symbol)),
    n_unique_Symbols = n_distinct(
      Symbol,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_tsv(
  go_counts,
  file.path(go_outdir, "GO_gene_list_counts.tsv")
)

print(go_counts, n = Inf)

cat(
  "\nGO gene lists saved to:\n",
  normalizePath(
    go_outdir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)


# ================================================================
# 按推断的祖先表达状态导出 GO Symbol 列表
# 依赖第 27 部分生成的 go_mapping
# ================================================================

ancestral_go_dir <- file.path(outdir, "GO", "ancestral_states")

dir.create(
  ancestral_go_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 1. 合并祖先状态与 Symbol
ancestral_go_mapping <- all_gene_classification %>%
  filter(
    eligible_for_analysis,
    !is.na(inferred_ancestral_state),
    final_category != "insufficient_data"
  ) %>%
  select(
    gene,
    inferred_ancestral_state,
    final_category
  ) %>%
  left_join(
    go_mapping %>%
      select(gene, Mbr_geneID, Symbol) %>%
      distinct(),
    by = "gene"
  ) %>%
  mutate(
    Symbol = str_trim(as.character(Symbol)),
    Symbol = if_else(
      is.na(Symbol) | Symbol %in% c("", "-", "NA", "<NA>"),
      NA_character_,
      Symbol
    )
  )

# 2. 内部状态名称与输出文件名
ancestral_file_names <- c(
  dual     = "ancestral_both_gonads",
  testis   = "ancestral_testis_biased",
  ovary    = "ancestral_ovary_biased",
  non_bias = "ancestral_non_biased"
)

# 3. 分别导出 Symbol 和详细映射表
for (state in names(ancestral_file_names)) {
  
  state_data <- ancestral_go_mapping %>%
    filter(inferred_ancestral_state == state)
  
  state_symbols <- state_data %>%
    filter(!is.na(Symbol)) %>%
    pull(Symbol) %>%
    unique() %>%
    sort()
  
  file_prefix <- unname(ancestral_file_names[state])
  
  writeLines(
    state_symbols,
    file.path(
      ancestral_go_dir,
      paste0(file_prefix, ".Symbols.txt")
    ),
    useBytes = TRUE
  )
  
  write_tsv(
    state_data,
    file.path(
      ancestral_go_dir,
      paste0(file_prefix, ".mapping.tsv")
    )
  )
}

# 4. 汇总各组数量，包括空组
ancestral_go_counts <- tibble(
  inferred_ancestral_state = names(ancestral_file_names),
  output_group = unname(ancestral_file_names)
) %>%
  left_join(
    ancestral_go_mapping %>%
      group_by(inferred_ancestral_state) %>%
      summarise(
        n_orthogroups = n_distinct(gene),
        n_orthogroups_with_symbol =
          n_distinct(gene[!is.na(Symbol)]),
        n_unique_symbols =
          n_distinct(Symbol, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "inferred_ancestral_state"
  ) %>%
  mutate(
    across(starts_with("n_"), ~ replace_na(.x, 0L))
  )

write_tsv(
  ancestral_go_counts,
  file.path(ancestral_go_dir, "ancestral_GO_counts.tsv")
)

print(ancestral_go_counts, n = Inf, width = Inf)

cat(
  "\n输出目录：\n",
  normalizePath(ancestral_go_dir, winslash = "/", mustWork = TRUE),
  "\n"
)


# ------------------------------------------------
# 7. Print one category directly for copying
# ------------------------------------------------

# Change this to any category in category_levels.
category_to_print <- "conserved_ovary_biased"

symbols_to_print <- go_mapping %>%
  filter(final_category == category_to_print) %>%
  pull(Symbol) %>%
  get_symbols()

cat(
  "\n",
  paste(symbols_to_print, collapse = "\n"),
  "\n",
  sep = ""
)


# 只保留符合正式分析标准的基因
analysis_genes <- all_gene_classification %>%
  filter(
    eligible_for_analysis,
    final_category != "insufficient_data"
  )

# 重新统计各类别及百分比
analysis_category_counts <- analysis_genes %>%
  count(final_category, name = "n_genes") %>%
  mutate(
    percentage = 100 * n_genes / sum(n_genes)
  ) %>%
  arrange(desc(n_genes))

print(analysis_category_counts, n = Inf)

write_tsv(
  analysis_category_counts,
  file.path(outdir, "analysis_only_category_counts.tsv")
)
