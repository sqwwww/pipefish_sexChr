library(tidyverse)
library(patchwork)


# ============================================================
# 1. Read single-copy ortholog table
# ============================================================

ortho_table <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.orthoFinder_result/OrthoTable_SingleCopy.table"
)

names(ortho_table) <- c(
  "Orthogroup",
  "Cin",
  "Dda",
  "Hhe",
  "Mbr",
  "Hab",
  "Her"
)

ortho_table <- lapply(
  ortho_table,
  function(x) gsub("Cint", "Cin", x)
) %>%
  as.data.frame()

ortho_table <- lapply(
  ortho_table,
  function(x) gsub("Ddac", "Dda", x)
) %>%
  as.data.frame()

# Remove transcript suffix
ortho_table$Hhe <- gsub(
  "\\.t[0-9]",
  "",
  ortho_table$Hhe
)

ortho_table$Hab <- gsub(
  "\\.t[0-9]",
  "",
  ortho_table$Hab
)

ortho_table2 <- ortho_table


# ============================================================
# 2. Add bulk RNA-seq expression
#
# IMPORTANT:
# Remove _Mp / _MP / _Ms / _Fa
# before calculating testis relative expression.
# ============================================================


# ------------------------------------------------------------
# Cin
# ------------------------------------------------------------

Cin_df <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/papa.Fo.tau.tsv"
)

new_df <- Cin_df %>%
  select(-contains("tau"))

# Remove pouch-related tissues
new_df <- new_df %>%
  select(
    -matches("_Mp|_MP|_Ms|_Fa")
  )

new_df$geneID <- rownames(new_df) %>%
  str_remove("gene-") %>%
  str_split("_") %>%
  map_chr(1) %>%
  paste0("_Cin")

table(ortho_table$Cin %in% new_df$geneID)

ortho_table2 <- ortho_table2 %>%
  left_join(
    new_df,
    by = c("Cin" = "geneID")
  )


# ------------------------------------------------------------
# Dda
# ------------------------------------------------------------

Dda_df <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/BJ.Fo.tau.tsv"
)

new_df <- Dda_df %>%
  select(-contains("tau"))

# Remove pouch-related tissues
new_df <- new_df %>%
  select(
    -matches("_Mp|_MP|_Ms|_Fa")
  )

new_df$geneID <- rownames(new_df) %>%
  str_remove("gene-") %>%
  str_split("_") %>%
  map_chr(1) %>%
  paste0("_Dda")

table(ortho_table$Dda %in% new_df$geneID)

ortho_table2 <- ortho_table2 %>%
  left_join(
    new_df,
    by = c("Dda" = "geneID")
  )


# ------------------------------------------------------------
# Hhe
# ------------------------------------------------------------

Hhe_df <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/260529/HF.Fo.tau.tsv"
)

new_df <- Hhe_df %>%
  select(-contains("tau"))

# Remove pouch-related tissues
new_df <- new_df %>%
  select(
    -matches("_Mp|_MP|_Ms|_Fa")
  )

new_df$geneID <- rownames(new_df)

table(ortho_table$Hhe %in% new_df$geneID)

ortho_table2 <- ortho_table2 %>%
  left_join(
    new_df,
    by = c("Hhe" = "geneID")
  )


# ------------------------------------------------------------
# Mbr
# ------------------------------------------------------------

Mbr_df <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/HJ.Fo.tau.tsv"
)

new_df <- Mbr_df %>%
  select(-contains("tau"))

# Remove pouch-related tissues
new_df <- new_df %>%
  select(
    -matches("_Mp|_MP|_Ms|_Fa")
  )

new_df$geneID <- rownames(new_df) %>%
  paste0("_HJ")

table(ortho_table$Mbr %in% new_df$geneID)

ortho_table2 <- ortho_table2 %>%
  left_join(
    new_df,
    by = c("Mbr" = "geneID")
  )


# ------------------------------------------------------------
# Hab
# ------------------------------------------------------------

Hab_df <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/Ha_Hav.Fo.tau.tsv"
)

new_df <- Hab_df %>%
  select(-contains("tau"))

# Remove pouch-related tissues
new_df <- new_df %>%
  select(
    -matches("_Mp|_MP|_Ms|_Fa")
  )

new_df$geneID <- rownames(new_df)

table(ortho_table$Hab %in% new_df$geneID)

ortho_table2 <- ortho_table2 %>%
  left_join(
    new_df,
    by = c("Hab" = "geneID")
  )


# ------------------------------------------------------------
# Her
# ------------------------------------------------------------

Her_df <- read.table(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/He_Hev.Fo.tau.tsv"
)

new_df <- Her_df %>%
  select(-contains("tau"))

# Remove pouch-related tissues
new_df <- new_df %>%
  select(
    -matches("_Mp|_MP|_Ms|_Fa")
  )

new_df$geneID <- rownames(new_df) %>%
  paste0("_He")

table(ortho_table$Her %in% new_df$geneID)

ortho_table2 <- ortho_table2 %>%
  left_join(
    new_df,
    by = c("Her" = "geneID")
  )


# ============================================================
# 3. Add Mbr chromosome / SDR / stratum information
# ============================================================

loc_data <- read.table(
  "ref/eggnog/03.HJ/HJ_withEggnogPannzer2SymbolInfo.tsv"
)

table(loc_data$type_new)


# ------------------------------------------------------------
# Broad classification
# ------------------------------------------------------------

loc_data$type1 <- case_when(
  loc_data$type_new == "PAR" ~ "PAR",
  loc_data$type_new %in% c("S0", "S1", "S2", "S3") ~ "SDR",
  TRUE ~ "auto"
)


# ------------------------------------------------------------
# Detailed strata
# ------------------------------------------------------------

loc_data$type2 <- case_when(
  loc_data$type_new == "PAR" ~ "PAR",
  loc_data$type_new == "S0" ~ "S0",
  loc_data$type_new == "S1" ~ "S1",
  loc_data$type_new == "S2" ~ "S2",
  loc_data$type_new == "S3" ~ "S3",
  TRUE ~ "auto"
)


table(loc_data$type1)
table(loc_data$type2)


loc_data$geneID <- paste0(
  loc_data$geneID,
  "_HJ"
)

table(
  loc_data$geneID %in% ortho_table2$Mbr
)


ortho_table3 <- ortho_table2 %>%
  left_join(
    loc_data,
    by = c("Mbr" = "geneID")
  )


# IMPORTANT:
# Do NOT use global na.omit()
ortho_table4 <- ortho_table3


# Mbr gene ID without _HJ suffix
ortho_table4$Mbr_id <- ortho_table4$Mbr %>%
  str_remove("_HJ")


# ============================================================
# 4. Add ancestral LG information using Dda
#
# LG5 / LG12 / LG13 / LG15
# ============================================================

Dda_loc <- read.table(
  "mid_data/newGene_pipeline/02.geneLocation/Ddac.gene_lcation_info2"
)

Dda_loc <- Dda_loc %>%
  select(
    V1,
    V4
  )

Dda_loc$V4 <- Dda_loc$V4 %>%
  str_replace(
    pattern = "Ddac",
    replacement = "Dda"
  )

names(Dda_loc) <- c(
  "DdaChrom",
  "DdaGene"
)


table(
  Dda_loc$DdaGene %in% ortho_table4$Dda
)


Dda_loc$Dda_type <- case_when(
  Dda_loc$DdaChrom == "5"  ~ "5",
  Dda_loc$DdaChrom == "12" ~ "12",
  Dda_loc$DdaChrom == "13" ~ "13",
  Dda_loc$DdaChrom == "15" ~ "15",
  TRUE ~ "auto"
)


table(Dda_loc$Dda_type)


ortho_table4 <- ortho_table4 %>%
  left_join(
    Dda_loc,
    by = c("Dda" = "DdaGene")
  )


ortho_table4$LG_type <- ortho_table4$Dda_type


# Same correction used in the original analysis:
# Dda Chr12 genes outside the Mbr sex chromosome
# are treated as autosomal background.
ortho_table4$LG_type[
  ortho_table4$LG_type == "12" &
    ortho_table4$type2 == "auto"
] <- "auto"


table(
  ortho_table4$LG_type,
  useNA = "ifany"
)


# ============================================================
# 5. Extract RAW expression
#
# NO quantile normalization.
#
# _Mp / _MP / _Ms / _Fa have already been removed.
# ============================================================

Hab_expr <- ortho_table4 %>%
  select(
    matches("^Hav_")
  )

Her_expr <- ortho_table4 %>%
  select(
    matches("^Hev_")
  )

Hhe_expr <- ortho_table4 %>%
  select(
    matches("^Hha_")
  )

Cin_expr <- ortho_table4 %>%
  select(
    matches("^Cha_")
  )

Dda_expr <- ortho_table4 %>%
  select(
    matches("^Dda_")
  ) %>%
  select(
    -any_of("Dda_type")
  )

Mbr_expr <- ortho_table4 %>%
  select(
    matches("^Mba_")
  )


# ============================================================
# 6. Check retained expression columns
#
# VERY IMPORTANT for this sensitivity analysis
# ============================================================

colnames(Hab_expr)
colnames(Her_expr)
colnames(Hhe_expr)
colnames(Cin_expr)
colnames(Dda_expr)
colnames(Mbr_expr)


# Confirm that none of the excluded tissue patterns remain
grep(
  "_Mp|_MP|_Ms|_Fa",
  colnames(ortho_table4),
  value = TRUE
)


# All expression columns should be numeric/integer
sapply(Hab_expr, class)
sapply(Her_expr, class)
sapply(Hhe_expr, class)
sapply(Cin_expr, class)
sapply(Dda_expr, class)
sapply(Mbr_expr, class)


# ============================================================
# 7. Calculate testis relative expression
#
# TR = testis expression /
#      total expression across RETAINED tissues
#
# NA in other tissues is ignored in rowSums().
#
# If testis itself is NA:
# TR remains NA.
#
# If total expression = 0:
# TR = NA.
# ============================================================

calc_ratio <- function(x, total) {
  
  ifelse(
    total > 0,
    x / total,
    NA_real_
  )
  
}


# ------------------------------------------------------------
# Total expression across retained tissues
# ------------------------------------------------------------

Hab_total <- rowSums(
  Hab_expr,
  na.rm = TRUE
)

Her_total <- rowSums(
  Her_expr,
  na.rm = TRUE
)

Hhe_total <- rowSums(
  Hhe_expr,
  na.rm = TRUE
)

Cin_total <- rowSums(
  Cin_expr,
  na.rm = TRUE
)

Dda_total <- rowSums(
  Dda_expr,
  na.rm = TRUE
)

Mbr_total <- rowSums(
  Mbr_expr,
  na.rm = TRUE
)


# ------------------------------------------------------------
# Testis relative expression
# ------------------------------------------------------------

ratio_df <- data.frame(
  
  Hab_TR = calc_ratio(
    Hab_expr$Hav_Mt,
    Hab_total
  ),
  
  Her_TR = calc_ratio(
    Her_expr$Hev_Mt,
    Her_total
  ),
  
  Hhe_TR = calc_ratio(
    Hhe_expr$Hha_Mt,
    Hhe_total
  ),
  
  Cin_TR = calc_ratio(
    Cin_expr$Cha_Mt,
    Cin_total
  ),
  
  Dda_TR = calc_ratio(
    Dda_expr$Dda_Mt,
    Dda_total
  ),
  
  Mbr_TR = calc_ratio(
    Mbr_expr$Mba_Mt,
    Mbr_total
  )
  
)


# Add genomic categories
ratio_df$type1 <- ortho_table4$type1
ratio_df$type2 <- ortho_table4$type2
ratio_df$LG_type <- ortho_table4$LG_type

# Mbr gene ID for ASE joining
ratio_df$geneID <- ortho_table4$Mbr_id


head(ratio_df)


# ============================================================
# 8. Define Mbr testis Y-biased genes (TY)
# ============================================================

ASE_data <- read.table(
  "mid_data/allele_specific/HJ/correctBiase_250427/XYbiased_info_geneUnit.v3.rmdup.table",
  header = TRUE,
  sep = "\t"
)


# ------------------------------------------------------------
# Select testis Y-biased genes
# ------------------------------------------------------------

TY_genes <- ASE_data %>%
  
  filter(
    type != "Not signif",
    tissue == "Mba_Mt",
    type == "Ybiased",
    abs(logYXratio) > 0.5
  ) %>%
  
  select(
    geneID,
    type,
    tissue,
    logYXratio
  )


nrow(TY_genes)

head(TY_genes)


# ============================================================
# 9. Join TY information
# ============================================================

plot_data <- ratio_df %>%
  left_join(
    TY_genes,
    by = "geneID"
  )


# ------------------------------------------------------------
# Keep TWO classifications:
#
# type1:
#   SDR / PAR / auto
#
# type4:
#   TY / SDR_others / PAR / auto
#
# This allows Panel A to contain BOTH whole SDR and
# its subdivision into TY + SDR_others.
# ------------------------------------------------------------

plot_data$type4 <- plot_data$type1


# Only SDR genes can become TY
plot_data$type4[
  plot_data$type1 == "SDR" &
    plot_data$type == "Ybiased" &
    plot_data$tissue == "Mba_Mt" &
    abs(plot_data$logYXratio) > 0.5
] <- "TY"


# Remaining SDR genes
plot_data$type4[
  plot_data$type4 == "SDR"
] <- "SDR_others"


table(
  plot_data$type1,
  useNA = "ifany"
)

table(
  plot_data$type4,
  useNA = "ifany"
)


# Sanity check:
# TY + SDR_others should equal whole SDR
sum(plot_data$type1 == "SDR", na.rm = TRUE)

sum(
  plot_data$type4 %in% c("TY", "SDR_others"),
  na.rm = TRUE
)


# ============================================================
# 10. Species order
# ============================================================

species_list <- c(
  "Hab",
  "Her",
  "Hhe",
  "Cin",
  "Dda",
  "Mbr"
)


# ============================================================
# 11. Check species-specific gene universe
# ============================================================

gene_universe_summary <- map_dfr(
  
  species_list,
  
  function(sp) {
    
    tr_col <- paste0(
      sp,
      "_TR"
    )
    
    tibble(
      
      species = sp,
      
      total_genes =
        nrow(plot_data),
      
      genes_with_TR =
        sum(
          !is.na(
            plot_data[[tr_col]]
          )
        ),
      
      missing_TR =
        sum(
          is.na(
            plot_data[[tr_col]]
          )
        )
      
    )
    
  }
  
)


gene_universe_summary


# ============================================================
# 12. Check species-specific top 10%
# ============================================================

top10_summary <- map_dfr(
  
  species_list,
  
  function(sp) {
    
    tr_col <- paste0(
      sp,
      "_TR"
    )
    
    x <- plot_data[[tr_col]]
    
    # Only valid genes for THIS species
    x <- x[
      !is.na(x)
    ]
    
    
    cutoff <- quantile(
      x,
      probs = 0.9,
      na.rm = TRUE
    )
    
    
    tibble(
      
      species = sp,
      
      N = length(x),
      
      cutoff = cutoff,
      
      N_top10 =
        sum(
          x >= cutoff
        ),
      
      proportion =
        sum(
          x >= cutoff
        ) /
        length(x)
      
    )
    
  }
  
)


top10_summary


# ============================================================
# 13. Fisher enrichment function
#
# Each target group is compared with autosome.
# ============================================================

fisher_enrichment <- function(
    data,
    species,
    group_col,
    test_group,
    background_group = "auto",
    cutoff_prob = 0.9
) {
  
  
  # Testis-ratio column
  tr_col <- paste0(
    species,
    "_TR"
  )
  
  
  # ----------------------------------------------------------
  # Species-specific gene universe
  # ----------------------------------------------------------
  
  tmp <- data %>%
    
    filter(
      !is.na(.data[[tr_col]]),
      !is.na(.data[[group_col]])
    )
  
  
  # ----------------------------------------------------------
  # Species-specific 90th percentile
  # ----------------------------------------------------------
  
  cutoff <- quantile(
    tmp[[tr_col]],
    probs = cutoff_prob,
    na.rm = TRUE
  )
  
  
  tmp <- tmp %>%
    
    mutate(
      testis_biased =
        .data[[tr_col]] >= cutoff
    )
  
  
  # ----------------------------------------------------------
  # Counts
  # ----------------------------------------------------------
  
  test_high <- sum(
    tmp[[group_col]] == test_group &
      tmp$testis_biased
  )
  
  
  test_low <- sum(
    tmp[[group_col]] == test_group &
      !tmp$testis_biased
  )
  
  
  auto_high <- sum(
    tmp[[group_col]] == background_group &
      tmp$testis_biased
  )
  
  
  auto_low <- sum(
    tmp[[group_col]] == background_group &
      !tmp$testis_biased
  )
  
  
  # ----------------------------------------------------------
  # Protect against a completely absent target group
  # ----------------------------------------------------------
  
  if (
    (test_high + test_low) == 0 ||
    (auto_high + auto_low) == 0
  ) {
    
    return(
      tibble(
        
        species = species,
        group = test_group,
        cutoff = cutoff,
        
        test_high = test_high,
        test_low = test_low,
        test_total = test_high + test_low,
        
        auto_high = auto_high,
        auto_low = auto_low,
        auto_total = auto_high + auto_low,
        
        testRatio = NA_real_,
        autosomeRatio = NA_real_,
        oddsRatio = NA_real_,
        Pvalue = NA_real_
        
      )
    )
    
  }
  
  
  # ----------------------------------------------------------
  # Fisher 2 x 2 table
  #
  #                  Testis top10     Other
  # target group          a             b
  # autosome              c             d
  # ----------------------------------------------------------
  
  fisher_matrix <- matrix(
    
    c(
      test_high,
      test_low,
      auto_high,
      auto_low
    ),
    
    nrow = 2,
    
    byrow = TRUE,
    
    dimnames = list(
      
      region = c(
        test_group,
        background_group
      ),
      
      expression = c(
        "Testis_top10",
        "Other"
      )
      
    )
    
  )
  
  
  fisher_result <- fisher.test(
    fisher_matrix
  )
  
  
  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------
  
  tibble(
    
    species = species,
    
    group = test_group,
    
    cutoff = cutoff,
    
    
    test_high = test_high,
    
    test_low = test_low,
    
    test_total =
      test_high +
      test_low,
    
    
    auto_high = auto_high,
    
    auto_low = auto_low,
    
    auto_total =
      auto_high +
      auto_low,
    
    
    testRatio =
      test_high /
      (
        test_high +
          test_low
      ),
    
    autosomeRatio =
      auto_high /
      (
        auto_high +
          auto_low
      ),
    
    
    oddsRatio =
      unname(
        fisher_result$estimate
      ),
    
    
    Pvalue =
      fisher_result$p.value
    
  )
  
}


# ============================================================
# 14. PANEL A
#
# Whole SDR + TY subdivision
#
# Rows:
# SDR
# TY
# SDR_others
# PAR
#
# All compared with autosome.
# ============================================================


# ------------------------------------------------------------
# A1. Whole SDR vs autosome
# Uses type1
# ------------------------------------------------------------

A_SDR <- map_dfr(
  
  species_list,
  
  function(sp) {
    
    fisher_enrichment(
      
      data = plot_data,
      
      species = sp,
      
      group_col = "type1",
      
      test_group = "SDR",
      
      background_group = "auto"
      
    )
    
  }
  
)


# ------------------------------------------------------------
# A2. TY / SDR_others / PAR vs autosome
# Uses type4
# ------------------------------------------------------------

A_subgroups <- expand_grid(
  
  species = species_list,
  
  group = c(
    "TY",
    "SDR_others",
    "PAR"
  )
  
) %>%
  
  pmap_dfr(
    
    function(species, group) {
      
      fisher_enrichment(
        
        data = plot_data,
        
        species = species,
        
        group_col = "type4",
        
        test_group = group,
        
        background_group = "auto"
        
      )
      
    }
    
  )


# Combine
final_out1 <- bind_rows(
  A_SDR,
  A_subgroups
)


# BH correction within Panel A
final_out1 <- final_out1 %>%
  
  mutate(
    
    Padj_BH =
      p.adjust(
        Pvalue,
        method = "BH"
      ),
    
    mlogP =
      -log10(
        Pvalue
      ),
    
    mlogPadj_BH =
      -log10(
        Padj_BH
      ),
    
    direction =
      case_when(
        is.na(testRatio) ~ NA_character_,
        testRatio > autosomeRatio ~ "Enriched",
        testRatio < autosomeRatio ~ "Depleted",
        TRUE ~ "Equal"
      ),
    
    species =
      factor(
        species,
        levels = species_list
      ),
    
    # Reverse factor order so SDR appears at the TOP
    group =
      factor(
        group,
        levels = c(
          "PAR",
          "SDR_others",
          "TY",
          "SDR"
        )
      )
    
  )


# ============================================================
# 15. PANEL B
#
# Ancestral LG5 / LG12 / LG13 / LG15
# versus autosome
# ============================================================

group_list2 <- c(
  "5",
  "12",
  "13",
  "15"
)


final_out2 <- expand_grid(
  
  species = species_list,
  
  group = group_list2
  
) %>%
  
  pmap_dfr(
    
    function(species, group) {
      
      fisher_enrichment(
        
        data = plot_data,
        
        species = species,
        
        group_col = "LG_type",
        
        test_group = group,
        
        background_group = "auto"
        
      )
      
    }
    
  )


final_out2 <- final_out2 %>%
  
  mutate(
    
    Padj_BH =
      p.adjust(
        Pvalue,
        method = "BH"
      ),
    
    mlogP =
      -log10(
        Pvalue
      ),
    
    mlogPadj_BH =
      -log10(
        Padj_BH
      ),
    
    direction =
      case_when(
        is.na(testRatio) ~ NA_character_,
        testRatio > autosomeRatio ~ "Enriched",
        testRatio < autosomeRatio ~ "Depleted",
        TRUE ~ "Equal"
      ),
    
    species =
      factor(
        species,
        levels = species_list
      ),
    
    # Desired top-to-bottom order:
    # LG5, LG12, LG13, LG15
    group =
      factor(
        group,
        levels = c(
          "15",
          "13",
          "12",
          "5"
        )
      )
    
  )


# ============================================================
# 16. PANEL C
#
# S3 / S2 / S1 / S0 / PAR
# versus autosome
# ============================================================

group_list3 <- c(
  "S3",
  "S2",
  "S1",
  "S0",
  "PAR"
)


final_out3 <- expand_grid(
  
  species = species_list,
  
  group = group_list3
  
) %>%
  
  pmap_dfr(
    
    function(species, group) {
      
      fisher_enrichment(
        
        data = plot_data,
        
        species = species,
        
        group_col = "type2",
        
        test_group = group,
        
        background_group = "auto"
        
      )
      
    }
    
  )


final_out3 <- final_out3 %>%
  
  mutate(
    
    Padj_BH =
      p.adjust(
        Pvalue,
        method = "BH"
      ),
    
    mlogP =
      -log10(
        Pvalue
      ),
    
    mlogPadj_BH =
      -log10(
        Padj_BH
      ),
    
    direction =
      case_when(
        is.na(testRatio) ~ NA_character_,
        testRatio > autosomeRatio ~ "Enriched",
        testRatio < autosomeRatio ~ "Depleted",
        TRUE ~ "Equal"
      ),
    
    species =
      factor(
        species,
        levels = species_list
      ),
    
    # Desired top-to-bottom order:
    # S3, S2, S1, S0, PAR
    group =
      factor(
        group,
        levels = c(
          "PAR",
          "S0",
          "S1",
          "S2",
          "S3"
        )
      )
    
  )


# ============================================================
# 17. Check result tables
# ============================================================

result_table1 <- final_out1 %>%
  
  select(
    species,
    group,
    test_high,
    test_total,
    auto_high,
    auto_total,
    testRatio,
    autosomeRatio,
    oddsRatio,
    Pvalue,
    Padj_BH,
    direction
  ) %>%
  
  arrange(
    group,
    species
  )


result_table2 <- final_out2 %>%
  
  select(
    species,
    group,
    test_high,
    test_total,
    auto_high,
    auto_total,
    testRatio,
    autosomeRatio,
    oddsRatio,
    Pvalue,
    Padj_BH,
    direction
  ) %>%
  
  arrange(
    group,
    species
  )


result_table3 <- final_out3 %>%
  
  select(
    species,
    group,
    test_high,
    test_total,
    auto_high,
    auto_total,
    testRatio,
    autosomeRatio,
    oddsRatio,
    Pvalue,
    Padj_BH,
    direction
  ) %>%
  
  arrange(
    group,
    species
  )


result_table1
result_table2
result_table3


# ============================================================
# 18. Significant results
# ============================================================

# Raw P < 0.05
final_out1 %>%
  filter(Pvalue < 0.05)

final_out2 %>%
  filter(Pvalue < 0.05)

final_out3 %>%
  filter(Pvalue < 0.05)


# BH-adjusted P < 0.05
final_out1 %>%
  filter(Padj_BH < 0.05)

final_out2 %>%
  filter(Padj_BH < 0.05)

final_out3 %>%
  filter(Padj_BH < 0.05)


# ============================================================
# 19. Rename LG labels for plotting
# ============================================================

final_out2 <- final_out2 %>%
  
  mutate(
    
    group_label = case_when(
      as.character(group) == "5" ~ "LG5",
      as.character(group) == "12" ~ "LG12",
      as.character(group) == "13" ~ "LG13",
      as.character(group) == "15" ~ "LG15"
    ),
    
    group_label = factor(
      group_label,
      levels = c(
        "LG15",
        "LG13",
        "LG12",
        "LG5"
      )
    )
    
  )


# ============================================================
# 20. Plot settings
#
# Size  = % testis-biased genes
# Color = -log10(raw Fisher P)
# *     = raw P < 0.05
# ============================================================


# ------------------------------------------------------------
# Panel A
# ------------------------------------------------------------

p1 <- ggplot(
  
  final_out1,
  
  aes(
    x = species,
    y = group,
    size = testRatio,
    color = mlogP
  )
  
) +
  
  geom_point() +
  
  geom_text(
    data = final_out1 %>%
      filter(
        Pvalue < 0.05
      ),
    aes(
      label = "*"
    ),
    color = "black",
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_size_continuous(
    name = "% Testis-biased genes",
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  scale_color_gradient2(
    low = "#0067a5",
    mid = "#f3f3f4",
    high = "#be0032",
    midpoint = -log10(0.05),
    name = expression(-log[10](P))
  ) +
  
  labs(
    title = "Testis-biased genes enrichment",
    x = NULL,
    y = NULL
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      size = 10,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )


# ------------------------------------------------------------
# Panel B
# ------------------------------------------------------------

p2 <- ggplot(
  
  final_out2,
  
  aes(
    x = species,
    y = group_label,
    size = testRatio,
    color = mlogP
  )
  
) +
  
  geom_point() +
  
  geom_text(
    data = final_out2 %>%
      filter(
        Pvalue < 0.05
      ),
    aes(
      label = "*"
    ),
    color = "black",
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_size_continuous(
    name = "% Testis-biased genes",
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  scale_color_gradient2(
    low = "#0067a5",
    mid = "#f3f3f4",
    high = "#be0032",
    midpoint = -log10(0.05),
    name = expression(-log[10](P))
  ) +
  
  labs(
    title = "Testis-biased genes enrichment",
    x = NULL,
    y = NULL
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      size = 10,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )


# ------------------------------------------------------------
# Panel C
# ------------------------------------------------------------

p3 <- ggplot(
  
  final_out3,
  
  aes(
    x = species,
    y = group,
    size = testRatio,
    color = mlogP
  )
  
) +
  
  geom_point() +
  
  geom_text(
    data = final_out3 %>%
      filter(
        Pvalue < 0.05
      ),
    aes(
      label = "*"
    ),
    color = "black",
    size = 3,
    fontface = "bold",
    show.legend = FALSE
  ) +
  
  scale_size_continuous(
    name = "% Testis-biased genes",
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  
  scale_color_gradient2(
    low = "#0067a5",
    mid = "#f3f3f4",
    high = "#be0032",
    midpoint = -log10(0.05),
    name = expression(-log[10](P))
  ) +
  
  labs(
    title = "Testis-biased genes enrichment",
    x = NULL,
    y = NULL
  ) +
  
  theme_bw() +
  
  theme(
    plot.title = element_text(
      size = 10,
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )


# ============================================================
# 21. Combine panels vertically
#
# A
# B
# C
# ============================================================

p_all <- p1 / p2 / p3

p_all <- p_all +
  
  plot_layout(
    ncol = 1,
    guides = "collect",
    heights = c(1, 1, 1.15)
  ) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    legend.position = "right",
    
    plot.tag = element_text(
      face = "bold",
      size = 12
    )
  )


p_all


# ============================================================
# 22. Save combined figure
# ============================================================

ggsave(
  
  filename =
    "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/testis_Specificity_removePouchTissues.combined.pdf",
  
  plot = p_all,
  
  width = 5,
  
  height = 7
  
)


# ============================================================
# 23. Save statistical tables
# ============================================================

write.table(
  
  final_out1,
  
  "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/testis_SDR_TY_fisher_removePouchTissues.tsv",
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)


write.table(
  
  final_out2,
  
  "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/testis_LG_fisher_removePouchTissues.tsv",
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)


write.table(
  
  final_out3,
  
  "plots/DNA_pipe/hj/hifi/final2/allele_specific/260124/testis_strata_fisher_removePouchTissues.tsv",
  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
  
)

