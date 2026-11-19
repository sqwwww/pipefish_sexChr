
library(tidyverse)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# 1. Read gene information
# ============================================================

my_genes <- readxl::read_xlsx(
  path = "mid_data/DNA_pipe/sex_related_genes/sex_related_genes6.xlsx",
  sheet = 1
)

my_genes$Cin_ID <- my_genes$Cin_ID %>%
  str_split(pattern = ", ") %>%
  map_chr(1) %>%
  str_remove("_Cint")

my_genes$Mbr_ID <- my_genes$Mbr_ID %>%
  str_remove("_HJ")

my_genes$Hab_ID <- my_genes$Hab_ID %>%
  str_remove("_Ha") %>%
  str_split(pattern = "\\.") %>%
  map_chr(1)

my_genes$Dre_ID <- my_genes$Dre_ID %>%
  str_remove("__?zbfe99")

my_genes$Dda_ID <- my_genes$Dda_ID %>%
  str_split(pattern = ", ") %>%
  map_chr(1) %>%
  str_remove("_Ddac")

my_genes$Hhe_ID <- my_genes$Hhe_ID %>%
  str_split(pattern = "\\.") %>%
  map_chr(1)

my_genes$Her_ID <- my_genes$Her_ID %>%
  str_remove("_He")

combined_df <- my_genes

# ============================================================
# 2. Read tissue expression data
# ============================================================

read_tissue_df <- function(file, species_id, remove_gene_prefix = FALSE) {
  
  tmp_data <- read.table(file)
  
  new_df <- tmp_data %>%
    select(-contains("tau"))
  
  if (remove_gene_prefix) {
    new_df$geneID <- rownames(new_df) %>%
      str_remove("gene-") %>%
      str_split("_") %>%
      map_chr(1)
  } else {
    new_df$geneID <- rownames(new_df)
  }
  
  new_df
}

# Cin
Cin_df <- read_tissue_df(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/papa.Fo.tau.tsv",
  "Cin",
  TRUE
)

combined_df2 <- combined_df %>%
  left_join(Cin_df, by = c("Cin_ID" = "geneID"))

# Dda
Dda_df <- read_tissue_df(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/BJ.Fo.tau.tsv",
  "Dda",
  TRUE
)

combined_df3 <- combined_df2 %>%
  left_join(Dda_df, by = c("Dda_ID" = "geneID"))

# Hhe
Hhe_df <- read_tissue_df(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/260529/HF.Fo.tau.tsv",
  "Hhe"
)

combined_df4 <- combined_df3 %>%
  left_join(Hhe_df, by = c("Hhe_ID" = "geneID"))

# Mbr
Mbr_df <- read_tissue_df(
  "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/250325/HJ.Fo.tau.tsv",
  "Mbr"
)

combined_df5 <- combined_df4 %>%
  left_join(Mbr_df, by = c("Mbr_ID" = "geneID"))

# Hab: adult, 3 months, 2 months, 1 month
Hab_files <- c(
  "250325/Ha_Hav.Fo.tau.tsv",
  "Ha_Ha3.Fo.tau.tsv",
  "Ha_Ha2.Fo.tau.tsv",
  "Ha_Ha1.Fo.tau.tsv"
)

base_dir <- "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau/"

combined_df6 <- combined_df5

for (f in Hab_files) {
  
  new_df <- read_tissue_df(
    paste0(base_dir, f),
    "Hab"
  )
  
  combined_df6 <- combined_df6 %>%
    left_join(new_df, by = c("Hab_ID" = "geneID"))
}

# Her: adult, 3 months, 2 months, 1 month
Her_files <- c(
  "250325/He_Hev.Fo.tau.tsv",
  "He_He3.Fo.tau.tsv",
  "He_He2.Fo.tau.tsv",
  "He_He1.Fo.tau.tsv"
)

combined_df7 <- combined_df6

for (f in Her_files) {
  
  new_df <- read_tissue_df(
    paste0(base_dir, f),
    "Her"
  )
  
  combined_df7 <- combined_df7 %>%
    left_join(new_df, by = c("Her_ID" = "geneID"))
}

# ============================================================
# 3. Average Hab and Her expression across ages
# ============================================================

mean_tissue_expression <- function(df, prefix, new_prefix) {
  
  tissue_codes <- c(
    "Fa", "Fb", "Fk", "Fo",
    "Mp", "Mb", "Mk", "Ms", "Mt"
  )
  
  result <- lapply(tissue_codes, function(tissue) {
    
    pattern <- paste0("^", prefix, "._", tissue, "$")
    
    cols <- grep(
      pattern,
      colnames(df),
      value = TRUE,
      ignore.case = TRUE
    )
    
    if (length(cols) == 0) {
      return(rep(NA_real_, nrow(df)))
    }
    
    x <- df[, cols, drop = FALSE] %>%
      mutate(across(everything(), as.numeric))
    
    y <- rowMeans(x, na.rm = TRUE)
    y[is.nan(y)] <- NA_real_
    
    y
  })
  
  result <- as.data.frame(result)
  
  colnames(result) <- paste0(
    new_prefix, "_", tissue_codes
  )
  
  rownames(result) <- combined_df7$geneName
  
  result
}

Ha_mean_df <- mean_tissue_expression(
  combined_df7, "Ha", "Ham"
)

He_mean_df <- mean_tissue_expression(
  combined_df7, "He", "Hem"
)

# Combine averaged Hab/Her and other species
other_df <- combined_df7 %>%
  select(-matches("^H[ae]._")) %>%
  select(-(1:9))

combined_df9 <- cbind(
  Ha_mean_df,
  He_mean_df,
  other_df
)

# Keep six tissues:
# Mt Testis, Fo Ovary, Mk Male kidney,
# Fk Female kidney, Mb Male brain, Fb Female brain

combined_df9 <- combined_df9 %>%
  select(matches("_(Mt|Fo|Mk|Fk|Mb|Fb)$"))

# ============================================================
# 4. Select genes
# ============================================================

selected_genes <- c(
  "dmrt1",
  "ar",
  "gsdf1",
  "gsdf2",
  "cyp11a",
  "amh"
)

expression_df <- combined_df9 %>%
  mutate(
    across(
      everything(),
      ~ suppressWarnings(as.numeric(as.character(.x)))
    )
  )

missing_genes <- setdiff(
  selected_genes,
  rownames(expression_df)
)

if (length(missing_genes) > 0) {
  warning(
    "Missing genes: ",
    paste(missing_genes, collapse = ", ")
  )
}

genes_to_plot <- selected_genes[
  selected_genes %in% rownames(expression_df)
]

if (length(genes_to_plot) == 0) {
  stop("No selected genes found.")
}

selected_df <- as.matrix(
  expression_df[
    genes_to_plot,
    ,
    drop = FALSE
  ]
)

selected_df[!is.finite(selected_df)] <- NA_real_

# ============================================================
# 5. Raw TPM matrix
# No log2 transformation
# No row Z-score
# No normalization
# ============================================================

heatmap_mat <- selected_df

# Check raw TPM distribution
print(summary(as.vector(heatmap_mat)))

print(
  quantile(
    as.vector(heatmap_mat),
    probs = c(
      0, 0.25, 0.5, 0.75,
      0.9, 0.95, 0.99, 1
    ),
    na.rm = TRUE
  )
)

# ============================================================
# 6. Species and tissue information
# ============================================================

column_species <- case_when(
  str_starts(colnames(heatmap_mat), "Ham_") ~ "Ha",
  str_starts(colnames(heatmap_mat), "Hem_") ~ "He",
  str_starts(colnames(heatmap_mat), "Hha_") ~ "Hh",
  str_starts(colnames(heatmap_mat), "Cha_") ~ "Ch",
  str_starts(colnames(heatmap_mat), "Dda_") ~ "Dd",
  str_starts(colnames(heatmap_mat), "Mba_") ~ "Mb",
  TRUE ~ "Other"
)

species_order <- c(
  "Ha", "He", "Hh", "Ch", "Dd", "Mb", "Other"
)

column_species <- factor(
  column_species,
  levels = species_order
)

tissue_code <- str_extract(
  colnames(heatmap_mat),
  "(Fb|Fk|Fo|Mb|Mk|Mt)$"
)

tissue_order <- c(
  "Mt", "Fo", "Mk", "Fk", "Mb", "Fb"
)

column_tissue_code <- factor(
  tissue_code,
  levels = tissue_order
)

column_order <- order(
  column_species,
  column_tissue_code
)

heatmap_mat <- heatmap_mat[
  ,
  column_order,
  drop = FALSE
]

column_species <- droplevels(
  column_species[column_order]
)

tissue_code <- tissue_code[column_order]

tissue_names <- c(
  Mt = "Testis",
  Fo = "Ovary",
  Mk = "Male kidney",
  Fk = "Female kidney",
  Mb = "Male brain",
  Fb = "Female brain"
)

column_tissue <- factor(
  unname(tissue_names[tissue_code]),
  levels = unname(tissue_names)
)

tissue_labels <- as.character(column_tissue)

# ============================================================
# 7. Top annotations
# ============================================================

species_colours <- c(
  Ha = "#5B8FF9",
  He = "#61DDAA",
  Hh = "#65789B",
  Ch = "#F6BD16",
  Dd = "#7262FD",
  Mb = "#E8684A",
  Other = "#A0A0A0"
)

tissue_colours <- c(
  "Female brain"  = "#9ECAE1",
  "Female kidney" = "#6BAED6",
  "Ovary"         = "#E78AC3",
  "Male brain"    = "#A1D99B",
  "Male kidney"   = "#31A354",
  "Testis"        = "#F46D43"
)

top_anno <- HeatmapAnnotation(
  Species = column_species,
  Tissue = column_tissue,
  col = list(
    Species = species_colours,
    Tissue = tissue_colours
  ),
  annotation_name_side = "left",
  show_legend = TRUE
)

# ============================================================
# 8. Raw TPM colour scale
# ============================================================

# Adjust these breaks after inspecting TPM quantiles.
# Values above 50 TPM will have the same maximum colour.

col_fun <- colorRamp2(
  c(0, 5, 50),
  c(
    "#0067a5",
    "#f3f3f4",
    "#be0032"
  )
)

# ============================================================
# 9. Draw heatmap
# ============================================================

sc_select <- Heatmap(
  heatmap_mat,
  
  name = "TPM",
  col = col_fun,
  na_col = "grey90",
  
  top_annotation = top_anno,
  
  column_split = column_species,
  column_labels = tissue_labels,
  
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "left",
  row_names_gp = gpar(
    fontface = "italic",
    fontsize = 10
  ),
  
  column_names_rot = 60,
  column_names_gp = gpar(fontsize = 8),
  
  column_gap = unit(2.5, "mm"),
  
  rect_gp = gpar(
    col = "white",
    lwd = 0.4
  ),
  
  heatmap_legend_param = list(
    title = "TPM",
    at = c(0, 5, 10, 20, 50)
  )
)

# Display heatmap
draw(
  sc_select,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# ============================================================
# 10. Save PDF
# ============================================================

output_dir <- "plots/DNA_pipe/sex_related_genes"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

pdf(
  file.path(
    output_dir,
    "selected_genes_all_tissues_heatmap.v4.pdf"
  ),
  width = max(
    10,
    ncol(heatmap_mat) * 0.22
  ),
  height = max(
    4.5,
    nrow(heatmap_mat) * 0.42 + 2.5
  ),
  useDingbats = FALSE
)

draw(
  sc_select,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

# # ============================================================
# # 11. Save raw TPM matrix
# # ============================================================
# 
# write.csv(
#   heatmap_mat,
#   file.path(
#     output_dir,
#     "selected_genes_all_tissues_heatmap_raw_TPM_matrix.csv"
#   ),
#   quote = FALSE
# )
