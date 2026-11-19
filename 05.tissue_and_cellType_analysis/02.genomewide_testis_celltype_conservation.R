# Complex turnover subtypes are pooled for plotting and category counts.
# Exact cell-type combinations remain in the solutions column.
# Complex and ambiguous turnover are pooled. Other is an informative cell type.
# Conservation and turnover both require at least three informative species.
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(grid)
})

outdir <- "C:/Users/suqianwei/Documents/project/seahorse/seahorse/plots/DNA_pipe/sex_related_genes/all_gene_cellType/separate_germ_subtypes"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
ortho_file <- paste0("mid_data/RNA_counts/tissue_specific_analysis/",
                     "02.analysis_241013/01.orthoFinder_result/OrthoTable_SingleCopy.table")
dre_idconvert_file <- "../ref/zbf/Lawson_zebrafish_v4_3_2geneinfo.txt"
idconvert_file <- "../seahorse_sc/ref/ortho_zbfens99_Ha_Cin_Mbr/Ha_Cint_Mbr_zbfens99_singleCopy.convert3"
sc_files <- c(
  Hab = "../seahorse_sc/mid_data/combined/pseudobulk/Hab_rmBias_pseudobulk.v2.tsv",
  Cin = "../seahorse_sc/mid_data/combined/pseudobulk/Cin_rmBias_pseudobulk.v3.tsv",
  Mbr = "../seahorse_sc/mid_data/combined/pseudobulk/Mbr_rmBias_pseudobulk.v2.tsv",
  Dre = "../seahorse_sc/mid_data/combined/pseudobulk/zbf_rmBias_pseudobulk.v2.tsv"
)
species_order <- names(sc_files)
min_conserved_species <- 3L
min_turnover_species <- 3L
sc_low_expr_quantile <- 0.50
sc_top1_top2_cutoff <- 1.5
# Group by biological subtype, then order similar categorical patterns together.
# This avoids assigning arbitrary numeric distances to nominal cell types.
# No pagination: the main PDF is one overview page.

save_tsv <- function(x, name) write_tsv(x, file.path(outdir, name), na = "NA")
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
clean_id <- function(x) {
  x <- str_trim(as.character(x))
  x <- str_remove(x, "^gene-")
  x <- str_remove(x, "_{1,2}(Cint|Cin|HJ|Mbr|Ha|Hab|Ddac|Dda|He|Her|HF|Hhe|zbfe99)$")
  x <- str_remove(x, "\\.t[0-9]+$")
  x <- str_replace(x, "^(ENS[A-Z]*G[0-9]+)\\.[0-9]+$", "\\1")
  x <- str_replace(x, "^(HaiMa[0-9]+)\\.[0-9]+$", "\\1")
  x[x %in% c("", "NA", "-", ".")] <- NA_character_
  x
}
missing_files <- c(ortho_file, idconvert_file, sc_files, dre_idconvert_file)
missing_files <- missing_files[!file.exists(missing_files)]
if (length(missing_files)) stop("Missing input files:\n", paste(missing_files, collapse = "\n"))

ortho <- read.table(ortho_file)
if (ncol(ortho) != 7L) stop("Expected 7 columns in OrthoTable_SingleCopy.table")
names(ortho) <- c("Orthogroup", "Cin", "Dda", "Hhe", "Mbr", "Hab", "Her")
ortho <- as_tibble(ortho) %>% mutate(across(everything(), as.character))
if (anyDuplicated(ortho$Orthogroup) || anyNA(ortho$Orthogroup)) stop("Invalid orthogroup IDs")
if (any(str_detect(unlist(ortho[-1]), ","), na.rm = TRUE)) stop("Expected single-copy IDs, found comma-separated copies")
all_genes <- ortho$Orthogroup
native_map <- ortho %>% select(Orthogroup, Hab, Cin, Mbr) %>%
  pivot_longer(-Orthogroup, names_to = "species", values_to = "gene_id") %>%
  mutate(gene_id = clean_id(gene_id))

# Preserve the original Mbr -> LLgeneSymbol bridge for ALL THREE syngnathids.
# Also retain native-ID candidates for features that were not renamed.
conv <- read.delim(idconvert_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("Mbr", "LLgeneSymbol") %in% names(conv))) stop("Conversion table requires Mbr and LLgeneSymbol")
conv_raw <- as_tibble(conv)
conv <- conv_raw %>% transmute(Mbr_key = clean_id(Mbr), symbol = clean_id(LLgeneSymbol)) %>%
  filter(!is.na(Mbr_key), !is.na(symbol)) %>% distinct()
bridge <- ortho %>% transmute(Orthogroup, Mbr_key = clean_id(Mbr)) %>%
  inner_join(conv, by = "Mbr_key") %>% select(Orthogroup, symbol)
ambiguous_bridge <- bridge %>% count(Orthogroup) %>% filter(n > 1)
save_tsv(ambiguous_bridge, "01_ambiguous_symbol_bridge.tsv")
bridge <- bridge %>% anti_join(ambiguous_bridge, by = "Orthogroup")
symbol_map <- crossing(bridge, species = c("Hab", "Cin", "Mbr")) %>%
  rename(gene_id = symbol)
# Original verified gsdf corrections: Hab gsdf refers to HaiMa04881;
# Mbr g4136 may have the feature name gsdf.
extra_map <- native_map %>% filter(species == "Mbr", gene_id == "g4136") %>% mutate(gene_id = "gsdf")
gene_map <- bind_rows(native_map, symbol_map, extra_map) %>%
  mutate(gene_id = if_else(species == "Hab" & gene_id == "gsdf", "HaiMa04881", gene_id))
# Add zebrafish using the existing orthology bridge, NOT a new symbol-based
# orthology inference. LLgeneSymbol must refer to the zebrafish ortholog in
# Ha_Cint_Mbr_zbfens99_singleCopy.convert3 (as in the original scRNA workflow).
# If that table contains an ENSDARG ID column, use it directly in preference
# to the symbol bridge. Otherwise resolve its zebrafish symbols with Lawson.
dre_info <- read.delim(dre_idconvert_file, check.names = FALSE,
                       stringsAsFactors = FALSE)
if (!all(c("Ens99geneIDversion", "LLgeneSymbol") %in% names(dre_info)))
  stop("Dre gene-info table requires Ens99geneIDversion and LLgeneSymbol")
dre_info <- as_tibble(dre_info) %>%
  transmute(dre_id = clean_id(Ens99geneIDversion), symbol = clean_id(LLgeneSymbol)) %>%
  filter(!is.na(dre_id), !is.na(symbol)) %>% distinct()
# Symbols associated with several distinct zebrafish genes are not safe aliases.
ambiguous_dre_symbols <- dre_info %>% count(symbol, name = "n_Dre_IDs") %>% filter(n_Dre_IDs > 1)
save_tsv(ambiguous_dre_symbols, "00_Dre_ambiguous_symbols.tsv")
dre_aliases <- dre_info %>% anti_join(ambiguous_dre_symbols, by = "symbol")
dre_id_cols <- names(conv_raw)[vapply(conv_raw, function(x)
  any(str_detect(as.character(x), "ENSDARG[0-9]+"), na.rm = TRUE), logical(1))]
if (length(dre_id_cols)) {
  # Restrict to simple single-gene entries; exclude multi-copy entries.
  dre_links <- conv_raw %>% transmute(Mbr_key = clean_id(Mbr),
                                      across(all_of(dre_id_cols), as.character)) %>%
    pivot_longer(-Mbr_key, names_to = "source_column", values_to = "dre_id") %>%
    mutate(dre_id = clean_id(dre_id)) %>%
    filter(!is.na(Mbr_key), str_detect(dre_id, "^ENSDARG[0-9]+$")) %>%
    select(Mbr_key, dre_id) %>% distinct()
  dre_og <- ortho %>% transmute(Orthogroup, Mbr_key = clean_id(Mbr)) %>%
    inner_join(dre_links, by = "Mbr_key") %>% select(Orthogroup, dre_id) %>% distinct()
  dre_mapping_source <- paste("Direct ENSDARG columns:", paste(dre_id_cols, collapse = ", "))
} else {
  dre_og <- bridge %>% inner_join(dre_aliases, by = "symbol") %>%
    select(Orthogroup, dre_id) %>% distinct()
  dre_mapping_source <- "Existing single-copy orthology LLgeneSymbol bridge + Lawson"
}
# Enforce unique OG <-> Dre gene associations among the mapped candidates.
bad_og <- dre_og %>% count(Orthogroup, name = "n_Dre_IDs") %>% filter(n_Dre_IDs > 1)
bad_dre <- dre_og %>% count(dre_id, name = "n_OGs") %>% filter(n_OGs > 1)
save_tsv(bad_og, "00_Dre_ambiguous_orthogroups.tsv")
save_tsv(bad_dre, "00_Dre_genes_in_multiple_orthogroups.tsv")
dre_og <- dre_og %>% anti_join(bad_og, by = "Orthogroup") %>% anti_join(bad_dre, by = "dre_id")
if (!nrow(dre_og)) stop("No unique Dre orthologs resolved; inspect conversion-table columns and IDs")
save_tsv(dre_og, "00_Dre_orthogroup_gene_mapping.tsv")
# Accept either original Ensembl IDs or renamed LLgeneSymbol pseudobulk rows.
dre_map <- bind_rows(
  dre_og %>% transmute(Orthogroup, species = "Dre", gene_id = dre_id),
  dre_og %>% inner_join(dre_aliases, by = "dre_id") %>%
    transmute(Orthogroup, species = "Dre", gene_id = symbol)
) %>% distinct()
gene_map <- bind_rows(gene_map, dre_map)
message("Dre mapping source: ", dre_mapping_source,
        "; uniquely linked OGs: ", nrow(dre_og))
gene_map <- gene_map %>% filter(!is.na(gene_id)) %>% distinct()
# A pseudobulk feature must not be copied into several orthogroups.
ambiguous_features <- gene_map %>% group_by(species, gene_id) %>%
  summarise(n_OG = n_distinct(Orthogroup), .groups = "drop") %>% filter(n_OG > 1)
save_tsv(ambiguous_features, "02_ambiguous_features_excluded.tsv")
gene_map <- gene_map %>% anti_join(ambiguous_features, by = c("species", "gene_id"))
save_tsv(gene_map, "03_gene_ID_mapping.tsv")

make_long_pb <- function(path, sp) {
  # Same reader convention as the original working script: row names are IDs.
  d <- read.table(path, sep = "\t", check.names = FALSE)
  if (identical(rownames(d), as.character(seq_len(nrow(d)))))
    stop("No gene row names detected in ", path, "; inspect the pseudobulk header")
  if (!all(vapply(d, is.numeric, logical(1)))) stop("Non-numeric pseudobulk columns in ", path)
  as_tibble(d, rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "celltype", values_to = "expr") %>%
    mutate(species = sp, gene_id = clean_id(gene_id),
           gene_id = if_else(sp == "Hab" & gene_id == "gsdf", "HaiMa04881", gene_id))
}
pb_long <- map2_dfr(unname(sc_files), names(sc_files), make_long_pb)
if (any(!is.finite(pb_long$expr) & !is.na(pb_long$expr)) || any(pb_long$expr < 0, na.rm = TRUE))
  stop("Expected finite nonnegative pseudobulk expression or NA")
# One threshold per species, calculated BEFORE ortholog selection.
gene_maxima <- pb_long %>% group_by(species, gene_id) %>%
  summarise(gene_max = safe_max(expr), .groups = "drop")
threshold_df <- gene_maxima %>% group_by(species) %>%
  summarise(threshold = quantile(gene_max, sc_low_expr_quantile, na.rm = TRUE, names = FALSE),
            n_genes_for_threshold = sum(!is.na(gene_max)), .groups = "drop")
if (any(!is.finite(threshold_df$threshold))) stop("Undefined low-expression threshold")
threshold_qc <- gene_maxima %>% left_join(threshold_df, by = "species") %>%
  group_by(species, threshold, n_genes_for_threshold) %>%
  summarise(n_below = sum(gene_max < threshold, na.rm = TRUE),
            n_equal = sum(gene_max == threshold, na.rm = TRUE),
            n_above = sum(gene_max > threshold, na.rm = TRUE),
            n_zero = sum(gene_max == 0, na.rm = TRUE), .groups = "drop")
save_tsv(threshold_qc, "04_threshold_QC_all_pseudobulk_genes.tsv")
print(threshold_qc)
standardize_celltype <- function(x) {
  # Merge only early/late spermatids; retain other cell types separately.
  x <- str_to_lower(str_trim(x))
  if_else(x %in% c("early.spermatid", "late.spermatid"),
          "spermatid", x)
}
save_tsv(tibble(original = unique(pb_long$celltype),
                standardized = standardize_celltype(unique(pb_long$celltype))), "05_celltype_labels.tsv")
pb_mapped <- pb_long %>% inner_join(gene_map, by = c("species", "gene_id")) %>%
  mutate(celltype = standardize_celltype(celltype)) %>%
  group_by(Orthogroup, species, celltype) %>% summarise(expr = safe_max(expr), .groups = "drop")
mapping_qc <- pb_mapped %>% filter(!is.na(expr)) %>% distinct(species, Orthogroup) %>%
  count(species, name = "n_mapped_ogs") %>%
  complete(species = species_order, fill = list(n_mapped_ogs = 0)) %>%
  mutate(total_ogs = length(all_genes), fraction_mapped = n_mapped_ogs / total_ogs)
save_tsv(mapping_qc, "06_mapping_QC.tsv")
print(mapping_qc)
if (any(mapping_qc$n_mapped_ogs == 0)) stop("Zero mapped OGs: inspect ID mapping and feature names")
if (any(mapping_qc$fraction_mapped < 0.5)) warning("Some species have <50% mapping; inspect 06_mapping_QC.tsv before interpretation")

rank_types <- function(d) {
  d <- d %>% filter(!is.na(expr)) %>% arrange(desc(expr), celltype)
  tibble(top1 = if (nrow(d)) d$celltype[1] else NA_character_,
         top2 = if (nrow(d) >= 2) d$celltype[2] else NA_character_,
         e1 = if (nrow(d)) d$expr[1] else NA_real_,
         e2 = if (nrow(d) >= 2) d$expr[2] else NA_real_)
}
ranked <- pb_mapped %>% group_by(Orthogroup, species) %>%
  group_modify(~rank_types(.x)) %>% ungroup()
sc_states <- crossing(Orthogroup = all_genes, species = species_order) %>%
  left_join(ranked, by = c("Orthogroup", "species")) %>%
  left_join(threshold_df, by = "species") %>% mutate(
    # >= threshold reproduces the original boundary; zero is never informative.
    valid1 = !is.na(e1) & e1 > 0 & e1 >= threshold,
    valid2 = !is.na(e2) & e2 > 0 & e2 >= threshold,
    ratio = if_else(valid2, e1 / e2, Inf),
    celltype_class = case_when(is.na(e1) ~ "missing", !valid1 ~ "low_expr",
                               valid2 & ratio <= sc_top1_top2_cutoff ~ "two_celltypes",
                               TRUE ~ "one_celltype"),
    plot_celltype1 = case_when(celltype_class == "missing" ~ "missing",
                               celltype_class == "low_expr" ~ "low_expr", TRUE ~ top1),
    plot_celltype2 = if_else(celltype_class == "two_celltypes", top2, NA_character_)) %>%
  rowwise() %>% mutate(celltype_set = list({
    if (celltype_class %in% c("one_celltype", "two_celltypes"))
      sort(unique(na.omit(c(plot_celltype1, plot_celltype2)))) else character(0)
  }), informative = length(celltype_set) > 0) %>% ungroup()
save_tsv(sc_states %>% select(-celltype_set), "07_species_celltype_states.tsv")
state_qc <- sc_states %>% count(species, celltype_class, informative, name = "n_orthogroups")
save_tsv(state_qc, "08_state_counts.tsv")
print(state_qc, n = Inf)

# This is the original minimum-distinct-state rule, not a phylogenetic
# ancestral-state reconstruction. Turnover pairs are UNDIRECTED.
solve_sets <- function(sets) {
  sets <- sets[lengths(sets) > 0]
  n <- length(sets)
  if (!n) return(tibble(n_valid_species = 0L, minimum_distinct_celltypes = 0L,
                        n_optimal_solutions = 0L, shared = "none", solutions = "none"))
  grid <- do.call(expand.grid, c(sets, stringsAsFactors = FALSE))
  ss <- lapply(seq_len(nrow(grid)), function(i) sort(unique(as.character(unlist(grid[i, ])))))
  k <- min(lengths(ss))
  opt <- sort(unique(vapply(ss[lengths(ss) == k], paste, collapse = " <-> ", FUN.VALUE = character(1))))
  shared <- Reduce(intersect, sets)
  tibble(n_valid_species = n, minimum_distinct_celltypes = k,
         n_optimal_solutions = length(opt), shared = if (length(shared)) paste(sort(shared), collapse = " + ") else "none",
         solutions = paste(opt, collapse = " / "))
}
celltype_gene_status <- sc_states %>% group_by(Orthogroup) %>%
  group_modify(~solve_sets(.x$celltype_set)) %>% ungroup() %>% mutate(
    status = case_when(
      n_valid_species >= min_conserved_species & minimum_distinct_celltypes == 1 ~ "Conserved",
      n_valid_species >= min_turnover_species & minimum_distinct_celltypes > 1 ~ "Turnover",
      TRUE ~ "Insufficient data"),
    subtype = case_when(status == "Insufficient data" ~ "Insufficient data",
                        status == "Conserved" ~ paste0("Conserved: ", shared),
                        n_optimal_solutions > 1 ~ "Complex/ambiguous turnover",
                        minimum_distinct_celltypes == 2 ~ paste0("Turnover: ", solutions),
                        TRUE ~ "Complex/ambiguous turnover"))
species_coverage <- sc_states %>% group_by(Orthogroup) %>%
  summarise(informative_species = paste(species[informative], collapse = ";"),
            Dre_informative = any(species == "Dre" & informative), .groups = "drop")
celltype_gene_status <- celltype_gene_status %>% left_join(species_coverage, by = "Orthogroup")
save_tsv(celltype_gene_status, "09_all_orthogroups_classification.tsv")
counts <- celltype_gene_status %>% count(status, subtype, name = "n_orthogroups") %>%
  mutate(percent_all = 100 * n_orthogroups / length(all_genes))
save_tsv(counts, "10_category_counts.tsv")
save_tsv(celltype_gene_status %>% count(n_valid_species, name = "n_orthogroups"), "11_valid_species_counts.tsv")
print(counts, n = Inf)


# ---------------- Enrichr symbol lists: both cell-type resolutions ----------------
# Recompute the grouped analysis BEFORE ranking; relabeling final states is not equivalent.
classify_grouped_for_GO <- function(pb_mapped) {
  pb_mapped <- pb_mapped %>%
    mutate(celltype = if_else(celltype %in% c("spermatogonia", "spermatocyte", "spermatid"),
                              "germ cell", celltype)) %>%
    group_by(Orthogroup, species, celltype) %>%
    summarise(expr = safe_max(expr), .groups = "drop")
  ranked <- pb_mapped %>% group_by(Orthogroup, species) %>%
    group_modify(~rank_types(.x)) %>% ungroup()
  sc_states <- crossing(Orthogroup = all_genes, species = species_order) %>%
    left_join(ranked, by = c("Orthogroup", "species")) %>%
    left_join(threshold_df, by = "species") %>% mutate(
      # >= threshold reproduces the original boundary; zero is never informative.
      valid1 = !is.na(e1) & e1 > 0 & e1 >= threshold,
      valid2 = !is.na(e2) & e2 > 0 & e2 >= threshold,
      ratio = if_else(valid2, e1 / e2, Inf),
      celltype_class = case_when(is.na(e1) ~ "missing", !valid1 ~ "low_expr",
                                 valid2 & ratio <= sc_top1_top2_cutoff ~ "two_celltypes",
                                 TRUE ~ "one_celltype"),
      plot_celltype1 = case_when(celltype_class == "missing" ~ "missing",
                                 celltype_class == "low_expr" ~ "low_expr", TRUE ~ top1),
      plot_celltype2 = if_else(celltype_class == "two_celltypes", top2, NA_character_)) %>%
    rowwise() %>% mutate(celltype_set = list({
      if (celltype_class %in% c("one_celltype", "two_celltypes"))
        sort(unique(na.omit(c(plot_celltype1, plot_celltype2)))) else character(0)
    }), informative = length(celltype_set) > 0) %>% ungroup()
  celltype_gene_status <- sc_states %>% group_by(Orthogroup) %>%
    group_modify(~solve_sets(.x$celltype_set)) %>% ungroup() %>% mutate(
      status = case_when(
        n_valid_species >= min_conserved_species & minimum_distinct_celltypes == 1 ~ "Conserved",
        n_valid_species >= min_turnover_species & minimum_distinct_celltypes > 1 ~ "Turnover",
        TRUE ~ "Insufficient data"),
      subtype = case_when(status == "Insufficient data" ~ "Insufficient data",
                          status == "Conserved" ~ paste0("Conserved: ", shared),
                          n_optimal_solutions > 1 ~ "Complex/ambiguous turnover",
                          minimum_distinct_celltypes == 2 ~ paste0("Turnover: ", solutions),
                          TRUE ~ "Complex/ambiguous turnover"))
  
  celltype_gene_status
}
grouped_GO_status <- classify_grouped_for_GO(pb_mapped)
go_dir <- file.path(outdir, "GO_Enrichr")
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(grouped_GO_status, file.path(go_dir, "grouped_classification.tsv"))

# Use the existing zebrafish LLgeneSymbol bridge, without changing case or
# inferring human orthologs. Missing/ambiguous symbols remain visible in detail tables.
go_symbol_map <- bridge %>%
  transmute(Orthogroup, symbol = str_trim(as.character(symbol))) %>%
  filter(!is.na(symbol), !symbol %in% c("", "NA", "-", ".")) %>% distinct()

go_groups <- list(
  grouped_conserved_germ_cell = filter(grouped_GO_status, subtype == "Conserved: germ cell"),
  grouped_conserved_sertoli = filter(grouped_GO_status, subtype == "Conserved: sertoli"),
  grouped_turnover_germ_cell_sertoli = filter(grouped_GO_status,
                                              status == "Turnover", minimum_distinct_celltypes == 2, n_optimal_solutions == 1,
                                              solutions %in% c("germ cell <-> sertoli", "sertoli <-> germ cell")),
  stage_conserved_all = filter(celltype_gene_status, status == "Conserved"),
  stage_turnover_all = filter(celltype_gene_status, status == "Turnover"),
  stage_conserved_spermatogonia = filter(celltype_gene_status, subtype == "Conserved: spermatogonia"),
  stage_conserved_spermatocyte = filter(celltype_gene_status, subtype == "Conserved: spermatocyte"),
  stage_conserved_spermatid = filter(celltype_gene_status, subtype == "Conserved: spermatid"),
  stage_conserved_sertoli = filter(celltype_gene_status, subtype == "Conserved: sertoli"),
  background_grouped_comparable = filter(grouped_GO_status, status != "Insufficient data"),
  background_stage_comparable = filter(celltype_gene_status, status != "Insufficient data")
)
go_export_summary <- imap_dfr(go_groups, function(d, label) {
  detail <- d %>% left_join(go_symbol_map, by = "Orthogroup") %>% arrange(Orthogroup, symbol)
  symbols <- sort(unique(detail$symbol[!is.na(detail$symbol)]))
  # One unique symbol per line; no header, quotes or row numbers.
  write_lines(symbols, file.path(go_dir, paste0(label, ".txt")))
  write_tsv(detail, file.path(go_dir, paste0(label, "_details.tsv")), na = "NA")
  tibble(group = label, n_orthogroups = n_distinct(d$Orthogroup),
         n_orthogroups_with_symbol = n_distinct(detail$Orthogroup[!is.na(detail$symbol)]),
         n_orthogroups_without_symbol = n_distinct(detail$Orthogroup[is.na(detail$symbol)]),
         n_unique_symbols = length(symbols))
})
# Reference counts are QC expectations, never used to select or force genes.
go_expected <- tibble(
  group = c("grouped_conserved_germ_cell", "grouped_conserved_sertoli",
            "grouped_turnover_germ_cell_sertoli", "stage_conserved_all", "stage_turnover_all"),
  reference_n = c(3669L, 155L, 434L, 1706L, 4767L))
go_export_summary <- go_export_summary %>% left_join(go_expected, by = "group") %>%
  mutate(matches_reference = if_else(is.na(reference_n), NA, n_orthogroups == reference_n))
write_tsv(go_export_summary, file.path(go_dir, "00_export_summary.tsv"))
print(go_export_summary, n = Inf, width = Inf)
if (any(go_export_summary$matches_reference == FALSE, na.rm = TRUE))
  warning("GO group counts differ from supplied reference counts; inspect 00_export_summary.tsv")
write_lines(c(
  "TXT files: unique zebrafish LLgeneSymbol entries, one per line, no header.",
  "Use an enrichment service/library that supports these identifiers and organism.",
  "Human symbols have NOT been inferred; uppercasing is not orthology mapping.",
  "grouped_*: all germ-cell stages merged by maximum expression BEFORE top-two ranking.",
  "stage_*: spermatogonia/spermatocyte/spermatid separate; early/late spermatids merged.",
  "Single conserved-cell-type lists exclude categories with multiple shared cell types.",
  "stage_turnover_all includes complex/ambiguous turnover.",
  "background_* lists contain all comparable genes for the corresponding resolution.",
  "Use the matching background where the chosen enrichment interface supports a custom background.",
  "Details retain unmapped genes; summary reports OG counts and unique symbol counts separately."
), file.path(go_dir, "README.txt"))
# ---------------- end Enrichr export ----------------

# Only comparable genes appear in the overview; all OGs remain in exported tables.
plot_info <- celltype_gene_status %>% filter(status != "Insufficient data")
if (!nrow(plot_info)) stop("No OGs with enough informative species; QC tables have been saved")
celltype_plot_long <- sc_states %>% filter(Orthogroup %in% plot_info$Orthogroup)
# Canonical categorical pattern: two-cell states are unordered for row sorting.
patterns <- celltype_plot_long %>% rowwise() %>%
  mutate(pattern = paste(sort(na.omit(c(plot_celltype1, plot_celltype2))), collapse = "+")) %>%
  ungroup() %>% select(Orthogroup, species, pattern) %>%
  pivot_wider(names_from = species, values_from = pattern)
subtype_levels <- plot_info %>% distinct(status, subtype) %>%
  mutate(rank = case_when(status == "Conserved" ~ 1L,
                          str_starts(subtype, "Turnover:") ~ 2L,
                          str_starts(subtype, "Complex") ~ 3L, TRUE ~ 4L)) %>%
  arrange(rank, subtype) %>% pull(subtype)
plot_info <- plot_info %>% left_join(patterns, by = "Orthogroup") %>%
  mutate(subtype = factor(subtype, levels = subtype_levels)) %>%
  arrange(subtype, across(all_of(species_order)), Orthogroup)
row_ids <- plot_info$Orthogroup
save_tsv(plot_info %>% mutate(plot_row = row_number()), "12_heatmap_row_order.tsv")
make_matrix <- function(field) {
  w <- celltype_plot_long %>% select(Orthogroup, species, all_of(field)) %>%
    pivot_wider(names_from = species, values_from = all_of(field)) %>% column_to_rownames("Orthogroup")
  as.matrix(w[row_ids, species_order, drop = FALSE])
}
mat1 <- make_matrix("plot_celltype1")
mat2 <- make_matrix("plot_celltype2")
celltype_col <- c(spermatogonia = "#0067a5", spermatocyte = "#f38400",
                  spermatid = "#008856", sertoli = "#E68FAC", leydig = "#B3446C",
                  immune = "#654522", erythroid = "#E9BD1B", other = "#BDBDBD",
                  low_expr = "#E5E5E5", missing = "#FFFFFF")
extra_types <- setdiff(na.omit(unique(c(mat1, mat2))), names(celltype_col))
if (length(extra_types)) celltype_col <- c(celltype_col, setNames(grDevices::hcl.colors(length(extra_types), "Dark 3"), extra_types))
# Two adjacent subcolumns per species preserve both retained major cell types,
# including in a rasterized overview containing thousands of rows.
# Single-type / low-expression / missing cells repeat their color in both halves.
plot_mat <- matrix(NA_character_, nrow(mat1), 2 * ncol(mat1),
                   dimnames = list(row_ids, as.vector(rbind(paste0(species_order, "_1"), paste0(species_order, "_2")))))
plot_mat[, seq(1, ncol(plot_mat), 2)] <- mat1
second <- mat2
second[is.na(second)] <- mat1[is.na(second)]
plot_mat[, seq(2, ncol(plot_mat), 2)] <- second
n_by_type <- table(plot_info$subtype)
labels <- paste0(subtype_levels, " (n = ", as.integer(n_by_type[subtype_levels]), ")")
split <- factor(labels[match(as.character(plot_info$subtype), subtype_levels)], levels = labels)
ht <- Heatmap(plot_mat, name = "Major cell type", col = celltype_col,
              cluster_rows = FALSE, cluster_columns = FALSE, cluster_row_slices = FALSE,
              cluster_column_slices = FALSE, row_split = split,
              column_split = factor(rep(species_order, each = 2), levels = species_order),
              column_title = species_order, column_gap = unit(2, "mm"),
              show_row_names = FALSE, show_column_names = FALSE,
              row_title_rot = 0, row_title_gp = gpar(fontsize = 8), row_gap = unit(1.5, "mm"),
              rect_gp = gpar(col = NA), use_raster = TRUE, raster_quality = 3,
              raster_resize_mat = FALSE,
              heatmap_legend_param = list(title = "Major cell type"))
pdf(file.path(outdir, "13_all_orthogroups_celltype_grouped.pdf"),
    width = 11, height = max(10, length(subtype_levels) * 0.45 + 3))
draw(ht, heatmap_legend_side = "right")
dev.off()
# Compact overview counts, independent of heatmap raster resolution.
p <- counts %>% filter(status != "Insufficient data") %>%
  mutate(subtype = factor(subtype, levels = rev(subtype_levels))) %>%
  ggplot(aes(n_orthogroups, subtype, fill = status)) + geom_col() +
  geom_text(aes(label = n_orthogroups), hjust = -0.15, size = 3) +
  scale_fill_manual(values = c(Conserved = "#2A9D8F", Turnover = "#BE0032")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Number of orthogroups", y = NULL) + theme_classic() + theme(legend.position = "bottom")
ggsave(file.path(outdir, "14_celltype_subtype_counts.pdf"), p,
       width = 10, height = max(4, length(subtype_levels) * 0.3 + 1.5), limitsize = FALSE)
save_tsv(tibble(parameter = c("species", "minimum_conserved_species", "minimum_turnover_species", "low_expression_quantile", "top1_top2_cutoff", "Dre_mapping_source"),
                value = c(paste(species_order, collapse = ","), min_conserved_species, min_turnover_species,
                          sc_low_expr_quantile, sc_top1_top2_cutoff, dre_mapping_source)), "15_settings.tsv")
capture.output(sessionInfo(), file = file.path(outdir, "16_sessionInfo.txt"))
cat("\nTotal single-copy orthogroups:", length(all_genes),
    "\nComparable / plotted orthogroups:", nrow(plot_info),
    "\nConserved:", sum(celltype_gene_status$status == "Conserved"),
    "\nTurnover:", sum(celltype_gene_status$status == "Turnover"),
    "\nInsufficient:", sum(celltype_gene_status$status == "Insufficient data"),
    "\nOutput:", outdir, "\n")
# Tissue-level testis enrichment within the cell-type Conserved/Turnover groups.
# Read genome-wide TPM tables; do not restrict to the curated sex-related genes.
# Includes testis-only and dual enrichment. Dre has no tissue data here.
tissue_base <- "mid_data/RNA_counts/tissue_specific_analysis/02.analysis_241013/01.tis_spec_genes/raw_tau"
tissue_files <- list(
  Hab = file.path(tissue_base, c("250325/Ha_Hav.Fo.tau.tsv", "Ha_Ha3.Fo.tau.tsv",
                                 "Ha_Ha2.Fo.tau.tsv", "Ha_Ha1.Fo.tau.tsv")),
  Cin = file.path(tissue_base, "250325/papa.Fo.tau.tsv"),
  Mbr = file.path(tissue_base, "250325/HJ.Fo.tau.tsv")
)
tissue_cutoff <- 1
tissue_min_TPM <- 1
tissue_pseudocount <- 0.1
missing_tissue_files <- unlist(tissue_files)[!file.exists(unlist(tissue_files))]
if (length(missing_tissue_files)) stop("Missing tissue files:\n", paste(missing_tissue_files, collapse = "\n"))

read_tissue_tpm <- function(path, sp) {
  d <- read.table(path, check.names = FALSE)
  if (identical(rownames(d), as.character(seq_len(nrow(d)))))
    stop("No gene row names detected in tissue file: ", path)
  cols <- names(d)[str_detect(names(d), "_(Mt|Fo|Mb|Mk|Fb|Fk)$")]
  if (!length(cols)) stop("No expected tissue columns in: ", path)
  if (!all(vapply(d[cols], is.numeric, logical(1)))) stop("Non-numeric TPM columns: ", path)
  x <- as_tibble(d[cols], rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample", values_to = "TPM") %>%
    mutate(gene_id = clean_id(gene_id), tissue = str_extract(sample, "[^_]+$"),
           species = sp, source_file = path)
  # Match the Cin tissue-ID convention in the supplied tissue script.
  if (sp == "Cin") x <- x %>% mutate(
    gene_id = if_else(str_starts(gene_id, "unm_"),
                      str_replace(gene_id, "^(.+)_\\1$", "\\1"),
                      str_extract(gene_id, "^[^_]+")))
  if (any(!is.finite(x$TPM) & !is.na(x$TPM)) || any(x$TPM < 0, na.rm = TRUE))
    stop("Invalid TPM values: ", path)
  if (anyDuplicated(x[c("gene_id", "tissue")]))
    stop("Duplicated normalized gene/tissue entries: ", path)
  x
}
tissue_raw <- imap_dfr(tissue_files, function(paths, sp) {
  map_dfr(paths, ~read_tissue_tpm(.x, sp))
})
# Hab: average each tissue across adult, 3m, 2m and 1m, as in the source.
# Cin/Mbr: one input file each. Somatic reference is max(brain, kidney).
tissue_mean <- tissue_raw %>% group_by(species, gene_id, tissue) %>%
  summarise(TPM = if (all(is.na(TPM))) NA_real_ else mean(TPM, na.rm = TRUE),
            n_values = sum(!is.na(TPM)), .groups = "drop")
save_tsv(tissue_mean, "17_tissue_TPM_and_stage_coverage.tsv")
tissue_grouped <- tissue_mean %>% mutate(tissue_group = case_when(
  tissue == "Mt" ~ "T", tissue == "Fo" ~ "O",
  tissue %in% c("Mb", "Mk") ~ "S_M", TRUE ~ "S_F")) %>%
  group_by(species, gene_id, tissue_group) %>%
  summarise(TPM = safe_max(TPM), .groups = "drop") %>%
  pivot_wider(names_from = tissue_group, values_from = TPM)
for (nm in c("T", "O", "S_M", "S_F")) {
  if (!nm %in% names(tissue_grouped)) tissue_grouped[[nm]] <- NA_real_
}
tissue_map <- native_map %>% filter(species %in% names(tissue_files), !is.na(gene_id)) %>% distinct()
if (anyDuplicated(tissue_map[c("species", "gene_id")]))
  stop("Non-unique native tissue gene to OG mapping")
tissue_og <- tissue_map %>% left_join(tissue_grouped, by = c("species", "gene_id")) %>%
  mutate(E_T = log2((T + tissue_pseudocount) / (S_M + tissue_pseudocount)),
         E_O = log2((O + tissue_pseudocount) / (S_F + tissue_pseudocount)),
         tissue_category = case_when(
           is.na(T) | is.na(O) | is.na(S_M) | is.na(S_F) ~ "missing_data",
           pmax(T, O) <= tissue_min_TPM ~ "low_expr",
           E_T > tissue_cutoff & T > tissue_min_TPM &
             E_O > tissue_cutoff & O > tissue_min_TPM ~ "gonad_biased",
           E_T > tissue_cutoff & T > tissue_min_TPM ~ "testis_biased",
           E_O > tissue_cutoff & O > tissue_min_TPM ~ "ovary_biased",
           TRUE ~ "other"))
tissue_celltype_detail <- crossing(
  Orthogroup = celltype_gene_status$Orthogroup,
  species = names(tissue_files)) %>%
  left_join(celltype_gene_status %>% select(Orthogroup, status, subtype), by = "Orthogroup") %>%
  left_join(tissue_og, by = c("Orthogroup", "species")) %>%
  mutate(tissue_category = replace_na(tissue_category, "missing_data"),
         tissue_informative = !tissue_category %in% c("missing_data", "low_expr"),
         testis_including_dual = tissue_category %in% c("testis_biased", "gonad_biased"))
save_tsv(tissue_celltype_detail, "18_celltype_tissue_gene_details.tsv")
tissue_percentages <- tissue_celltype_detail %>%
  filter(status %in% c("Conserved", "Turnover")) %>%
  group_by(status, species) %>% summarise(
    n_group_total = n(),
    n_tissue_informative = sum(tissue_informative),
    n_testis_only = sum(tissue_category == "testis_biased"),
    n_dual = sum(tissue_category == "gonad_biased"),
    n_testis_including_dual = sum(testis_including_dual),
    n_low_expr = sum(tissue_category == "low_expr"),
    n_missing = sum(tissue_category == "missing_data"), .groups = "drop") %>%
  complete(status = c("Conserved", "Turnover"), species = names(tissue_files),
           fill = list(n_group_total = 0L, n_tissue_informative = 0L,
                       n_testis_only = 0L, n_dual = 0L, n_testis_including_dual = 0L,
                       n_low_expr = 0L, n_missing = 0L)) %>%
  mutate(percent_testis_including_dual = if_else(n_tissue_informative > 0,
                                                 100 * n_testis_including_dual / n_tissue_informative, NA_real_))
save_tsv(tissue_percentages, "19_testis_bias_percent_by_celltype_status_species.tsv")
save_tsv(tissue_percentages %>% select(status, species, percent_testis_including_dual) %>%
           pivot_wider(names_from = species, values_from = percent_testis_including_dual),
         "20_testis_bias_percent_wide.tsv")
save_tsv(tibble(parameter = c("species", "enrichment_cutoff", "minimum_gonad_TPM", "pseudocount",
                              "numerator", "denominator", "Hab_stages", "somatic_reference"),
                value = c("Hab,Cin,Mbr", tissue_cutoff, tissue_min_TPM, tissue_pseudocount,
                          "testis_biased + gonad_biased", "Tissue-informative OGs within each cell-type status and species; excludes low_expr and missing_data",
                          "Mean per tissue across adult,3m,2m,1m; available values", "Sex-matched max(brain,kidney)")),
         "21_tissue_percentage_settings.tsv")
print(tissue_percentages, n = Inf, width = Inf)



# 各生殖细胞亚型、Sertoli及其明确两类型turnover的数量
celltype_gene_status %>%
  filter(
    subtype %in% c(
      "Conserved: spermatogonia",
      "Conserved: spermatocyte",
      "Conserved: spermatid",
      "Conserved: sertoli",
      "Turnover: sertoli <-> spermatogonia",
      "Turnover: sertoli <-> spermatocyte",
      "Turnover: sertoli <-> spermatid"
    )
  ) %>%
  count(subtype)

# 提取上述类别基因，并添加基因名称和各物种细胞类型
example_candidates <- celltype_gene_status %>%
  filter(
    subtype %in% c(
      "Conserved: spermatogonia",
      "Conserved: spermatocyte",
      "Conserved: spermatid",
      "Conserved: sertoli",
      "Turnover: sertoli <-> spermatogonia",
      "Turnover: sertoli <-> spermatocyte",
      "Turnover: sertoli <-> spermatid"
    )
  ) %>%
  select(Orthogroup, subtype, n_valid_species) %>%
  left_join(
    bridge %>% distinct(Orthogroup, symbol),
    by = "Orthogroup"
  ) %>%
  left_join(
    sc_states %>%
      select(
        Orthogroup, species,
        plot_celltype1, plot_celltype2
      ) %>%
      mutate(
        celltype = if_else(
          is.na(plot_celltype2),
          plot_celltype1,
          paste(plot_celltype1, plot_celltype2, sep = " + ")
        )
      ) %>%
      select(Orthogroup, species, celltype) %>%
      pivot_wider(
        names_from = species,
        values_from = celltype
      ),
    by = "Orthogroup"
  )

print(example_candidates, n = 50, width = Inf)

write_tsv(
  example_candidates,
  file.path(outdir, "22_examples_germ_subtypes_sertoli.tsv")
)
n_conserved <- sum(celltype_gene_status$status == "Conserved")
n_turnover  <- sum(celltype_gene_status$status == "Turnover")


test_result <- binom.test(
  x = n_turnover,
  n = n_conserved + n_turnover,
  p = 0.5,
  alternative = "two.sided"
)

print(test_result)

# not include Complex/ambiguous turnover
binom.test(
  x = 1980,
  n = 1706 + 1980,
  p = 0.5,
  alternative = "two.sided"
)
