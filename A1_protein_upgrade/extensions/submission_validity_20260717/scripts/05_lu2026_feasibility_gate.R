#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
extension_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
upgrade_root <- normalizePath(file.path(extension_root, "..", ".."), winslash = "/", mustWork = TRUE)

source_file <- file.path(extension_root, "data_raw", "lu_2026", "extracted", "43587_2026_1123_MOESM9_ESM.xlsx")
inventory_file <- file.path(upgrade_root, "data_processed", "syn51365303_inventory.tsv")
table_dir <- file.path(extension_root, "tables")
config_dir <- file.path(extension_root, "config")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(source_file), file.exists(inventory_file))

support <- as.data.table(read_excel(source_file, sheet = "Fig_6b"))
support_columns <- setdiff(names(support), c("protein_id", "protein_label"))
support[, n_supporting_dataset_outcomes := rowSums(as.matrix(.SD) == 1, na.rm = TRUE), .SDcols = support_columns]
dataset_support <- data.table(
  GNPC = support[["GNPC_APOE4_protein_AD"]] == 1,
  BF2 = support[["BF2_APOE4_protein_AD"]] == 1 | support[["BF2_APOE4_protein_Aβ"]] == 1,
  ADNI = support[["ADNI_APOE4_protein_Aβ"]] == 1
)
for (column in names(dataset_support)) set(dataset_support, which(is.na(dataset_support[[column]])), column, FALSE)
support[, n_supporting_datasets := rowSums(dataset_support)]
support[, supporting_datasets := apply(dataset_support, 1L, function(x) paste(names(which(x)), collapse = ";"))]
support[, SomaScan_aptamer := sub("^.*__", "", protein_id)]
support[, Lu_candidate_status := fifelse(
  n_supporting_dataset_outcomes >= 2L & n_supporting_datasets >= 2L,
  "candidate_supported_in_at_least_two_dataset_outcomes_and_two_independent_datasets",
  "below_result_blinded_post_publication_support_gate"
)]
support[, source_candidate_definition := paste(
  "APOE4-related mediator in Lu et al. 2026 Fig. 6b,",
  "derived from the source study's observational upstream-downstream mediation framework"
)]
support[, `:=`(
  source_study_first_author = "Lu",
  publication_year = 2026L,
  journal = "Nature Aging",
  DOI = "10.1038/s43587-026-01123-0",
  PMCID = "PMC13190297",
  source_data_location = "Fig. 6b source data (43587_2026_1123_MOESM9_ESM.xlsx)"
)]

aptamer_map <- as.data.table(read_excel(source_file, sheet = "Fig_6a", skip = 1))[, 1:4]
setnames(aptamer_map, c("SomaScan_aptamer", "UniProt_ID", "gene_symbol", "Lu_plot_label"))
aptamer_map <- unique(aptamer_map, by = "SomaScan_aptamer")
support <- merge(support, aptamer_map, by = "SomaScan_aptamer", all.x = TRUE)

inventory <- fread(inventory_file, na.strings = c("NA", ""))
parsed <- tstrsplit(sub("_v[0-9]+_.*$", "", inventory$name), "_", fixed = TRUE)
inventory[, `:=`(
  inventory_gene_symbol = parsed[[1]],
  inventory_UniProt_ID = parsed[[2]],
  Olink_target_ID = parsed[[3]]
)]
inventory <- inventory[grepl("^[A-Za-z0-9]+_[A-Z0-9]+_OID[0-9]+_v", name)]
inventory_map <- unique(inventory[, .(
  inventory_gene_symbol,
  inventory_UniProt_ID,
  Olink_target_ID,
  synapse_id,
  Synapse_file_name = name,
  Synapse_relative_path = relative_path
)])

support <- merge(
  support,
  inventory_map,
  by.x = c("gene_symbol", "UniProt_ID"),
  by.y = c("inventory_gene_symbol", "inventory_UniProt_ID"),
  all.x = TRUE
)
support[, exact_SomaScan_UniProt_to_Olink_inventory_match := !is.na(synapse_id)]
support[, eligible_for_Lu2026_exact_assay_extension :=
          n_supporting_dataset_outcomes >= 2L & n_supporting_datasets >= 2L &
          exact_SomaScan_UniProt_to_Olink_inventory_match]
support[, mapping_decision := fifelse(
  n_supporting_dataset_outcomes < 2L | n_supporting_datasets < 2L,
  "excluded_below_two_dataset_outcome_or_two_independent_dataset_support_gate",
  fifelse(
    is.na(UniProt_ID) | UniProt_ID == "",
    "excluded_SomaScan_aptamer_UniProt_unresolved",
    fifelse(
      exact_SomaScan_UniProt_to_Olink_inventory_match,
      "eligible_exact_SomaScan_aptamer_UniProt_to_Olink_assay",
      "excluded_no_exact_UniProt_matched_Olink_assay_in_UKB_PPP_inventory"
    )
  )
)]
support[, evidence_boundary := paste(
  "Lu et al. 2026 is an observational plasma/CSF proteomic study.",
  "The source mediation classification is not genetic causal evidence,",
  "and this candidate status does not alter the primary 25-protein panel."
)]
setorder(support, -n_supporting_datasets, -n_supporting_dataset_outcomes, gene_symbol, SomaScan_aptamer)
fwrite(support, file.path(table_dir, "Lu2026_cross_dataset_candidate_feasibility.tsv"), sep = "\t", na = "NA")

candidate_universe <- support[n_supporting_dataset_outcomes >= 2L & n_supporting_datasets >= 2L]
exact_candidates <- candidate_universe[eligible_for_Lu2026_exact_assay_extension == TRUE]
gate_passed <- nrow(exact_candidates) >= 5L

gate <- data.table(
  criterion = c(
    "Lu_candidates_supported_in_at_least_two_dataset_outcomes_and_two_independent_datasets",
    "exact_SomaScan_aptamer_UniProt_to_UKB_PPP_Olink_candidates",
    "minimum_exact_assay_extension_gate",
    "primary_panel_modified"
  ),
  observed = c(nrow(candidate_universe), nrow(exact_candidates), as.character(gate_passed), "FALSE"),
  threshold = c(">=1 for source-universe documentation", ">=5", "TRUE", "FALSE"),
  passed = c(nrow(candidate_universe) >= 1L, nrow(exact_candidates) >= 5L, gate_passed, TRUE)
)
fwrite(gate, file.path(table_dir, "Lu2026_exact_assay_extension_gate.tsv"), sep = "\t")

selection <- exact_candidates[, .(
  synapse_id,
  selection_status = "selected_for_download",
  gene_symbol,
  UniProt_ID,
  Olink_target_ID,
  SomaScan_aptamer,
  Lu_plot_label,
  n_supporting_dataset_outcomes,
  n_supporting_datasets,
  supporting_datasets,
  selection_reason = "Lu_2026_Fig6b_support_in_at_least_two_dataset_outcomes_and_two_independent_datasets_plus_exact_aptamer_UniProt_to_Olink_inventory_match"
)]
fwrite(selection, file.path(config_dir, "lu2026_exact_assay_download_selection.tsv"), sep = "\t", na = "NA")

if (!gate_passed) {
  stop("Lu 2026 exact-assay gate failed; candidate MR extension must not run.")
}

message("Lu 2026 exact-assay feasibility gate passed.")
message("Cross-dataset candidates: ", nrow(candidate_universe))
message("Exact-assay candidates selected: ", nrow(exact_candidates))
message("Selected genes: ", paste(exact_candidates$gene_symbol, collapse = ", "))
