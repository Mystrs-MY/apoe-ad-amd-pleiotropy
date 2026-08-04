library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
extension_dir <- file.path(
  root, "A1_protein_upgrade", "extensions", "PWAS2026_crosswalk_extension", "tables"
)
submission_dir <- file.path(root, "tables_submission", "supplementary_tables")
submission_root <- file.path(root, "tables_submission")
global_table_dir <- file.path(root, "tables")
upgrade_table_dir <- file.path(
  root, "A1_protein_upgrade", "tables", "supplementary_name_match_revision"
)
dir.create(submission_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(global_table_dir, recursive = TRUE, showWarnings = FALSE)

clean_text <- function(x) {
  if (!is.character(x)) return(x)
  x <- gsub("PWAS5", "five_protein_crosswalk", x, fixed = TRUE)
  x <- gsub("frozen_primary", "prespecified_primary", x, fixed = TRUE)
  x <- gsub("prespecified 25-protein", "prespecified 25-protein", x, fixed = TRUE)
  x <- gsub("frozen", "prespecified", x, fixed = TRUE)
  x <- gsub("not_reported_in_prespecified_crosswalk", "not_reported", x, fixed = TRUE)
  x
}

clean_table <- function(dt) {
  for (column in names(dt)) set(dt, j = column, value = clean_text(dt[[column]]))
  setnames(
    dt,
    old = intersect(
      c("frozen_reported_total", "recomputed_from_frozen_protein_rows"),
      names(dt)
    ),
    new = c(
      frozen_reported_total = "prespecified_reported_total",
      recomputed_from_frozen_protein_rows = "recomputed_from_prespecified_protein_rows"
    )[intersect(
      c("frozen_reported_total", "recomputed_from_frozen_protein_rows"),
      names(dt)
    )]
  )
  dt
}

# Table S29: submission-facing panel-to-evidence linkage without local paths or QA metadata.
s29 <- data.table(
  panel = c("Fig. 5b", "Fig. 5c", "Fig. 5d", "Fig. 5d"),
  displayed_content = c(
    "APOE variant-to-protein alpha versus protein-to-outcome beta linkage map",
    "Protein-level indirect effects and mediated proportions",
    "Literature-prioritized versus biology-guided panel boundary comparison",
    "Literature-prioritized versus biology-guided panel membership"
  ),
  primary_evidence_tables = c(
    "Tables S22, S23",
    "Table S25",
    "Tables S11, S27",
    "Tables S19, S27"
  ),
  source_data_object = c(
    "Figure_5b_APOE_linkage_scatter.csv",
    "Figure_5c_protein_mediation_forest.csv",
    "Figure_5d_panel_boundary_comparison.csv",
    "Figure_5d_panel_membership.csv"
  ),
  interpretation_role = c(
    "Alpha-beta triangulation; not de novo proteome-wide discovery",
    "Two-step mediation estimates for alpha- and beta-eligible proteins",
    "Sensitivity comparison across panel definitions and covariance assumptions",
    "Descriptive overlap and membership comparison"
  )
)

s29_name <- "TableS29_Figure5_Source_Data_and_QA_Manifest.tsv"
for (destination in c(submission_dir, submission_root, global_table_dir, upgrade_table_dir)) {
  fwrite(s29, file.path(destination, s29_name), sep = "\t", na = "NA")
}

# Table S30: scientific provenance and mapping fields only.
crosswalk <- fread(file.path(extension_dir, "PWAS5_crosswalk_mapping.tsv"), na.strings = c("", "NA"))
extraction <- fread(file.path(extension_dir, "PWAS5_extraction_QA.tsv"), na.strings = c("", "NA"))
integrity <- fread(file.path(extension_dir, "PWAS5_integrity_QA.tsv"), na.strings = c("", "NA"))
integrity <- integrity[, .(gene_symbol, data_integrity_pass = overall_pass)]
s30 <- merge(crosswalk, extraction, by = "gene_symbol", all.x = TRUE)
s30 <- merge(s30, integrity, by = "gene_symbol", all.x = TRUE)
s30[, UKB_PPP_OID := assay_target_ID]
s30[, prespecified_before_analysis := frozen_before_results]
s30 <- s30[, .(
  analysis_order,
  gene_symbol,
  literature_protein_name,
  literature_platform,
  literature_assay_id,
  literature_source,
  literature_source_sheet,
  literature_source_row,
  literature_P_value,
  assay_target_name,
  UKB_PPP_OID,
  UniProt_ID,
  Olink_panel,
  synapse_id,
  mapping_confidence,
  strict_mapping_eligible,
  selection_rule,
  prespecified_before_analysis,
  candidate_assay_count,
  protein_form,
  cross_platform_correlation,
  exact_assay_replication,
  interpretation_boundary,
  rows_scanned_autosomes,
  rs429358_match_count,
  rs7412_match_count,
  n_pqtl_candidates_after_APOE_exclusion,
  candidate_QC_rule,
  alpha_unique_pass,
  data_integrity_pass
)]
s30 <- clean_table(s30)[order(analysis_order)]
for (destination in c(submission_dir, submission_root, global_table_dir)) {
  fwrite(
    s30,
    file.path(destination, "TableS30_Five_Protein_PWAS_Data_Sources_and_Crosswalk.tsv"),
    sep = "\t",
    na = "NA"
  )
}

combine_sections <- function(files, labels, drop_columns = character()) {
  stopifnot(length(files) == length(labels))
  parts <- Map(function(file, label) {
    path <- file.path(extension_dir, file)
    if (!file.exists(path)) stop(sprintf("Missing extension table: %s", path))
    dt <- fread(path, na.strings = c("", "NA"))
    dt[, source_section := label]
    setcolorder(dt, c("source_section", setdiff(names(dt), "source_section")))
    dt
  }, files, labels)
  result <- rbindlist(parts, fill = TRUE, use.names = TRUE)
  result[, (intersect(drop_columns, names(result))) := NULL]
  clean_table(result)
}

s31 <- combine_sections(
  c(
    "PWAS5_APOE_alpha.tsv",
    "PWAS5_beta_main.tsv",
    "PWAS5_beta_cis.tsv",
    "PWAS5_harmonized_instruments_main.tsv",
    "PWAS5_harmonized_instruments_cis.tsv"
  ),
  c(
    "APOE_alpha",
    "protein_outcome_beta_main",
    "protein_outcome_beta_cis",
    "harmonized_instruments_main",
    "harmonized_instruments_cis"
  ),
  drop_columns = c("source_tar", "archive_path", "file_path")
)
for (destination in c(submission_dir, submission_root, global_table_dir)) {
  fwrite(
    s31,
    file.path(destination, "TableS31_Five_Protein_PWAS_APOE_Alpha_and_Beta.tsv"),
    sep = "\t",
    na = "NA"
  )
}

s32 <- combine_sections(
  c(
    "PWAS5_two_step_mediation_main.tsv",
    "PWAS5_two_step_mediation_cis.tsv",
    "PWAS5_two_step_mediation_strict.tsv",
    "PWAS5_incremental_aggregate_mediation.tsv",
    "PWAS5_covariance_mapping_sensitivity.tsv",
    "PWAS5_primary_reproduction_QA.tsv"
  ),
  c(
    "mediation_main",
    "mediation_cis",
    "mediation_strict",
    "incremental_aggregate_mediation",
    "covariance_mapping_sensitivity",
    "primary_reproduction_check"
  ),
  drop_columns = c("source_tar", "archive_path", "file_path")
)
for (destination in c(submission_dir, submission_root, global_table_dir)) {
  fwrite(
    s32,
    file.path(destination, "TableS32_Five_Protein_PWAS_Mediation_and_Aggregate_Sensitivity.tsv"),
    sep = "\t",
    na = "NA"
  )
}

redact_local_paths <- function(dt) {
  for (column in names(dt)) {
    values <- dt[[column]]
    if (!is.character(values)) next
    local_path <- !is.na(values) & grepl("^[A-Za-z]:[\\\\/]", values)
    resource_placeholder <- !is.na(values) & grepl(
      "^<(external_resource_root|project_root)>",
      values,
      ignore.case = TRUE
    )
    redact <- local_path | resource_placeholder
    normalized_paths <- gsub("\\\\", "/", values[redact])
    values[redact] <- sub("^.*/", "", normalized_paths)
    set(dt, j = column, value = values)
  }
  dt
}

is_p_value_column <- function(column) {
  normalized <- tolower(gsub("[. ]", "_", column))
  normalized %in% c("p", "p_value", "pvalue") ||
    grepl("(^p_|_p$|_p_value$|_pvalue$)", normalized)
}

format_underflow_p_values <- function(dt) {
  zero_tokens <- c("0", "0.0", "0.00", "0e+00", "0.00e+00")
  for (column in names(dt)) {
    if (!is_p_value_column(column)) next
    values <- as.character(dt[[column]])
    underflow <- !is.na(values) & trimws(values) %in% zero_tokens
    values[underflow] <- "<1e-300"
    set(dt, j = column, value = values)
  }
  dt
}

# Preserve source file names while removing workstation-specific absolute paths.
path_redaction_tables <- c(
  "TableS19_Protein_Provenance_Master.tsv",
  "TableS22_APOE_Variant_to_Protein_Alpha.tsv",
  "TableS25_Expanded_Primary_Two_Step_Mediation.tsv",
  "TableS26_Cis_Only_Two_Step_Mediation.tsv"
)

for (file in path_redaction_tables) {
  source <- file.path(submission_dir, file)
  packaged_source <- file.path(upgrade_table_dir, file)
  if (file.exists(packaged_source)) {
    file.copy(packaged_source, source, overwrite = TRUE)
  } else if (!file.exists(source)) {
    stop(sprintf("Missing packaged and submission table: %s", file))
  }
  if (!file.exists(source)) stop(sprintf("Missing submission table: %s", source))
  submission_table <- redact_local_paths(fread(source, na.strings = c("", "NA")))
  fwrite(submission_table, source, sep = "\t", na = "NA")

  root_copy <- file.path(submission_root, file)
  if (file.exists(root_copy)) file.copy(source, root_copy, overwrite = TRUE)
}

# A numerical P value cannot equal zero. Preserve analytic source tables and
# format underflow only in the submission-facing copies.
submission_files <- list.files(
  submission_dir,
  pattern = "\\.(tsv|csv)$",
  full.names = TRUE,
  ignore.case = TRUE
)
for (source in submission_files) {
  separator <- if (grepl("\\.tsv$", source, ignore.case = TRUE)) "\t" else ","
  submission_table <- format_underflow_p_values(
    fread(source, sep = separator, na.strings = c("", "NA"))
  )
  fwrite(submission_table, source, sep = separator, na = "NA")

  root_copy <- file.path(submission_root, basename(source))
  if (file.exists(root_copy)) file.copy(source, root_copy, overwrite = TRUE)
}

cat("Submission-facing Tables S18-S39 prepared; exact-zero P values were formatted as <1e-300.\n")
