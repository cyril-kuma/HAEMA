#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
opts <- list()
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[[i]])
  opts[[key]] <- args[[i + 1]]
}

flag <- function(name, default = FALSE) {
  value <- opts[[name]]
  if (is.null(value)) {
    return(default)
  }
  tolower(as.character(value)) %in% c("1", "true", "yes", "y", "on")
}

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(check.names = FALSE))
  }
  read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

pkg_available <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

fail_if_strict <- function(condition, message, strict) {
  if (condition && strict) {
    stop(message, call. = FALSE)
  }
}

as_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(x))
}

make_feature_tax_table <- function(master, asv_counts) {
  if (nrow(asv_counts) == 0 || !all(c("feature_id", "marker", "sequence") %in% names(asv_counts))) {
    return(data.frame(check.names = FALSE))
  }
  feature_meta <- unique(asv_counts[, c("feature_id", "marker", "sequence"), drop = FALSE])
  tax_cols <- intersect(
    c(
      "marker", "sequence", "host_assignment", "taxon_rank", "assignment_status",
      "confidence", "pident", "coverage", "evalue", "bitscore", "blast_source",
      "fallback_used", "primary_assignment_status", "assignment_method"
    ),
    names(master)
  )
  tax <- unique(master[, tax_cols, drop = FALSE])
  merged <- merge(feature_meta, tax, by = intersect(c("marker", "sequence"), names(tax)), all.x = TRUE)
  merged <- merged[!duplicated(merged$feature_id), , drop = FALSE]
  rownames(merged) <- merged$feature_id
  merged
}

fallback_decontam <- function(otu, sample_data, threshold) {
  if (nrow(otu) == 0 || nrow(sample_data) == 0 || !"sample_type" %in% names(sample_data)) {
    return(data.frame(
      feature_id = rownames(otu),
      contaminant = FALSE,
      method = "fallback_prevalence_no_negative_controls",
      p = NA_real_,
      neg_prevalence = 0,
      pos_prevalence = 0,
      max_negative_count = 0,
      stringsAsFactors = FALSE
    ))
  }
  sample_types <- tolower(as.character(sample_data[colnames(otu), "sample_type"]))
  neg <- grepl("negative", sample_types)
  pos <- !neg
  rows <- lapply(seq_len(nrow(otu)), function(i) {
    counts <- as.numeric(otu[i, ])
    neg_counts <- counts[neg]
    pos_counts <- counts[pos]
    neg_prev <- if (length(neg_counts)) mean(neg_counts > 0) else 0
    pos_prev <- if (length(pos_counts)) mean(pos_counts > 0) else 0
    max_neg <- if (length(neg_counts)) max(neg_counts, na.rm = TRUE) else 0
    contaminant <- isTRUE(neg_prev > 0 && (pos_prev <= neg_prev || max_neg >= max(pos_counts, 0, na.rm = TRUE)) && neg_prev >= threshold)
    data.frame(
      feature_id = rownames(otu)[i],
      contaminant = contaminant,
      method = "fallback_negative_control_prevalence",
      p = NA_real_,
      neg_prevalence = neg_prev,
      pos_prevalence = pos_prev,
      max_negative_count = max_neg,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

master <- read_tsv(opts[["master-endpoint"]])
asv_counts <- read_tsv(opts[["asv-count-table"]])
sample_summary <- read_tsv(opts[["sample-summary"]])
marker_summary <- read_tsv(opts[["marker-summary"]])
qc_summary <- read_tsv(opts[["qc-summary"]])
contamination_flags <- read_tsv(opts[["contamination-flags"]])
host_calls <- read_tsv(opts[["host-calls"]])

outdir <- opts[["output-dir"]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

enable_phyloseq <- flag("enable-phyloseq", TRUE)
enable_decontam <- flag("enable-decontam", TRUE)
strict_bioc <- flag("strict-bioconductor", FALSE)
contam_threshold <- as_numeric_safe(opts[["decontam-threshold"]])
if (is.na(contam_threshold)) {
  contam_threshold <- 0.5
}

metadata_cols <- intersect(
  c(
    "sample_uid", "run_id", "sample_id", "barcode_id", "specimen_id", "sample_type",
    "control_type", "expected_host_scientific_name", "expected_host_taxid", "expected_marker_result",
    "species", "sibling_species", "feeding_status", "collection_date", "collection_time",
    "collection_location", "bioclimatic_zone", "collection_region", "collection_cordinates",
    "latitude", "longitude", "collection_context", "collection_method", "specimen_sex",
    "extraction_batch", "pcr_batch", "library_batch", "barcode_kit", "flowcell", "basecalling_model"
  ),
  names(master)
)
sample_data <- unique(master[, metadata_cols, drop = FALSE])
if ("sample_uid" %in% names(sample_data)) {
  rownames(sample_data) <- sample_data$sample_uid
}

tax_table <- make_feature_tax_table(master, asv_counts)

sample_cols <- setdiff(names(asv_counts), c("feature_id", "marker", "sequence"))
otu <- as.matrix(data.frame())
if (nrow(asv_counts) > 0 && length(sample_cols) > 0) {
  otu <- as.matrix(asv_counts[, sample_cols, drop = FALSE])
  mode(otu) <- "numeric"
  rownames(otu) <- asv_counts$feature_id
  common_samples <- intersect(colnames(otu), rownames(sample_data))
  otu <- otu[, common_samples, drop = FALSE]
  sample_data <- sample_data[common_samples, , drop = FALSE]
}

negative_rows <- data.frame(check.names = FALSE)
if ("control_status" %in% names(master)) {
  negative_rows <- master[grepl("negative", master$control_status, ignore.case = TRUE), , drop = FALSE]
}

thresholds <- data.frame(
  marker = character(),
  host_assignment = character(),
  sequence = character(),
  max_negative_control_count = numeric()
)
if (nrow(negative_rows) > 0) {
  key_cols <- intersect(c("marker", "host_assignment", "sequence"), names(negative_rows))
  if (length(key_cols) > 0 && "count" %in% names(negative_rows)) {
    negative_rows$count <- as_numeric_safe(negative_rows$count)
    thresholds <- aggregate(
      negative_rows$count,
      by = negative_rows[, key_cols, drop = FALSE],
      FUN = max,
      na.rm = TRUE
    )
    names(thresholds)[names(thresholds) == "x"] <- "max_negative_control_count"
  }
}
write_tsv(thresholds, file.path(outdir, "qc_background_thresholds.tsv"))

phyloseq_status <- "disabled"
phyloseq_object <- NULL
if (enable_phyloseq) {
  if (pkg_available("phyloseq")) {
    ps_otu <- phyloseq::otu_table(otu, taxa_are_rows = TRUE)
    ps_tax <- phyloseq::tax_table(as.matrix(tax_table[rownames(otu), setdiff(names(tax_table), c("feature_id", "marker", "sequence")), drop = FALSE]))
    ps_sample <- phyloseq::sample_data(sample_data)
    phyloseq_object <- phyloseq::phyloseq(ps_otu, ps_tax, ps_sample)
    phyloseq_status <- "formal_phyloseq"
  } else {
    fail_if_strict(TRUE, "enable_phyloseq was requested but the R package 'phyloseq' is not installed.", strict_bioc)
    phyloseq_object <- list(
      otu_table = otu,
      tax_table = tax_table,
      sample_data = sample_data,
      note = "Fallback phyloseq-like object because package 'phyloseq' is unavailable."
    )
    class(phyloseq_object) <- c("haema_phyloseq_fallback", "list")
    phyloseq_status <- "fallback_phyloseq_like"
  }
}
saveRDS(phyloseq_object, file = file.path(outdir, "bloodmeal_phyloseq.rds"))

decontam_status <- "disabled"
decontam_results <- data.frame(
  feature_id = rownames(otu),
  contaminant = FALSE,
  method = "disabled",
  p = NA_real_,
  neg_prevalence = NA_real_,
  pos_prevalence = NA_real_,
  max_negative_count = NA_real_
)
if (enable_decontam) {
  if (pkg_available("decontam") && pkg_available("phyloseq") && inherits(phyloseq_object, "phyloseq")) {
    sample_df <- as.data.frame(phyloseq::sample_data(phyloseq_object))
    sample_df$is_negative_control <- grepl("negative", sample_df$sample_type, ignore.case = TRUE)
    phyloseq::sample_data(phyloseq_object)$is_negative_control <- sample_df$is_negative_control
    decontam_df <- decontam::isContaminant(
      phyloseq_object,
      method = "prevalence",
      neg = "is_negative_control",
      threshold = contam_threshold
    )
    decontam_results <- data.frame(
      feature_id = rownames(decontam_df),
      contaminant = decontam_df$contaminant,
      method = "decontam_prevalence",
      p = decontam_df$p,
      neg_prevalence = NA_real_,
      pos_prevalence = NA_real_,
      max_negative_count = NA_real_,
      stringsAsFactors = FALSE
    )
    decontam_status <- "formal_decontam_prevalence"
  } else {
    fail_if_strict(
      TRUE,
      "enable_decontam was requested but packages 'decontam' and/or 'phyloseq' are unavailable.",
      strict_bioc
    )
    decontam_results <- fallback_decontam(otu, sample_data, contam_threshold)
    decontam_status <- "fallback_negative_control_prevalence"
  }
}
write_tsv(decontam_results, file.path(outdir, "decontam_results.tsv"))

contaminant_feature_ids <- character()
if (nrow(decontam_results) > 0 && all(c("feature_id", "contaminant") %in% names(decontam_results))) {
  contaminant_feature_ids <- as.character(decontam_results$feature_id[decontam_results$contaminant %in% TRUE])
}

clean_asv_counts <- asv_counts
if (nrow(clean_asv_counts) > 0 && "feature_id" %in% names(clean_asv_counts)) {
  clean_asv_counts <- clean_asv_counts[!(clean_asv_counts$feature_id %in% contaminant_feature_ids), , drop = FALSE]
}
write_tsv(clean_asv_counts, file.path(outdir, "asv_count_table_decontaminated.tsv"))

feature_lookup <- data.frame(check.names = FALSE)
if (nrow(asv_counts) > 0 && all(c("feature_id", "marker", "sequence") %in% names(asv_counts))) {
  feature_lookup <- unique(asv_counts[, c("feature_id", "marker", "sequence"), drop = FALSE])
}
master_with_feature <- master
contaminant_asv_ids <- character()
if (nrow(master) > 0 && nrow(feature_lookup) > 0 && all(c("marker", "sequence") %in% names(master))) {
  master_with_feature <- merge(master, feature_lookup, by = c("marker", "sequence"), all.x = TRUE, sort = FALSE)
  if (all(c("feature_id", "asv_id") %in% names(master_with_feature))) {
    contaminant_asv_ids <- unique(as.character(master_with_feature$asv_id[master_with_feature$feature_id %in% contaminant_feature_ids]))
  }
}

clean_host_calls <- host_calls
if (nrow(clean_host_calls) > 0 && "best_feature_id" %in% names(clean_host_calls)) {
  clean_host_calls <- clean_host_calls[!(clean_host_calls$best_feature_id %in% contaminant_asv_ids), , drop = FALSE]
}
write_tsv(clean_host_calls, file.path(outdir, "host_calls_decontaminated.tsv"))

contam_summary <- data.frame(
  metric = c(
    "phyloseq_status",
    "decontam_status",
    "strict_bioconductor",
    "decontam_threshold",
    "n_features",
    "n_contaminants"
  ),
  value = c(
    phyloseq_status,
    decontam_status,
    strict_bioc,
    contam_threshold,
    nrow(decontam_results),
    sum(decontam_results$contaminant %in% TRUE, na.rm = TRUE)
  )
)
write_tsv(contam_summary, file.path(outdir, "contamination_model_summary.tsv"))

ecology_object <- list(
  otu_table = asv_counts,
  tax_table = tax_table,
  sample_data = sample_data,
  phyloseq = phyloseq_object,
  decontam = decontam_results,
  master_endpoint = master,
  qc = list(
    sample_summary = sample_summary,
    marker_summary = marker_summary,
    qc_summary = qc_summary,
    contamination_flags = contamination_flags,
    background_thresholds = thresholds,
    contamination_model_summary = contam_summary
  ),
  provenance = list(
    object_type = "haema_bloodmeal_ecology",
    phyloseq_status = phyloseq_status,
    decontam_status = decontam_status,
    note = "Formal phyloseq/decontam objects are generated when packages are installed; otherwise fallback objects/results are emitted unless strict mode is enabled."
  )
)
class(ecology_object) <- c("haema_bloodmeal_ecology", "list")
saveRDS(ecology_object, file = file.path(outdir, "bloodmeal_ecology_data.rds"))

clean_sample_cols <- setdiff(names(clean_asv_counts), c("feature_id", "marker", "sequence"))
clean_otu <- as.matrix(data.frame())
clean_sample_data <- sample_data
clean_tax_table <- data.frame(check.names = FALSE)
if (nrow(clean_asv_counts) > 0 && length(clean_sample_cols) > 0) {
  clean_otu <- as.matrix(clean_asv_counts[, clean_sample_cols, drop = FALSE])
  mode(clean_otu) <- "numeric"
  rownames(clean_otu) <- clean_asv_counts$feature_id
  clean_common_samples <- intersect(colnames(clean_otu), rownames(sample_data))
  clean_otu <- clean_otu[, clean_common_samples, drop = FALSE]
  clean_sample_data <- sample_data[clean_common_samples, , drop = FALSE]
  if (nrow(tax_table) > 0) {
    clean_tax_table <- tax_table[intersect(rownames(clean_otu), rownames(tax_table)), , drop = FALSE]
  }
}

phyloseq_object_decontaminated <- NULL
if (enable_phyloseq) {
  if (pkg_available("phyloseq")) {
    ps_otu_clean <- phyloseq::otu_table(clean_otu, taxa_are_rows = TRUE)
    ps_tax_clean <- phyloseq::tax_table(as.matrix(clean_tax_table[rownames(clean_otu), setdiff(names(clean_tax_table), c("feature_id", "marker", "sequence")), drop = FALSE]))
    ps_sample_clean <- phyloseq::sample_data(clean_sample_data)
    phyloseq_object_decontaminated <- phyloseq::phyloseq(ps_otu_clean, ps_tax_clean, ps_sample_clean)
  } else {
    phyloseq_object_decontaminated <- list(
      otu_table = clean_otu,
      tax_table = clean_tax_table,
      sample_data = clean_sample_data,
      note = "Fallback decontaminated phyloseq-like object because package 'phyloseq' is unavailable."
    )
    class(phyloseq_object_decontaminated) <- c("haema_phyloseq_fallback", "list")
  }
}

ecology_object_decontaminated <- list(
  otu_table = clean_asv_counts,
  host_calls = clean_host_calls,
  tax_table = clean_tax_table,
  sample_data = clean_sample_data,
  phyloseq = phyloseq_object_decontaminated,
  decontam = decontam_results,
  contaminants_removed = contaminant_feature_ids,
  raw_master_endpoint = master,
  provenance = list(
    object_type = "haema_bloodmeal_ecology_decontaminated",
    phyloseq_status = phyloseq_status,
    decontam_status = decontam_status,
    note = "Rows flagged as contaminants are removed only from this derived object and decontaminated endpoint tables; raw endpoint files are preserved."
  )
)
class(ecology_object_decontaminated) <- c("haema_bloodmeal_ecology_decontaminated", "list")
saveRDS(ecology_object_decontaminated, file = file.path(outdir, "bloodmeal_ecology_data_decontaminated.rds"))

manifest <- data.frame(
  output = c(
    "bloodmeal_ecology_data.rds",
    "bloodmeal_ecology_data_decontaminated.rds",
    "bloodmeal_phyloseq.rds",
    "decontam_results.tsv",
    "host_calls_decontaminated.tsv",
    "asv_count_table_decontaminated.tsv",
    "qc_background_thresholds.tsv",
    "contamination_model_summary.tsv"
  ),
  description = c(
    "RDS object containing otu_table, tax_table, sample_data, phyloseq/decontam slots, endpoint, and QC slots",
    "Derived RDS object after formal/fallback contaminant feature removal; raw endpoint rows are preserved elsewhere",
    paste("Formal phyloseq object or fallback phyloseq-like object:", phyloseq_status),
    paste("Formal decontam results or fallback negative-control prevalence model:", decontam_status),
    "Host-call table after removing calls supported only by contaminant-flagged features",
    "ASV/consensus count table after removing contaminant-flagged features",
    "Maximum observed feature/background counts in negative controls",
    "Status and counts for the R contamination/phyloseq stage"
  )
)
write_tsv(manifest, file.path(outdir, "r_outputs_manifest.tsv"))
