#!/usr/bin/env Rscript
# Phyloseq-native manuscript figures from the HÆMA bloodmeal_phyloseq.rds object.
#
# These complement the Python/matplotlib figures (which read the flat TSVs) by visualising the
# phyloseq object directly with phyloseq + ggplot2: host-community composition, alpha diversity,
# beta-diversity ordination, a host x sample heatmap, and a decontam prevalence diagnostic.
#
# Caveats baked into the figures: this is a community-ecology toolkit applied to a small blood-meal
# host set (~5 vertebrate hosts, many single-host mosquitoes), so alpha diversity is "host breadth"
# (low resolution) and ordination reflects majority-host identity rather than a diversity gradient.
# Read abundances are evidence summaries, not validated quantitative diet. Controls are excluded
# from the ecological views. Outputs: vector PDF + 300 dpi PNG (SVG is produced for the Python
# figures). If the .rds is not a formal phyloseq object (R fallback path) the script no-ops cleanly.

suppressMessages({
  ok_ps <- requireNamespace("phyloseq", quietly = TRUE)
  library(methods)
})
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = "") {
  i <- match(flag, args); if (!is.na(i) && i < length(args)) args[i + 1] else default
}
phyloseq_rds <- getarg("--phyloseq")
decontam_tsv <- getarg("--decontam-results")
outdir       <- getarg("--outdir", "phyloseq_figures")
formats      <- strsplit(getarg("--formats", "pdf,png"), ",")[[1]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

note_skip <- function(msg) {
  cat(sprintf("  [skip] %s\n", msg))
  writeLines(msg, file.path(outdir, "phyloseq_figures_SKIPPED.txt"))
}

if (!ok_ps) { note_skip("phyloseq not installed in this container (R fallback path); no phyloseq figures."); quit(save = "no", status = 0) }
suppressMessages({ library(phyloseq); library(ggplot2) })

ps_raw <- tryCatch(readRDS(phyloseq_rds), error = function(e) NULL)
if (is.null(ps_raw) || !inherits(ps_raw, "phyloseq")) {
  note_skip(sprintf("%s is not a formal phyloseq object; no phyloseq figures.", phyloseq_rds)); quit(save = "no", status = 0)
}

# ---- shared style: colour-blind-safe host palette (matches the Python figures) ----
HOST_COLORS <- c("Homo sapiens" = "#0072B2", "Ovis aries" = "#E69F00", "Bos taurus" = "#009E73",
                 "Capra hircus" = "#CC79A7", "Canis lupus familiaris" = "#D55E00",
                 "unassigned" = "#BBBBBB")
theme_set(theme_bw(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey92", colour = NA),
                plot.title = element_text(face = "bold", size = 12),
                plot.caption = element_text(size = 7, hjust = 0, colour = "grey30")))

manifest <- list()
save_fig <- function(p, name, w, h, caption, inputs) {
  files <- c()
  for (fmt in formats) {
    f <- file.path(outdir, paste0(name, ".", fmt))
    tryCatch({
      if (fmt == "pdf") {
        if (capabilities("cairo")) grDevices::cairo_pdf(f, width = w, height = h) else pdf(f, width = w, height = h)
      } else {
        ttype <- if (capabilities("cairo")) "cairo" else NULL
        png(f, width = w, height = h, units = "in", res = 300, type = ttype)
      }
      print(p); dev.off(); files <- c(files, basename(f))
    }, error = function(e) { try(dev.off(), silent = TRUE); cat(sprintf("  [warn] %s.%s: %s\n", name, fmt, conditionMessage(e))) })
  }
  manifest[[length(manifest) + 1]] <<- data.frame(
    figure = name, files = paste(files, collapse = ";"), inputs = paste(inputs, collapse = ";"),
    caption = gsub("\n", " ", caption), stringsAsFactors = FALSE)
  cat(sprintf("  [ok] %s: %s\n", name, paste(files, collapse = ", ")))
}
host_scale <- function(aes_fn = scale_fill_manual)
  aes_fn(values = HOST_COLORS, na.value = "#BBBBBB", name = "Host")

# ---- agglomerate ASVs to the host_assignment rank (374 noise ASVs -> ~6 host taxa) ----
psh <- tryCatch(tax_glom(ps_raw, taxrank = "host_assignment", NArm = FALSE), error = function(e) ps_raw)
taxa_names(psh) <- as.character(tax_table(psh)[, "host_assignment"])
# field samples only for ecology (controls excluded); drop empty taxa
field <- tryCatch(prune_samples(sample_data(psh)$sample_type == "sample", psh), error = function(e) psh)
field <- prune_taxa(taxa_sums(field) > 0, field)
# host-only view (drop the unassigned noise pseudo-taxon) for diversity / ordination
hosts <- tryCatch(prune_taxa(taxa_names(field) != "unassigned", field), error = function(e) field)
hosts <- prune_samples(sample_sums(hosts) > 0, hosts)
cat(sprintf("phyloseq: %d taxa x %d samples; field=%d samples; host taxa=%d\n",
            ntaxa(ps_raw), nsamples(ps_raw), nsamples(field), ntaxa(hosts)))

ABUND_CAP <- "Read abundances are evidence summaries, not validated quantitative diet."

# ---- Fig 11: host-community composition by ecological zone ----
tryCatch({
  rel <- transform_sample_counts(field, function(x) if (sum(x) > 0) x / sum(x) else x)
  p <- plot_bar(rel, fill = "host_assignment") +
    facet_grid(~ collection_region, scales = "free_x", space = "free") +
    host_scale() + labs(y = "Relative read abundance", x = "Field sample",
      title = "Blood-meal host composition (phyloseq) by ecological zone",
      caption = ABUND_CAP) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
          legend.position = "right") + geom_bar(stat = "identity", width = 0.9)
  save_fig(p, "figure_11_phyloseq_composition", 9.5, 5.2,
           paste("Figure 11. Host-community composition from the phyloseq object: relative read",
                 "abundance of each vertebrate host per field mosquito (ASVs agglomerated to",
                 "host_assignment; 'unassigned' = unresolved/noise), faceted by ecological zone.",
                 ABUND_CAP), c("bloodmeal_phyloseq.rds"))
}, error = function(e) note_skip(paste("figure 11:", conditionMessage(e))))

# ---- Fig 12: alpha diversity (host breadth) by zone ----
tryCatch({
  if (ntaxa(hosts) < 1 || nsamples(hosts) < 2) stop("too few host taxa/samples")
  p <- plot_richness(hosts, x = "collection_region", color = "collection_region",
                     measures = c("Observed", "Shannon")) +
    geom_boxplot(alpha = 0.25, outlier.shape = NA) +
    labs(x = "Ecological zone", title = "Alpha diversity of the blood-meal host community",
         caption = paste("'Diversity' here = host breadth across <=5 vertebrate hosts (low",
                         "resolution); Observed = number of host taxa fed upon. Controls excluded.")) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")
  save_fig(p, "figure_12_phyloseq_alpha_diversity", 8.5, 5.0,
           paste("Figure 12. Alpha diversity (Observed host richness and Shannon index) of the",
                 "blood-meal host community per field mosquito, by ecological zone. Diversity is",
                 "host breadth over at most five vertebrate hosts (low-resolution); points are",
                 "mosquitoes, boxes summarise each zone. Controls excluded."),
           c("bloodmeal_phyloseq.rds"))
}, error = function(e) note_skip(paste("figure 12:", conditionMessage(e))))

# ---- Fig 13: beta-diversity ordination (PCoA, Bray-Curtis) — guarded/caveated ----
tryCatch({
  if (ntaxa(hosts) < 2 || nsamples(hosts) < 4) stop("ordination not meaningful (<2 host taxa or <4 samples)")
  ord <- ordinate(hosts, method = "PCoA", distance = "bray")
  p <- plot_ordination(hosts, ord, color = "collection_region") +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_brewer(palette = "Dark2", name = "Ecological zone") +
    labs(title = "Blood-meal host community ordination (PCoA, Bray–Curtis)",
         caption = paste("With so few host taxa and many single-host mosquitoes, ordination mostly",
                         "separates samples by their majority host — it is not a diversity gradient.",
                         "Descriptive only; controls excluded.")) +
    theme(legend.position = "right")
  save_fig(p, "figure_13_phyloseq_ordination", 7.8, 5.6,
           paste("Figure 13. Principal-coordinates ordination (Bray–Curtis) of field mosquitoes by",
                 "blood-meal host community, coloured by ecological zone. Because the host set is",
                 "small and many meals are single-host, separation chiefly reflects majority-host",
                 "identity rather than a community gradient; shown descriptively. Controls excluded."),
           c("bloodmeal_phyloseq.rds"))
}, error = function(e) note_skip(paste("figure 13 (ordination):", conditionMessage(e))))

# ---- Fig 14: host x sample heatmap (robust beta view) ----
tryCatch({
  if (ntaxa(hosts) < 1) stop("no host taxa")
  p <- plot_heatmap(hosts, method = NULL, sample.label = "sample_id", taxa.label = "host_assignment",
                    sample.order = "collection_region", low = "#deebf7", high = "#08306b",
                    na.value = "grey95") +
    labs(title = "Host x sample read-abundance heatmap (phyloseq)",
         caption = paste("Samples ordered by ecological zone; colour = reads (white = absent).", ABUND_CAP)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
          axis.text.y = element_text(face = "italic"))
  save_fig(p, "figure_14_phyloseq_heatmap", 9.5, 4.6,
           paste("Figure 14. Host x sample read-abundance heatmap from the phyloseq object (field",
                 "mosquitoes, ordered by ecological zone). A robust low-diversity companion to the",
                 "ordination.", ABUND_CAP), c("bloodmeal_phyloseq.rds"))
}, error = function(e) note_skip(paste("figure 14 (heatmap):", conditionMessage(e))))

# ---- Fig 15: decontam prevalence diagnostic (reconstructed from the object) ----
tryCatch({
  if (nchar(decontam_tsv) == 0 || !file.exists(decontam_tsv)) stop("no decontam_results.tsv")
  dec <- read.delim(decontam_tsv, stringsAsFactors = FALSE)
  # per-ASV prevalence in negative controls vs field samples, from the (un-agglomerated) object
  otu <- as(otu_table(ps_raw), "matrix"); if (!taxa_are_rows(ps_raw)) otu <- t(otu)
  smp <- data.frame(sample_data(ps_raw)); st <- smp$sample_type
  neg <- colnames(otu) %in% rownames(smp)[st == "negative_control"]
  pos <- colnames(otu) %in% rownames(smp)[st == "sample"]
  prev <- data.frame(
    feature_id = rownames(otu),
    neg_prev = if (sum(neg) > 0) rowMeans(otu[, neg, drop = FALSE] > 0) else 0,
    pos_prev = if (sum(pos) > 0) rowMeans(otu[, pos, drop = FALSE] > 0) else 0,
    reads = rowSums(otu))
  prev <- merge(prev, dec[, c("feature_id", "contaminant")], by = "feature_id", all.x = TRUE)
  prev$contaminant <- ifelse(is.na(prev$contaminant) | prev$contaminant == "FALSE",
                             "retained", "flagged contaminant")
  n_flag <- sum(prev$contaminant == "flagged contaminant")
  p <- ggplot(prev, aes(pos_prev, neg_prev, colour = contaminant, size = reads)) +
    geom_jitter(width = 0.01, height = 0.01, alpha = 0.8) +
    scale_colour_manual(values = c("retained" = "#4C72B0", "flagged contaminant" = "#C44E52"), name = NULL) +
    scale_size_continuous(name = "total reads", trans = "log10") +
    labs(x = "Prevalence in field samples", y = "Prevalence in negative controls",
         title = "decontam prevalence diagnostic",
         caption = sprintf(paste("Each point is an ASV feature; decontam flagged %d as contaminant(s)",
                                 "(prevalence threshold). Reconstructed from the phyloseq object +",
                                 "decontam_results.tsv."), n_flag))
  save_fig(p, "figure_15_phyloseq_decontam", 7.6, 5.2,
           sprintf(paste("Figure 15. decontam prevalence diagnostic: per-ASV prevalence in field",
                         "samples versus negative controls, point size = total reads, colour =",
                         "decontam contaminant call (%d flagged). Reconstructed from the phyloseq",
                         "object and decontam_results.tsv."), n_flag),
           c("bloodmeal_phyloseq.rds", "decontam_results.tsv"))
}, error = function(e) note_skip(paste("figure 15 (decontam):", conditionMessage(e))))

# ---- captions + manifest ----
if (length(manifest) > 0) {
  man <- do.call(rbind, manifest)
  write.table(man, file.path(outdir, "phyloseq_figure_manifest.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE)
  cap <- file(file.path(outdir, "phyloseq_figure_captions.md"), "w")
  writeLines("# HÆMA phyloseq figure captions (draft)\n", cap)
  for (i in seq_len(nrow(man))) {
    writeLines(sprintf("**%s** — files: `%s`\n\n%s\n\n_Input files:_ %s\n\n---\n",
                       man$figure[i], man$files[i], man$caption[i], man$inputs[i]), cap)
  }
  close(cap)
  cat(sprintf("\nWrote %d phyloseq figures + captions + manifest to %s\n", nrow(man), outdir))
} else {
  note_skip("no phyloseq figures produced")
}
