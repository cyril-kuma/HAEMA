# Literature-guided redesign of mixed-host detection / denoising

**Status:** design plan (not yet implemented). **Date:** 2026-07-02.
**Motivation:** wet-lab known-ratio controls are not available in time (see
[`mixed_host_control_protocol.md`](mixed_host_control_protocol.md)); the denoising thresholds are
therefore calibrated and the mixed-detection architecture chosen **by adapting the validated
literature**, and validated **in silico** using the single-host controls we already have.

---

## 1. The core deduction (Logue et al. 2016)

Logue 2016 (the strongest HTS mixed-blood-meal paper in the corpus) detected and *quantified* mixed
meals **without any clustering, denoising-by-geometry, or assembly**. Their method:

1. Merge paired reads; strip primer/barcode; drop reads < 50 bp.
2. Collapse to **exact unique sequences**; **discard any unique observed < 10× across all samples**
   (their entire denoiser — true sequences recur, errors do not).
3. BLAST uniques vs `nt`, best match(es) at **> 90 % identity over full length**; ties → all species.
4. **Per mosquito, count read proportion per species.** Single host if one species > 90 % of reads;
   **mixed if a second species carries > 10 %** of reads.
5. Require **≥ 1000 reads/mosquito** (contamination / barcode-crosstalk guard).

**Principle: _classify-then-count_.** Mixed detection is an emergent property of counting
independently-classified reads per taxon at high depth (their avg 82,528 reads/mosquito made a 10 %
minor host = thousands of reads), not the output of a clustering algorithm. Logue explicitly warn
that read proportion ≠ blood proportion (amplification and mtDNA-copy bias) — the same caveat HÆMA
already encodes as `host_fractions_benchmarked: false`.

## 2. What the rest of the corpus says about detecting multiple hosts

| Paper | Platform | How multiple hosts are detected | What we adapt |
|-------|----------|--------------------------------|---------------|
| **Logue 2016** | Illumina 16S | classify-then-count; unique-seq + global ≥10 filter; >10 % minor | **the whole architecture** |
| **Reeves 2018** | Illumina COI | multiple ASV clusters per sample = multiple hosts; COI <200 bp → NUMT risk | ASV/OTU units; marker-length caution |
| **Kipp 2023** | **ONT** mitogenome | de-novo assembly yields multiple consensus contigs = multiple hosts; reference completeness limiting | ONT read behaviour; broad DB; consensus-per-taxon |
| **Channumsin 2021** | Sanger + cloning | chromatogram peak multiplicity; cloning confirmed **up to 4 hosts**; "only confident calls" | conservative calling; report uncertainty |
| **Santos 2019** | nested Sanger | double chromatogram peaks; **digestion (Sella) dominates ID success** | depth/quality gating by digestion |
| **Ogola 2017** | HRM, 3 markers | multiple melt peaks = mixed; multi-marker | multi-marker concordance |
| **Kent & Norris 2005** | multiplex PCR | multiple bands; detectable at **30:70** | a detection floor exists |
| **Townzen 2008 / Hadj-Henni 2015 / Altahir 2022 / Gueye 2023** | Sanger / PCR / ELISA | signal multiplicity; CytB multicopy; availability needed for preference | marker choice; framing |

**Unifying insight:** every method detects mixed meals as **multiplicity of independent
classification units per taxon** — reads (Logue), ASVs (Reeves), assembled contigs (Kipp),
chromatogram peaks (Sanger), melt peaks (HRM), gel bands (multiplex PCR). None depend on an
unsupervised geometric clustering of noisy reads. High depth + per-taxon counting is the HTS way.

## 3. Critical appraisal of the current HÆMA approach

Current path: `DENOISE_MIXED_TEMPLATES` (UMAP + HDBSCAN on k-mer spectra, or greedy identity
clustering) → `CLUSTER_CONSENSUS` (one consensus per cluster) → Medaka → BLAST **per cluster
consensus** → aggregate → RAMBO counts host reads per (sample, marker), calls mixed at ≥3 reads AND
≥1 % fraction.

**Where the fragility actually is:** the **pre-BLAST unsupervised clustering**.
- *Over-splitting:* ONT error can scatter one host's reads across several k-mer clusters → inflated
  cluster counts (already warned in `limitations.md` caveat 2) and unstable minor-host fractions.
- *Under-splitting / merging:* a rare minor host can be absorbed into a dominant cluster and never
  get its own consensus to BLAST → missed (this is the 50 % mixed-recovery gap seen in the earlier
  demonstration run).
- *Unbenchmarked knobs:* `mixed_denoise_min_cluster_size`, `mixed_denoise_min_cluster_fraction`
  have no principled value — they are the exact "not benchmarked" limitation.
- *k-mer geometry is host-blind:* clusters are formed by sequence composition, not by which host a
  read matches, so cluster boundaries need not correspond to host boundaries.

Crucially, **RAMBO already does Logue-style classify-then-count** — but on *cluster consensuses*, not
reads. The literature says: make the classification unit the **read (or a taxon-bin of reads)**, and
the fragile clustering step largely disappears.

## 4. Proposed redesign — reference-guided classify-then-count (ONT-adapted Logue)

Add a new, config-selectable mixed-detection strategy; **keep the existing one** (non-destructive,
mirrors how `--reference_mode` was added):

```
--mixed_detection_mode  {cluster_consensus | read_classification | hybrid}
                         default: cluster_consensus (unchanged) until validated
```

### Mode `read_classification` (recommended; direct Logue adaptation)

Per (sample, marker), after QC trimming:
1. **Classify every read** against the reference database with **minimap2** (ONT preset
   `map-ont`) — fast enough for full ONT depth (BLAST-per-read is too slow; minimap2 is the ONT
   idiom, and Kipp 2023 shows read-to-reference mapping is appropriate for ONT).
2. Assign each read to a host taxon by best hit passing identity/coverage floors; **ties → LCA**
   (reuse the existing conservative-LCA logic in `parse_blast_assignments.py`).
3. **Denoise by abundance, not geometry** — the Logue analogue for noisy ONT reads:
   - drop taxa supported by `< min_reads_per_host` (default 3, as RAMBO) **and**
     `< min_host_fraction` of marker reads (calibrated in §6);
   - optional **global filter**: drop a (taxon) seen in only one read across the whole run
     (Logue's ≥10 idea, scaled to ONT depth and per-run size).
4. **Count reads per taxon** → this is exactly the RAMBO input contract
   (`sample_uid, marker, host_assignment, count, fraction`), so **RAMBO, ecological indices,
   concordance and comparisons are unchanged downstream.**
5. Mixed if ≥ 2 host taxa pass the thresholds; single if one dominates; unresolved if none.

### Mode `hybrid` (best sequence quality + robust detection)

1. Classify reads (as above) and **bin reads by assigned taxon** (reference-guided binning replaces
   host-blind k-mer clustering — this is the key fix for over/under-splitting).
2. Build **one consensus per taxon-bin**, Medaka-polish it, and re-BLAST for a high-quality
   confirmatory species call (retains Kipp-style consensus quality and NUMT cross-checks).
3. Mixed detection still comes from **read counts per taxon-bin** (step 1), not from cluster counts.

This hybrid keeps everything good about the current pipeline (Medaka polishing, conservative LCA,
per-marker confirmation) while removing the unsupervised, unbenchmarked, host-blind clustering.

## 5. Nextflow implementation plan (staged, non-destructive)

**New/changed files (mirrors the `bin/ → module → subworkflow` convention in
`architecture_review.md`):**

- `bin/classify_reads.py` — minimap2 driver + per-read taxon assignment + abundance denoise +
  per-taxon read-count table (RAMBO-compatible schema). Reuses `parse_blast_assignments.py` LCA
  helpers; `--non_host_genera` vector filter applies here too.
- `modules/local/read_classification/main.nf` — process wrapping it (label `process_medium`,
  minimap2 container; emits the count table + a per-read assignment TSV for QC).
- `subworkflows/local/mixed_detection/main.nf` — `take` reads + reference; `if
  (params.mixed_detection_mode == 'read_classification') …` routes to the new module, else the
  existing `denoise → consensus → taxonomy` path; `emit` the same host-count channel so `main.nf`
  downstream is untouched.
- `nextflow.config` + `nextflow_schema.json` — add `mixed_detection_mode`,
  `min_host_fraction` (calibrated), `min_reads_per_host`, `read_class_min_identity`,
  `read_class_min_coverage`, `read_class_global_min_count`; add a `minimap2_container` (pinned).
- Docs: `methods.md` (both strategies + Logue provenance), `parameters.md`, `output.md`,
  and flip the relevant `limitations.md` caveat once validated.

**Feature-gated and reversible:** default stays `cluster_consensus`; switch to
`read_classification`/`hybrid` only after the in-silico validation (§6) supports it.

## 6. Validation **without** wet-lab — in-silico known-ratio mixtures

We cannot make wet-lab known-ratio controls in time, but we **already have single-host positive
controls** (each host sequenced alone). We can therefore build **synthetic known-ratio mixtures by
computationally pooling reads** from single-host control barcodes:

1. Take single-host control read sets (e.g. Human-only, Cattle-only, Chicken-only barcodes).
2. Sub-sample and combine at **defined read ratios** — 99:1, 95:5, 90:10, 75:25, 50:50 and
   reciprocals (a new helper `bin/make_insilico_mixtures.py`), producing synthetic barcodes with a
   **known minor-host read fraction**.
3. Run each mode (`cluster_consensus` vs `read_classification` vs `hybrid`) on these.
4. Metrics per mode × ratio × marker: **LOD** (minor host detected in ≥95 % of replicate draws),
   false-positive rate, and observed-vs-input minor fraction. Pick thresholds
   (`min_host_fraction`, `min_reads_per_host`, cluster knobs) from the ROC/recovery curves.

**Result (implemented — see [`denoising_calibration/README.md`](denoising_calibration/README.md)):**
The harness (`bin/make_insilico_mixtures.py`, `bin/classify_reads.py`,
`bin/benchmark_mixed_detection.py`) was run on 180 synthetic mixtures (Human/Cattle/Goat/Chicken,
99:1→50:50, ×3). Outcome: **LOD ≈ 5 % minor host**, and **`min_host_fraction = 0.02` is the operating
point (0 % false positives, 100 % detection ≥5 %)** — the current RAMBO 1 % default sits in the noise
floor (11–17 % false positives), and Logue's 10 % is over-conservative. Observed fractions track true
ratios monotonically (mild under-estimate) → qualitative only.

**What this validates and what it does not (critical-thinking honesty):**
- ✅ Calibrates the **bioinformatic** behaviour: detection floor, threshold values, and whether
  read_classification beats cluster_consensus for minor-host recovery — *exactly the denoising gap.*
- ❌ Does **not** capture extraction/PCR-amplification bias, mtDNA-copy differences, digestion, or
  the read-fraction↔blood-volume relationship. Those still require the wet-lab known-ratio controls
  in `mixed_host_control_protocol.md`. So this justifies threshold choice and architecture, but
  **`host_fractions_benchmarked` stays `false`** until wet-lab ratios exist.

This is the defensible "adapt from literature + validate in silico" path: the *architecture* is
taken from Logue/Reeves/Kipp, and the *thresholds* are chosen on synthetic mixtures built from our
own controls rather than guessed.

## 7. Scientific-framing safeguards (unchanged, re-affirmed)

- Read fractions remain **support evidence only**, never blood volume (Logue 2016).
- Vector self-hits excluded via `--non_host_genera`; controls excluded from indices.
- Conservative LCA on ties; species not overcalled from weak hits.
- Both strategies feed the **same** RAMBO / indices / concordance layers, so all existing
  denominators, confidence labels and uncertainty reporting are preserved.

## 8. Staged checklist

1. `bin/make_insilico_mixtures.py` + build synthetic mixtures from existing single-host controls.
2. `bin/classify_reads.py` (minimap2 + LCA + abundance denoise) — unit-tested against a tiny ref.
3. `read_classification` module + `mixed_detection` subworkflow (gated; default unchanged).
4. Run all three modes on synthetic mixtures → LOD/threshold curves → choose defaults.
5. If read_classification/hybrid wins, flip the default; update methods/limitations; add nf-test.
6. Re-run the subset (real controls) to confirm parity/improvement, then the full run.

**Estimated effort:** ~2–3 focused days; no wet-lab dependency for steps 1–5.
