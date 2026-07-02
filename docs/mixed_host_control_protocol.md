# Wet-lab protocol: known-ratio mixed-host controls for HÆMA calibration

**Purpose.** Produce mixed-host blood-meal controls of *known composition and known ratio* so that
the HÆMA ONT metabarcoding pipeline can be calibrated and its limits quantified:

1. **Detection** — the limit of detection (LOD) for a minor host as a function of its input
   fraction, sequencing depth, and marker.
2. **Quantitation** — whether, and how well, observed **read fractions** track **true template
   ratios** (and therefore whether host fractions may ever be reported as more than qualitative
   support evidence). This is what unlocks the `host_fractions_benchmarked` flag in the manifest.
3. **Threshold calibration** — evidence-based values for `mixed_denoise_min_cluster_size`,
   `mixed_denoise_min_cluster_fraction`, `rambo_min_host_reads`, and `rambo_min_host_fraction`.

> **Framing (do not violate).** This protocol characterises the read-fraction↔template relationship
> under controlled conditions. It does **not** license interpreting field read fractions as ingested
> blood volumes (Logue et al. 2016), and it says nothing about host preference (no availability data).

---

## 1. The central design principle — define the mixing basis

"Proportion" must be defined and measured, not assumed. Three bases differ:

| Basis | What it is | Problem |
|-------|-----------|---------|
| Blood volume | µL of each host blood pipetted | Ignores DNA-per-µL differences (nucleated avian/reptilian RBCs ≫ mammalian) |
| gDNA mass | ng of extracted DNA mixed | Ignores mtDNA copies per ng (varies by tissue/genome) |
| **mtDNA copies** | template molecules the primers actually see | **The correct denominator** for read fraction |

**Rule:** prepare mixtures gravimetrically/volumetrically *and* **quantify each single-host stock by
marker-specific qPCR or ddPCR** so every mixture has a *measured* mtDNA-copy ratio, not just a
nominal pipetted ratio. Report all three where possible; calibrate against the mtDNA-copy ratio.

---

## 2. Hosts (tailored to Ghanaian *An. gambiae* s.l. blood meals)

Prioritise the realistic peridomestic hosts plus one nucleated-RBC host to expose the DNA-vs-blood
confound:

- *Homo sapiens*, *Bos taurus* (cattle), *Capra hircus* (goat), *Ovis aries* (sheep),
  *Sus scrofa* (pig), *Canis lupus familiaris* (dog).
- **≥1 bird:** *Gallus gallus* (chicken) — nucleated RBCs; essential for testing quantitation bias.

Source from an abattoir / veterinary clinic / consented human donor under the approvals in §10.
Keep single-host stocks as the reference set; they double as the curated-panel/positive controls.

---

## 3. Experimental design

Tiered, so you can start cheap and controlled and add realism:

- **Tier 1 — gDNA mixtures (start here).** Extract each host separately, quantify, mix extracted
  gDNA at defined mtDNA-copy ratios. Maximum control; isolates the *sequencing/bioinformatics*
  behaviour from extraction/digestion variance. Best for threshold calibration.
- **Tier 2 — whole-blood mixtures.** Mix host bloods, then extract as a pool. Adds extraction and
  RBC-nucleation effects.
- **Tier 3 — membrane-fed mosquitoes (most realistic).** Feed *An. gambiae* (colony) through a
  membrane on defined blood mixtures; harvest abdomens; extract. Captures real inhibitors, the
  vector background (relevant to the vector self-amplification caveat), and — with a time course —
  digestion. Use for final validation, not initial calibration.

**Ratio gradient (per host pair), by measured mtDNA-copy ratio:**
100:0, 99:1, 95:5, 90:10, 75:25, 50:50 — and the reciprocals (so each host is both major and minor).
The 99:1–90:10 rungs define the LOD; 50:50 checks symmetry/primer bias.

**Three-host mixtures:** at least equal (33:33:33) and a skewed set (80:15:5) using a
human+cattle+goat/chicken combination reflecting real multi-host meals.

**Replication:** ≥3 independent replicates per condition (biological where feasible, else technical);
distribute across ≥2 barcodes and ideally ≥2 flow cells to capture batch/index-hop effects.

**Mandatory controls on every run:**
- No-template (PCR blank) and extraction blank (reagent contamination baseline → decontam).
- Single-host-only barcodes adjacent to mixtures (detect index hopping / crosstalk).
- A digestion time course (Tier 3): 0, 12, 24, 36, 48 h post-feed (Sella-scale analogue).

---

## 4. Quantify inputs (make the ratio "known")

For each single-host stock, before mixing:

1. **Total gDNA** — Qubit dsDNA (fluorometric, not NanoDrop).
2. **Marker-specific mtDNA copies** — qPCR (SYBR/probe) or, preferably, **ddPCR** for absolute copies
   per µL, using each of the three HÆMA markers (CytB, COI_short, COI_long) so per-marker primer
   efficiency is measured. Build standard curves from cloned/gBlock amplicons.
3. Compute the **measured mtDNA-copy ratio** for every mixture from these values. This is the
   x-axis for all calibration.

---

## 5. Library prep — match the production pipeline exactly

The calibration only transfers if the chemistry matches field runs:

- Primers: the tri-marker HÆMA set (`assets/primers.csv`) at production concentrations/cycling.
- Kit/chemistry: **SQK-RBK114-96**, **R10.4.1**, **Dorado SUP** (`r1041_e82_400bps_sup_v4.3.0`).
- Unique barcode per control; where possible place controls on the **same flow cells as field
  samples**, or run a dedicated calibration flow cell with identical settings.
- Depth: target ≥ field median depth per barcode; also over-sequence a subset to enable
  **downsampling** (§7) to find the minimum depth for reliable minor-host recovery.

---

## 6. Extraction

Use the same kit/protocol as field specimens (consistency > absolute yield). For Tier 3, extract
whole abdomens. Record digestion class per specimen (feeds the proposed `digestion_class` column).

---

## 7. Bioinformatics analysis (feeds directly back into HÆMA)

Run all controls through HÆMA with the **same parameters as field data**, then:

1. **Samplesheet.** Use the standard schema plus (proposed) `expected_host_fractions`,
   `mixing_basis` (blood/gDNA/mtDNA_copies), and `input_quant_method` columns (see §9). Set
   `sample_type = positive_control` with `expected_host_scientific_name` as `A;B[;C]`.
2. **Detection metrics.** Per mixing ratio × marker × depth: sensitivity (minor host detected in
   ≥95 % of replicates → **LOD**), specificity (no false hosts), and the existing
   `positive_control_check.tsv` recovery fields.
3. **Quantitation.** Regress observed `host_fraction` on measured mtDNA-copy fraction; report slope,
   intercept, R², and a **Bland–Altman** bias plot, **per marker**. A monotonic, low-scatter
   relationship → derive a calibration curve; otherwise conclude fractions are qualitative only.
4. **Threshold sweep.** Re-run (or post-process) across a grid —
   `mixed_denoise_min_cluster_size` ∈ {5,10,20,30}, `mixed_denoise_min_cluster_fraction` ∈
   {0.01,0.02,0.05,0.10}, `rambo_min_host_fraction` ∈ {0.005,0.01,0.02,0.05} — and plot
   recovery/false-positive (ROC-style) curves to choose an operating point that maximises minor-host
   recovery at acceptable false-positive rate.
5. **Depth.** Rarefy the over-sequenced subset to find the minimum reads/barcode for the chosen LOD.
6. **Per-marker.** Report CytB vs COI_short vs COI_long separately — primer bias and NUMT risk differ.

**Decision rule for `host_fractions_benchmarked`:** set `true` only if (3) shows a reproducible,
monotonic read-fraction↔copy-ratio relationship across ≥2 runs with an explicit, documented
uncertainty; otherwise keep `false` and report fractions as support evidence only.

---

## 8. Pitfalls and artifact controls

- **Index hopping** — unique/dual barcodes; single-host and no-template barcodes flanking mixtures.
- **Pipetting error at extremes** — reach 99:1 by serial dilution, not one 1 µL pipetting step.
- **Primer/amplification bias** — measured per marker by ddPCR (§4); it is *expected*, and is
  precisely what the calibration curve corrects for.
- **Chimeras / NUMTs** — inspect COI_short (moderate NUMT risk); cross-check CytB↔COI concordance.
- **Contamination** — extraction/PCR blanks + decontam; the vector's own mtDNA will appear (broad
  DB) and is handled by `--non_host_genera`.
- **Stochastic dropout** — low-input minor hosts drop out probabilistically; hence ≥3 replicates and
  the ≥95 % detection criterion for LOD.

---

## 9. Pipeline integration (what I can add on request)

To make these controls first-class inputs:

- **Samplesheet schema:** add optional `expected_host_fractions`, `mixing_basis`,
  `input_quant_method` columns (validated by `validate_inputs.py`), so `positive_control_check.tsv`
  can score *quantitative* recovery, not just presence.
- **Benchmarking script:** `bin/benchmark_mixed_controls.py` to compute LOD, the read-fraction↔copy
  regression + Bland–Altman, and the threshold-sweep curves from a completed run.
- These are additive; say the word and I will implement them and wire a `benchmarking` report
  section.

---

## 10. Ethics, permits, biosafety (Ghana context — confirm locally)

- **Human blood:** institutional review board approval + informed consent (e.g. Noguchi Memorial
  Institute for Medical Research / KCCR IRB).
- **Animal blood:** institutional animal ethics / veterinary approval; abattoir sourcing where
  appropriate reduces live-animal use.
- **Mosquito colony / membrane feeding:** insectary biosafety approval; follow institutional
  arthropod-containment SOPs.
- **Data/sample governance:** material transfer and data-sharing agreements as required.

*Confirm all approvals with your institution before any sampling; the above is a checklist, not
regulatory advice.*

---

## 11. Minimal viable version (if resources are tight)

If a full factorial is infeasible, the smallest experiment that still yields a defensible answer:

- 2 host pairs (e.g. Human+Cattle, Human+Chicken — one mammal-mammal, one mammal-bird),
- 4 ratios (95:5, 90:10, 75:25, 50:50) + reciprocals,
- 3 replicates, Tier 1 (gDNA) mixtures, ddPCR-quantified,
- one flow cell, standard depth.

This gives a per-marker LOD and a first read-fraction↔copy-ratio curve — enough to justify (or
refuse) quantitative reporting and to set denoising thresholds.
