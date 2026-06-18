# Vector–host ecological indices

HÆMA summarises the blood-meal host calls into the **standard entomological indices** used in
malaria-vector ecology. They are produced by the `ECOLOGICAL_INDICES` step (whenever the RAMBO
host-call model runs) and written to:

- `05_endpoint_files/ecological_indices.tsv` — tidy table (`stratum_type`, `stratum`, `metric`,
  `value`, `ci_low`, `ci_high`, `n`, `detail`)
- `06_reports/ecological_indices_summary.json` — nested machine-readable summary
- `07_figures/figure_09_ecological_indices.{pdf,svg,png}` — forest plot + feeding-type partition +
  diversity ([figures.md](figures.md))

All indices are computed **per mosquito** (a sample), with hosts **unioned across the three markers**
and **controls excluded**. Proportions carry **Wilson score 95% confidence intervals** (valid at the
small n typical of stratified blood-meal data, unlike the normal/Wald interval). Estimates are
reported **overall** and stratified by **ecological zone** (`--eco_index_zone_column`), **vector
sibling species** (`--eco_index_species_column`), **collection period** (year-month, from
`--eco_index_date_column`), and **season** (wet/dry, `--eco_index_wet_months`, default Apr–Oct).

## Indices computed

| Index | Definition | Notes / citation |
|---|---|---|
| **Human Blood Index (HBI)** | proportion of host-identified mosquitoes whose meal **contains human** blood | The cornerstone malaria metric; vectorial capacity ∝ HBI². A mixed human+animal meal **counts as human-positive** (it is a human bite). Garrett-Jones 1964, *Bull WHO* 30:241–261. |
| **Animal blood index / zoophily** | proportion containing ≥1 **non-human** host | Overlaps HBI for mixed meals, so HBI + zoophily can exceed 1. |
| **Feeding-type partition** | **human-only / mixed (human+animal) / animal-only** | These three **do** sum to 1 — the unambiguous composition. |
| **Mixed-feeding (multiple-blood-meal) rate** | proportion feeding on ≥2 distinct host taxa | Mixed feeding underpins this pipeline's purpose. |
| **Host-specific blood indices** | proportion of meals **containing** each host taxon (bovine, ovine, caprine, canine, …) | Do **not** sum to 1 (a mixed meal is counted under each host it contains). |
| **Host-community diversity** | richness *S*, Shannon *H′*, Gini–Simpson (1−*D*), Pielou evenness *J′* | Computed on **mosquito incidence** per host, not read fractions. Shannon 1948; Simpson 1949; Pielou 1966; Levins 1968 (niche breadth). |

## How indices are defined here (the choices that matter)

- **Unit = mosquito.** A host is "fed upon" if detected on **any** of the three markers
  (`host_call_table.tsv`, control rows excluded). Cluster/marker counts are *not* the unit.
- **Denominator = host-identified mosquitoes.** Both `n_tested` and `n_identified` (and the
  `identification_rate` with CI) are reported so the denominator is transparent. Using identified
  meals assumes identification success is **not differential by host** — stated, not assumed silently.
- **Mixed meals.** Counted as positive for **every** host they contain (so HBI and host-specific
  indices behave as in the precipitin-era literature). The human-only/mixed/animal-only partition is
  reported alongside for an unambiguous, summing-to-one view.
- **Diversity uses incidence** (number of mosquitoes feeding on each host), which is robust to the
  fact that read fractions are **evidence, not validated quantitative diet**.

## Caveats (read before interpreting)

1. **Detection ≠ quantitative diet.** Indices are bounded by curated-panel coverage and denoising
   sensitivity; minority hosts (<~5% of reads) can be missed, which **biases mixed-feeding and
   zoophily downward**. See [limitations.md](limitations.md).
2. **Small per-stratum n.** Zone/species strata have n≈5–14, so confidence intervals are wide and
   strata flagged `small_n`. These are **descriptive**, not formal hypothesis tests; no pairwise
   significance testing is performed (it would be subgroup-fishing at this n). Mosquitoes from the
   same site/campaign are not fully independent (shared local conditions).
3. **Ecological framing.** A meta-regression of the three major African vectors found HBI tracks
   **collection location more strongly than sibling species** (Orsborne et al. 2018, *Malar J* 17:479,
   location R²=0.29 vs species R²=0.11) — so read the zone stratification as primary and the
   species stratification as exploratory.
4. **Temporal strata are confounded.** Sampling was opportunistic across a few campaigns spanning
   several years, and **each campaign is also a different site and sequencing batch**. The
   collection-period and season strata (and Figure 10) therefore describe *variation across
   campaigns* — they are **not** a temporal trend or an isolated seasonal effect, and no trend test
   is performed. The wet/dry split is a coarse label (Ghana's rainfall is regionally bimodal).

## Demonstration run (Ghana, 26 host-identified field mosquitoes)

Overall **HBI = 0.81 (95% CI 0.62–0.91)** with a **15% mixed-feeding rate**. The zone gradient is the
clear signal: **Coastal Savannah HBI = 1.0** (all human-only), **Forest HBI = 0.56** (dogs, cattle),
and **Northern Savannah** shows heavy livestock use (zoophily 0.71, mixed-feeding 0.57; sheep, goat).
The per-species differences are within overlapping wide intervals (n=5–14) — consistent with the
location-over-species pattern above. Numbers are in `ecological_indices.tsv`.

## Regenerate

`ECOLOGICAL_INDICES` runs automatically with the host-call model. Standalone:

```bash
python3 pipeline/bin/compute_ecological_indices.py \
  --host-calls     results/<run>/05_endpoint_files/host_call_table.tsv \
  --master-endpoint results/<run>/05_endpoint_files/bloodmeal_master_endpoint.tsv \
  --output-tsv     results/<run>/05_endpoint_files/ecological_indices.tsv \
  --output-json    results/<run>/06_reports/ecological_indices_summary.json
```

(stdlib-only — no special container needed.) Then Figure 9 is rendered by `build_figures.py`.

## Key references
- Garrett-Jones C. The human blood index of malaria vectors… *Bull World Health Organ*. 1964;30:241–261.
- Orsborne J, et al. Using the human blood index to investigate host biting plasticity… *Malar J*. 2018;17:479.
- Levins R. *Evolution in Changing Environments*. Princeton; 1968 (niche breadth).
- Shannon CE (1948); Simpson EH (1949); Pielou EC (1966) — diversity indices.
