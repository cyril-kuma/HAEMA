#!/usr/bin/env python3
"""Compute host ecology statistical comparisons for HÆMA blood-meal data.

Takes the final host-call table and metadata as input and produces:
  - pairwise_hbi_comparisons.tsv: Fisher's exact tests for HBI between strata
  - pairwise_mixed_feeding_comparisons.tsv: Fisher's exact tests for mixed-feeding rates
  - host_use_statistical_tests_summary.tsv: Summary of all tests

Recommended tests:
  1. Fisher's exact test for pairwise HBI comparisons between zones
  2. Fisher's exact test for pairwise HBI comparisons between sibling species (if n permits)
  3. Fisher's exact test for mixed-feeding rate comparisons
  4. Holm correction within each test family
  5. Small-n warnings when any stratum has n < 5
  6. Optional Kruskal-Wallis test for host richness between zones (if meaningful)

All analyses are labeled as exploratory unless pre-specified as primary.
p-values are reported with Holm correction; raw p-values are also retained.
"""
import argparse
import csv
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path


def read_tsv(path):
    """Read a TSV file and return list of dicts."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


HUMAN = 'Homo sapiens'
UNRESOLVED_HOSTS = {'', 'unassigned', 'ambiguous', 'no_host_signal', 'unresolved'}
# Vector / non-host genera whose mtDNA can co-amplify against a broad reference database
# and must not be counted as blood-meal hosts (matched by genus / first token).
DEFAULT_NON_HOST_GENERA = {
    'Anopheles', 'Culex', 'Aedes', 'Ochlerotatus', 'Culiseta', 'Mansonia',
    'Culicidae', 'Culicinae', 'Anophelinae', 'Culicoides', 'Phlebotomus',
    'Lutzomyia', 'Simulium', 'Glossina', 'Ixodes', 'Rhipicephalus', 'Amblyomma',
}


def is_host(name, non_host_genera=DEFAULT_NON_HOST_GENERA):
    """True if a host_assignment names a real (vertebrate) host, not a vector/non-host."""
    if norm_taxon(name) in UNRESOLVED_HOSTS:
        return False
    genus = (name or '').replace('_', ' ').split()[0] if (name or '').strip() else ''
    return genus not in non_host_genera


def norm_taxon(name):
    """Normalise a taxon name for comparison: lowercase, strip, collapse whitespace."""
    return ' '.join((name or '').strip().lower().replace('_', ' ').split())


def parse_year_month(value):
    """Parse a YYYY-MM(-DD) collection date into (year, month); None if unparseable."""
    text = (value or '').strip()
    if len(text) < 7:
        return None
    try:
        year = int(text[0:4])
        month = int(text[5:7])
    except ValueError:
        return None
    if 1 <= month <= 12:
        return year, month
    return None


def build_host_sets(host_call_rows):
    """Aggregate the RAMBO host-call table to {sample_uid: set(host_taxa)}.

    Mirrors compute_ecological_indices.py: keeps only control_status == 'sample',
    drops unresolved/unassigned host calls, and collapses the (possibly multi-row,
    mixed-host, multi-marker) table to one host SET per specimen so a mosquito is
    counted once regardless of how many markers or hosts it carries. Also returns the
    set of all field-sample specimen ids (including tested-but-unidentified ones).
    """
    host_sets = defaultdict(set)
    all_uids = set()
    for r in host_call_rows:
        if (r.get('control_status') or '').strip() != 'sample':
            continue
        uid = r.get('sample_uid', '')
        if not uid:
            continue
        all_uids.add(uid)
        host = (r.get('host_assignment') or '').strip()
        if is_host(host):
            host_sets[uid].add(host)
    return host_sets, all_uids


def build_strata_maps(master_rows, zone_column, species_column, date_column, wet_months):
    """Map sample_uid -> stratum level for each stratification, from the master endpoint.

    Mirrors compute_ecological_indices.py exactly (sample_type == 'sample', one value
    per uid, season derived from the collection month). Returns an ordered dict of
    {stratum_type: {sample_uid: level}}.
    """
    zone_of, species_of, season_of = {}, {}, {}
    for r in master_rows:
        if r.get('sample_type') != 'sample':
            continue
        uid = r.get('sample_uid', '')
        if not uid or uid in zone_of:
            continue
        zone_of[uid] = (r.get(zone_column) or '').strip()
        species_of[uid] = (r.get(species_column) or '').strip()
        ym = parse_year_month(r.get(date_column, ''))
        if ym:
            _, month = ym
            season_of[uid] = 'wet' if month in wet_months else 'dry'
    return {
        'ecological_zone': zone_of,
        'sibling_species': species_of,
        'season': season_of,
    }


def counts_by_stratum(host_sets, all_uids, uid_to_level):
    """Count identified specimens, human-positive, mixed and host richness per stratum.

    Denominator is host-identified field specimens (non-empty host set); tested-but-
    unidentified specimens are excluded from these proportions, consistent with the
    Human Blood Index definition used in compute_ecological_indices.py.
    """
    agg = defaultdict(lambda: {'n': 0, 'n_human': 0, 'n_mixed': 0, 'hosts': set()})
    for uid in all_uids:
        level = uid_to_level.get(uid, '')
        if not level:
            continue
        hs = host_sets.get(uid)
        if not hs:
            continue
        d = agg[level]
        d['n'] += 1
        if HUMAN in hs:
            d['n_human'] += 1
        if len(hs) >= 2:
            d['n_mixed'] += 1
        d['hosts'] |= hs
    result = {}
    for level, d in agg.items():
        n = d['n']
        result[level] = {
            'n': n,
            'n_human': d['n_human'],
            'n_mixed': d['n_mixed'],
            'hbi': round(d['n_human'] / n, 6) if n else 0.0,
            'mixed_rate': round(d['n_mixed'] / n, 6) if n else 0.0,
            'richness': len(d['hosts']),
            'hosts': '|'.join(sorted(d['hosts'])),
            'small_n_warning': n < 5,
        }
    return result


def fisher_exact_2x2(a, b, c, d):
    """Compute two-tailed Fisher's exact test for 2x2 contingency table.
    
    | a  b  |
    | c  d  |
    
    Returns (odds_ratio, p_value).
    Uses scipy if available, otherwise a simplified implementation.
    """
    try:
        from scipy.stats import fisher_exact
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]])
        return float(odds_ratio), float(p_value)
    except ImportError:
        # Simplified implementation without scipy
        # For small tables, use the hypergeometric probability directly
        total = a + b + c + d
        if total < 2:
            return 1.0, 1.0
        # Use a simple approximation for very small tables
        # This is NOT a substitute for scipy for production use
        row1 = a + b
        col1 = a + c
        expected = (row1 * col1) / total if total > 0 else 0
        # Simple chi-square approximation (only for larger tables)
        if total >= 10 and expected > 0 and total - expected > 0:
            chi2 = ((a - expected) ** 2) / expected
            chi2 += ((b - (row1 - expected)) ** 2) / (row1 - expected) if row1 > expected else 0
            chi2 += ((c - (col1 - expected)) ** 2) / (col1 - expected) if col1 > expected else 0
            chi2 += ((d - (total - row1 - col1 + expected)) ** 2) / (total - row1 - col1 + expected) if (total - row1 - col1 + expected) > 0 else 0
            # Approximate p-value from chi-square with 1 df
            p_value = _chi2_survival(chi2, 1)
            return None, p_value
        return None, 1.0


def _chi2_survival(x, df):
    """Approximate chi-square survival function (1 - CDF)."""
    # Simple approximation for df=1
    if x <= 0:
        return 1.0
    # Use the relationship between chi-square(1) and normal distribution
    z = (x ** 0.5)
    # Approximation of 1 - Phi(z) using rational approximation
    p = 1.0 / (1.0 + 0.2316419 * z)
    p = p * (0.319381530 + p * (-0.356563782 + p * (1.781477937 + p * (-1.821255978 + p * 1.330274429))))
    p = p * (1.0 / (2.506628274654002)) * (2.718281828 ** (-z * z / 2))
    return p * 2  # Two-tailed


def holm_correction(p_values):
    """Apply Holm-Bonferroni correction to a list of p-values.
    
    Args:
        p_values: List of (index, p_value) tuples
        
    Returns:
        List of (index, raw_p, corrected_p) tuples
    """
    if not p_values:
        return []
    
    # Sort by p-value
    sorted_pvs = sorted(p_values, key=lambda x: x[1])
    
    corrected = []
    m = len(sorted_pvs)
    max_corrected = 1.0
    
    for i, (idx, raw_p) in enumerate(sorted_pvs):
        # Holm correction: multiply by (m - i), capped at 1.0
        corr_p = min(raw_p * (m - i), max_corrected)
        corrected.append((idx, raw_p, corr_p))
        max_corrected = corr_p
    
    # Restore original order
    corrected.sort(key=lambda x: x[0])
    return corrected


def pairwise_comparisons(stratum_data, value_key, label):
    """Compute pairwise Fisher's exact tests between all strata pairs.
    
    Args:
        stratum_data: dict from compute_hbi_by_stratum or compute_mixed_feeding_by_stratum
        value_key: 'hbi' -> use n_human; 'mixed_rate' -> use n_mixed
        label: Description for output
        
    Returns:
        List of comparison dicts
    """
    strata = list(stratum_data.keys())
    comparisons = []
    
    for s1, s2 in combinations(strata, 2):
        d1 = stratum_data[s1]
        d2 = stratum_data[s2]
        
        n1 = d1['n']
        x1 = d1[value_key]
        
        n2 = d2['n']
        x2 = d2[value_key]
        
        # 2x2 contingency table:
        # | x1    n1-x1  |  (event / no event in stratum 1)
        # | x2    n2-x2  |  (event / no event in stratum 2)
        a, b = x1, n1 - x1
        c, d = x2, n2 - x2
        
        odds_ratio, raw_p = fisher_exact_2x2(a, b, c, d)
        
        comparisons.append({
            'comparison': f'{s1} vs {s2}',
            'stratum_1': s1,
            'stratum_2': s2,
            'n_1': n1,
            'n_2': n2,
            'event_1': x1,
            'event_2': x2,
            'rate_1': round(x1 / n1, 6) if n1 > 0 else 0.0,
            'rate_2': round(x2 / n2, 6) if n2 > 0 else 0.0,
            'odds_ratio': round(odds_ratio, 4) if odds_ratio else None,
            'raw_p_value': round(raw_p, 6) if raw_p is not None else None,
            'corrected_p_value': None,  # Filled by Holm correction
            'significant_holm_005': None,  # Filled by Holm correction (corrected p < 0.05)
            'small_n_warning': (n1 < 5) or (n2 < 5),
            'test': 'fishers_exact',
            'test_family': label,
        })
    
    return comparisons


def kruskal_richness(host_sets, all_uids, uid_to_level):
    """Kruskal-Wallis test of per-specimen host richness across strata levels (C3).

    Per-specimen richness = number of distinct host taxa in that specimen's host set.
    Tests whether the distribution of per-specimen host richness differs across ecological
    zones. Host-unidentified specimens are excluded (consistent with the other indices).
    Exploratory; requires scipy (returns a null result with a note otherwise).
    """
    groups = defaultdict(list)
    for uid in all_uids:
        level = uid_to_level.get(uid, '')
        if not level:
            continue
        hs = host_sets.get(uid)
        if not hs:
            continue
        groups[level].append(len(hs))
    result = {
        'test': 'kruskal_wallis_host_richness',
        'n_groups': len(groups),
        'group_sizes': {k: len(v) for k, v in sorted(groups.items())},
        'statistic': None,
        'p_value': None,
        'small_n_warning': any(len(v) < 5 for v in groups.values()),
        'note': 'exploratory; per-specimen host richness across ecological zones',
    }
    if len(groups) < 2:
        result['note'] = 'need >=2 zones with host-identified specimens'
        return result
    try:
        from scipy.stats import kruskal
        stat, p = kruskal(*groups.values())
        result['statistic'] = round(float(stat), 4)
        result['p_value'] = round(float(p), 6)
    except ImportError:
        result['note'] = 'scipy not available; Kruskal-Wallis not computed'
    except ValueError as exc:
        result['note'] = f'Kruskal-Wallis not computed: {exc}'
    return result


def apply_holm(comparisons):
    """Fill corrected_p_value / significance flags in place via Holm-Bonferroni."""
    holm = holm_correction([(i, c['raw_p_value']) for i, c in enumerate(comparisons)
                            if c['raw_p_value'] is not None])
    for orig_idx, _raw_p, corr_p in holm:
        if orig_idx < len(comparisons):
            comparisons[orig_idx]['corrected_p_value'] = round(corr_p, 6)
            comparisons[orig_idx]['significant_holm_005'] = corr_p < 0.05
    return comparisons


COMPARISON_FIELDS = [
    'comparison', 'stratum_1', 'stratum_2', 'n_1', 'n_2', 'event_1', 'event_2',
    'rate_1', 'rate_2', 'odds_ratio', 'raw_p_value', 'corrected_p_value',
    'significant_holm_005', 'small_n_warning', 'test', 'test_family',
]


def write_comparisons(path, comparisons):
    """Write a comparison table, always emitting a header (even when empty)."""
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=COMPARISON_FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in comparisons:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description='Compute host ecology statistical comparisons for HÆMA blood-meal data'
    )
    parser.add_argument('--host-call-table', required=True,
                        help='Path to rambo_host_call_table.tsv')
    parser.add_argument('--master-endpoint', required=True,
                        help='Path to bloodmeal_master_endpoint.tsv (source of per-specimen strata)')
    parser.add_argument('--zone-column', default='collection_region',
                        help='master-endpoint column for ecological zone stratification')
    parser.add_argument('--species-column', default='sibling_species',
                        help='master-endpoint column for sibling species stratification')
    parser.add_argument('--date-column', default='collection_date',
                        help='master-endpoint column for collection date (used to derive season)')
    parser.add_argument('--wet-months', default='4,5,6,7,8,9,10',
                        help='Comma-separated month numbers classified as wet season')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for comparison tables')
    parser.add_argument('--exploratory-note', default='All statistical comparisons are exploratory unless pre-specified as primary in the study protocol.',
                        help='Note to include in output about exploratory nature')
    args = parser.parse_args()

    # Read input data
    host_call_rows = read_tsv(args.host_call_table)
    master_rows = read_tsv(args.master_endpoint)
    wet_months = {int(m) for m in str(args.wet_months).split(',') if m.strip().isdigit()}

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Aggregate host calls to one host-set per field specimen (controls excluded), and map
    # each specimen to its zone / sibling-species / season stratum from the master endpoint.
    host_sets, all_uids = build_host_sets(host_call_rows)
    strata = build_strata_maps(master_rows, args.zone_column, args.species_column,
                               args.date_column, wet_months)

    zone_stats = counts_by_stratum(host_sets, all_uids, strata['ecological_zone'])
    species_stats = counts_by_stratum(host_sets, all_uids, strata['sibling_species'])
    season_stats = counts_by_stratum(host_sets, all_uids, strata['season'])

    # === HBI comparisons: zones, sibling species, seasons ===
    hbi_zone_comparisons = apply_holm(pairwise_comparisons(zone_stats, 'n_human', 'hbi_zone'))
    hbi_species_comparisons = apply_holm(pairwise_comparisons(species_stats, 'n_human', 'hbi_species'))
    hbi_season_comparisons = apply_holm(pairwise_comparisons(season_stats, 'n_human', 'hbi_season'))
    write_comparisons(Path(args.output_dir) / 'pairwise_hbi_comparisons.tsv',
                      hbi_zone_comparisons + hbi_season_comparisons)
    write_comparisons(Path(args.output_dir) / 'pairwise_hbi_species_comparisons.tsv',
                      hbi_species_comparisons)

    # === Mixed-feeding comparisons between zones ===
    mixed_zone_comparisons = apply_holm(pairwise_comparisons(zone_stats, 'n_mixed', 'mixed_feeding_zone'))
    write_comparisons(Path(args.output_dir) / 'pairwise_mixed_feeding_comparisons.tsv',
                      mixed_zone_comparisons)

    # === Host richness by zone (descriptive; Kruskal-Wallis optional) ===
    richness_fields = ['stratum', 'n_specimens', 'host_richness_S', 'hosts', 'small_n_warning']
    with open(Path(args.output_dir) / 'host_richness_by_zone.tsv', 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=richness_fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for stratum, data in sorted(zone_stats.items()):
            writer.writerow({
                'stratum': stratum,
                'n_specimens': data['n'],
                'host_richness_S': data['richness'],
                'hosts': data['hosts'],
                'small_n_warning': data['small_n_warning'],
            })

    # === Kruskal-Wallis: per-specimen host richness across zones (C3, exploratory) ===
    kruskal_zone = kruskal_richness(host_sets, all_uids, strata['ecological_zone'])

    # === Write summary ===
    n_field = sum(1 for r in host_call_rows if (r.get('control_status') or '').strip() == 'sample')
    summary = {
        'host_call_rows': len(host_call_rows),
        'field_specimen_rows': n_field,
        'field_specimens_identified': sum(1 for u in all_uids if host_sets.get(u)),
        'hbi_zone_comparisons': len(hbi_zone_comparisons),
        'hbi_season_comparisons': len(hbi_season_comparisons),
        'hbi_species_comparisons': len(hbi_species_comparisons),
        'mixed_feeding_comparisons': len(mixed_zone_comparisons),
        'strata_zone': sorted(zone_stats.keys()),
        'strata_species': sorted(species_stats.keys()),
        'strata_season': sorted(season_stats.keys()),
        'kruskal_host_richness_zone': kruskal_zone,
        'exploratory_note': args.exploratory_note,
        'notes': [
            "All tests are Fisher's exact (two-tailed) on host-identified field specimens.",
            'A mosquito is counted once (host-set per sample_uid); controls are excluded.',
            'Holm-Bonferroni correction applied within each test family.',
            'Small-n warnings (n < 5 per stratum) indicate low power.',
            'scipy provides the exact Fisher test; without it a chi-square approximation is used (n >= 10 only) and odds_ratio is null.',
        ],
    }
    summary_path = Path(args.output_dir) / 'host_use_statistical_tests_summary.json'
    with open(summary_path, 'w') as fh:
        json.dump(summary, fh, indent=2)
        fh.write('\n')

    print('Host ecology comparisons completed:')
    print(f"  Field specimens identified: {summary['field_specimens_identified']}")
    print(f"  HBI zone comparisons: {len(hbi_zone_comparisons)}")
    print(f"  HBI species comparisons: {len(hbi_species_comparisons)}")
    print(f"  Mixed-feeding comparisons: {len(mixed_zone_comparisons)}")
    print(f"  Output dir: {args.output_dir}")


if __name__ == '__main__':
    main()
