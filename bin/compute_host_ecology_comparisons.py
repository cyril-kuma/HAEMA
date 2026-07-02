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


def read_metadata(path):
    """Read samplesheet metadata and return dict keyed by specimen_id."""
    rows = read_tsv(path)
    meta = {}
    for row in rows:
        key = row.get('specimen_id') or row.get('sample_id') or row.get('barcode_id', '')
        if key:
            meta[key] = row
    return meta


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


def compute_hbi_by_stratum(host_call_rows, metadata, stratum_column):
    """Compute HBI (Human Blood Index) by stratum.
    
    HBI = proportion of host-identified mosquitoes with human-containing blood meal.
    Controls are excluded.
    
    Returns:
        dict: {stratum_value: {'n': total, 'n_human': human_count, 'hbi': proportion}}
    """
    hbi_by_stratum = defaultdict(lambda: {'n': 0, 'n_human': 0})
    
    for row in host_call_rows:
        specimen_id = row.get('specimen_id') or row.get('sample_id', '')
        control_status = row.get('control_status', '') or ''
        
        # Exclude controls
        if control_status in ('extraction_blank', 'pcr_negative', 'positive_control'):
            continue
        
        # Get stratum value from metadata
        meta = metadata.get(specimen_id, {})
        stratum_value = meta.get(stratum_column, 'unknown')
        
        # Check if host-identified
        host_call = row.get('host_call', '') or row.get('scientific_name', '') or ''
        if not host_call or host_call in ('', 'unassigned', 'no_host_signal'):
            continue
        
        hbi_by_stratum[stratum_value]['n'] += 1
        
        # Check for human in host call (may be mixed)
        if 'homo sapiens' in host_call.lower():
            hbi_by_stratum[stratum_value]['n_human'] += 1
    
    # Compute HBI proportions
    result = {}
    for stratum, counts in hbi_by_stratum.items():
        n = counts['n']
        n_human = counts['n_human']
        hbi = n_human / n if n > 0 else 0.0
        result[stratum] = {
            'n': n,
            'n_human': n_human,
            'hbi': round(hbi, 6),
            'small_n_warning': n < 5,
        }
    
    return result


def compute_mixed_feeding_by_stratum(host_call_rows, metadata, stratum_column):
    """Compute mixed-feeding rate by stratum.
    
    Mixed feeding = proportion of host-identified mosquitoes with >= 2 distinct host taxa.
    Controls are excluded.
    
    Returns:
        dict: {stratum_value: {'n': total, 'n_mixed': mixed_count, 'rate': proportion}}
    """
    mixed_by_stratum = defaultdict(lambda: {'n': 0, 'n_mixed': 0})
    
    for row in host_call_rows:
        specimen_id = row.get('specimen_id') or row.get('sample_id', '')
        control_status = row.get('control_status', '') or ''
        
        if control_status in ('extraction_blank', 'pcr_negative', 'positive_control'):
            continue
        
        meta = metadata.get(specimen_id, {})
        stratum_value = meta.get(stratum_column, 'unknown')
        
        host_call = row.get('host_call', '') or row.get('scientific_name', '') or ''
        if not host_call or host_call in ('', 'unassigned', 'no_host_signal'):
            continue
        
        mixed_by_stratum[stratum_value]['n'] += 1
        
        # Check for mixed host (may contain multiple hosts)
        mixed_status = row.get('mixed_status', '') or ''
        if mixed_status == 'mixed_host':
            mixed_by_stratum[stratum_value]['n_mixed'] += 1
    
    result = {}
    for stratum, counts in mixed_by_stratum.items():
        n = counts['n']
        n_mixed = counts['n_mixed']
        rate = n_mixed / n if n > 0 else 0.0
        result[stratum] = {
            'n': n,
            'n_mixed': n_mixed,
            'mixed_rate': round(rate, 6),
            'small_n_warning': n < 5,
        }
    
    return result


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
            'raw_p_value': round(raw_p, 6) if raw_p else None,
            'corrected_p_value': None,  # Filled by Holm correction
            'significant_005': None,  # Filled by Holm correction
            'significant_050': None,  # Filled by Holm correction
            'small_n_warning': (n1 < 5) or (n2 < 5),
            'test': 'fishers_exact',
            'test_family': label,
        })
    
    return comparisons


def compute_host_richness_by_stratum(host_call_rows, metadata, stratum_column):
    """Compute host richness (S) by stratum.
    
    Returns:
        dict: {stratum_value: {'n': n_mosquitoes, 'richness': S}}
    """
    richness_by_stratum = defaultdict(lambda: {'n': 0, 'hosts': set()})
    
    for row in host_call_rows:
        specimen_id = row.get('specimen_id') or row.get('sample_id', '')
        control_status = row.get('control_status', '') or ''
        
        if control_status in ('extraction_blank', 'pcr_negative', 'positive_control'):
            continue
        
        meta = metadata.get(specimen_id, {})
        stratum_value = meta.get(stratum_column, 'unknown')
        
        host_call = row.get('host_call', '') or row.get('scientific_name', '') or ''
        if not host_call or host_call in ('', 'unassigned', 'no_host_signal'):
            continue
        
        richness_by_stratum[stratum_value]['n'] += 1
        richness_by_stratum[stratum_value]['hosts'].add(host_call.lower())
    
    result = {}
    for stratum, data in richness_by_stratum.items():
        result[stratum] = {
            'n': data['n'],
            'richness': len(data['hosts']),
            'hosts': '|'.join(sorted(data['hosts'])),
            'small_n_warning': data['n'] < 5,
        }
    
    return result


def main():
    parser = argparse.ArgumentParser(
        description='Compute host ecology statistical comparisons for HÆMA blood-meal data'
    )
    parser.add_argument('--host-call-table', required=True,
                        help='Path to rambo_host_call_table.tsv or final host-call table')
    parser.add_argument('--metadata', required=True,
                        help='Path to samplesheet metadata CSV/TSV')
    parser.add_argument('--zone-column', default='collection_region',
                        help='Samplesheet column for ecological zone stratification')
    parser.add_argument('--species-column', default='sibling_species',
                        help='Samplesheet column for sibling species stratification')
    parser.add_argument('--season-column', default='season',
                        help='Samplesheet column for season stratification (wet/dry)')
    parser.add_argument('--output-dir', required=True,
                        help='Output directory for comparison tables')
    parser.add_argument('--exploratory-note', default='All statistical comparisons are exploratory unless pre-specified as primary in the study protocol.',
                        help='Note to include in output about exploratory nature')
    args = parser.parse_args()
    
    # Read input data
    host_call_rows = read_tsv(args.host_call_table)
    metadata = read_metadata(args.metadata)
    
    # Ensure output directory exists
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    # === HBI comparisons by ecological zone ===
    hbi_by_zone = compute_hbi_by_stratum(host_call_rows, metadata, args.zone_column)
    hbi_zone_comparisons = pairwise_comparisons(hbi_by_zone, 'n_human', 'hbi_zone')
    
    # Apply Holm correction
    holm_hbi = holm_correction([(i, c['raw_p_value']) for i, c in enumerate(hbi_zone_comparisons) if c['raw_p_value'] is not None])
    for idx, (orig_idx, raw_p, corr_p) in enumerate(holm_hbi):
        if orig_idx < len(hbi_zone_comparisons):
            hbi_zone_comparisons[orig_idx]['corrected_p_value'] = round(corr_p, 6)
            hbi_zone_comparisons[orig_idx]['significant_005'] = corr_p < 0.05
            hbi_zone_comparisons[orig_idx]['significant_050'] = corr_p < 0.050
    
    # Write HBI zone comparisons
    if hbi_zone_comparisons:
        fieldnames = list(hbi_zone_comparisons[0].keys())
        outpath = Path(args.output_dir) / 'pairwise_hbi_comparisons.tsv'
        with open(outpath, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in hbi_zone_comparisons:
                writer.writerow(row)
        print(f"Wrote HBI zone comparisons: {outpath}")
    
    # === HBI comparisons by sibling species ===
    hbi_by_species = compute_hbi_by_stratum(host_call_rows, metadata, args.species_column)
    hbi_species_comparisons = pairwise_comparisons(hbi_by_species, 'n_human', 'hbi_species')
    
    holm_hbi_sp = holm_correction([(i, c['raw_p_value']) for i, c in enumerate(hbi_species_comparisons) if c['raw_p_value'] is not None])
    for idx, (orig_idx, raw_p, corr_p) in enumerate(holm_hbi_sp):
        if orig_idx < len(hbi_species_comparisons):
            hbi_species_comparisons[orig_idx]['corrected_p_value'] = round(corr_p, 6)
            hbi_species_comparisons[orig_idx]['significant_005'] = corr_p < 0.05
            hbi_species_comparisons[orig_idx]['significant_050'] = corr_p < 0.050
    
    if hbi_species_comparisons:
        fieldnames = list(hbi_species_comparisons[0].keys())
        outpath = Path(args.output_dir) / 'pairwise_hbi_species_comparisons.tsv'
        with open(outpath, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in hbi_species_comparisons:
                writer.writerow(row)
        print(f"Wrote HBI species comparisons: {outpath}")
    
    # === Mixed-feeding comparisons by zone ===
    mixed_by_zone = compute_mixed_feeding_by_stratum(host_call_rows, metadata, args.zone_column)
    mixed_zone_comparisons = pairwise_comparisons(mixed_by_zone, 'n_mixed', 'mixed_feeding_zone')
    
    holm_mixed = holm_correction([(i, c['raw_p_value']) for i, c in enumerate(mixed_zone_comparisons) if c['raw_p_value'] is not None])
    for idx, (orig_idx, raw_p, corr_p) in enumerate(holm_mixed):
        if orig_idx < len(mixed_zone_comparisons):
            mixed_zone_comparisons[orig_idx]['corrected_p_value'] = round(corr_p, 6)
            mixed_zone_comparisons[orig_idx]['significant_005'] = corr_p < 0.05
            mixed_zone_comparisons[orig_idx]['significant_050'] = corr_p < 0.050
    
    if mixed_zone_comparisons:
        fieldnames = list(mixed_zone_comparisons[0].keys())
        outpath = Path(args.output_dir) / 'pairwise_mixed_feeding_comparisons.tsv'
        with open(outpath, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in mixed_zone_comparisons:
                writer.writerow(row)
        print(f"Wrote mixed-feeding comparisons: {outpath}")
    
    # === Host richness by zone ===
    richness_by_zone = compute_host_richness_by_stratum(host_call_rows, metadata, args.zone_column)
    
    # Write richness summary
    richness_rows = []
    for stratum, data in richness_by_zone.items():
        richness_rows.append({
            'stratum': stratum,
            'n_mosquitoes': data['n'],
            'host_richness_S': data['richness'],
            'hosts': data['hosts'],
            'small_n_warning': data['small_n_warning'],
        })
    
    if richness_rows:
        outpath = Path(args.output_dir) / 'host_richness_by_zone.tsv'
        with open(outpath, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=list(richness_rows[0].keys()), delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in richness_rows:
                writer.writerow(row)
        print(f"Wrote host richness summary: {outpath}")
    
    # === Write summary ===
    summary = {
        'total_samples_analyzed': len(host_call_rows),
        'samples_excluded_controls': sum(1 for r in host_call_rows if r.get('control_status', '') in ('extraction_blank', 'pcr_negative', 'positive_control')),
        'hbi_zone_comparisons': len(hbi_zone_comparisons),
        'hbi_species_comparisons': len(hbi_species_comparisons),
        'mixed_feeding_comparisons': len(mixed_zone_comparisons),
        'strata_hbi_zone': list(hbi_by_zone.keys()),
        'strata_hbi_species': list(hbi_by_species.keys()),
        'strata_mixed_zone': list(mixed_by_zone.keys()),
        'strata_richness_zone': list(richness_by_zone.keys()),
        'exploratory_note': args.exploratory_note,
        'notes': [
            'All tests are Fisher\'s exact (two-tailed).',
            'Holm-Bonferroni correction applied within each test family.',
            'Small-n warnings (n < 5 per stratum) indicate low power.',
            'Kruskal-Wallis test for host richness is planned but not yet implemented.',
            'scipy is required for exact Fisher\'s exact test; without it, chi-square approximation is used (n >= 10 only).',
        ]
    }
    
    summary_path = Path(args.output_dir) / 'host_use_statistical_tests_summary.json'
    with open(summary_path, 'w') as fh:
        json.dump(summary, fh, indent=2)
        fh.write('\n')
    print(f"Wrote summary: {summary_path}")
    
    # Print summary to stdout
    print(f"\nHost ecology comparisons completed:")
    print(f"  HBI zone comparisons: {len(hbi_zone_comparisons)}")
    print(f"  HBI species comparisons: {len(hbi_species_comparisons)}")
    print(f"  Mixed-feeding comparisons: {len(mixed_zone_comparisons)}")
    print(f"  Host richness strata: {len(richness_by_zone)}")


if __name__ == '__main__':
    main()
