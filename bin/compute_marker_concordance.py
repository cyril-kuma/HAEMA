#!/usr/bin/env python3
"""Compute multi-marker concordance for HÆMA blood-meal host calls.

For each sample/specimen with host calls from two or more markers, assess whether
host assignments agree at species or genus level. Discordant calls between CytB
and COI for the same sample may indicate NUMT co-amplification, reference database
gaps, or genuine within-sample complexity.

Output concordance_status values:
  - full_species_agreement: All markers agree at species level
  - genus_agreement: All markers agree at genus level but not species
  - discordant: Markers disagree at genus or species level
  - single_marker_only: Only one marker has a host call
  - no_marker_signal: No markers have host calls
  - ambiguous_lca_only: All calls are LCA-resolved (genus-level or higher)

This feeds into:
  1. Final host-call confidence
  2. QC report
  3. Visualisation (concordance heatmap)
  4. Final thesis/manuscript output
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_tsv(path):
    """Read a TSV file and return list of dicts."""
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return []
    with p.open(newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def norm_taxon(name):
    """Normalize taxon name for comparison: lowercase, strip, collapse whitespace."""
    return ' '.join((name or '').strip().lower().replace('_', ' ').split())


def compute_concordance(host_call_rows):
    """Compute concordance across markers for each sample.
    
    Args:
        host_call_rows: List of dicts from rambo_host_call_table.tsv or similar
        
    Returns:
        List of concordance dicts keyed by specimen_id
    """
    # Group host calls by specimen_id
    calls_by_specimen = defaultdict(list)
    for row in host_call_rows:
        specimen_id = row.get('specimen_id') or row.get('sample_id') or row.get('sample_uid', '')
        if not specimen_id:
            continue
        marker = row.get('marker', '')
        host_call = row.get('host_call', '') or row.get('scientific_name', '') or ''
        if host_call and marker:
            calls_by_specimen[specimen_id].append({
                'marker': marker,
                'host_call': host_call,
                'confidence': row.get('confidence', ''),
                'pident': float(row.get('pident', 0) or 0),
                'lca_taxon': row.get('lca_taxon', '') or '',
                'top_hit': row.get('top_hit', '') or '',
                'assignment_method': row.get('assignment_method', ''),
                'numt_risk': row.get('numt_risk', ''),
            })
    
    concordance_results = []
    
    for specimen_id, calls in sorted(calls_by_specimen.items()):
        # Filter to markers with confident host calls
        confident_calls = [c for c in calls if c['host_call'] and c['confidence'] in ('high', 'medium', 'low')]
        
        if len(confident_calls) == 0:
            concordance_results.append({
                'specimen_id': specimen_id,
                'sample_id': calls[0].get('sample_id', '') if calls else '',
                'markers_with_signal': [],
                'host_calls_by_marker': {},
                'species_level_concordance': False,
                'genus_level_concordance': False,
                'concordance_status': 'no_marker_signal',
                'discordance_reason': '',
                'possible_numt_flag': False,
                'possible_mixed_meal_flag': False,
            })
            continue
        
        markers_with_signal = [c['marker'] for c in confident_calls]
        host_calls_by_marker = {c['marker']: c['host_call'] for c in confident_calls}
        
        # Check species-level concordance
        species_names = [norm_taxon(c['host_call']) for c in confident_calls]
        species_level_agree = len(set(species_names)) == 1
        
        # Check genus-level concordance
        genus_names = [norm_taxon(' '.join(c['host_call'].split()[:1])) for c in confident_calls]
        genus_level_agree = len(set(genus_names)) == 1
        
        # Check for ambiguous LCA-only calls
        lca_only = all(c['lca_taxon'] and c['assignment_method'] in ('taxid_lca', 'conservative_lca') for c in confident_calls)
        
        # Determine concordance status
        if len(confident_calls) == 1:
            concordance_status = 'single_marker_only'
            discordance_reason = ''
        elif species_level_agree:
            concordance_status = 'full_species_agreement'
            discordance_reason = ''
        elif genus_level_agree:
            concordance_status = 'genus_agreement'
            discordance_reason = 'All markers agree at genus level but not species'
        elif lca_only:
            concordance_status = 'ambiguous_lca_only'
            discordance_reason = 'All calls are LCA-resolved (genus-level or higher)'
        else:
            concordance_status = 'discordant'
            discordance_reason = 'Markers disagree at genus or species level'
        
        # Check for possible NUMT flag
        # NUMT risk: discordance between COI markers and CytB
        has_cytb = any(c['marker'] == 'cyt_b' for c in confident_calls)
        has_coi = any(c['marker'].startswith('co1') for c in confident_calls)
        possible_numt = has_cytb and has_coi and not genus_level_agree
        
        # Check for possible mixed meal flag
        # Mixed meal: different hosts called by different markers for same sample
        possible_mixed = len(set(species_names)) > 1 and not species_level_agree
        
        concordance_results.append({
            'specimen_id': specimen_id,
            'sample_id': calls[0].get('sample_id', '') if calls else '',
            'markers_with_signal': '|'.join(markers_with_signal),
            'host_calls_by_marker': json.dumps(host_calls_by_marker),
            'species_level_concordance': species_level_agree,
            'genus_level_concordance': genus_level_agree,
            'concordance_status': concordance_status,
            'discordance_reason': discordance_reason,
            'possible_numt_flag': possible_numt,
            'possible_mixed_meal_flag': possible_mixed,
        })
    
    return concordance_results


def compute_cohens_kappa(concordance_results):
    """Compute weighted Cohen's Kappa across all markers.
    
    Simplified: compares species-level agreement vs chance agreement.
    """
    # Count agreements and disagreements
    full_agree = sum(1 for r in concordance_results if r['concordance_status'] == 'full_species_agreement')
    genus_agree = sum(1 for r in concordance_results if r['concordance_status'] == 'genus_agreement')
    discordant = sum(1 for r in concordance_results if r['concordance_status'] == 'discordant')
    ambiguous = sum(1 for r in concordance_results if r['concordance_status'] in ('ambiguous_lca_only', 'single_marker_only'))
    no_signal = sum(1 for r in concordance_results if r['concordance_status'] == 'no_marker_signal')
    
    total = len(concordance_results)
    if total == 0:
        return 0.0, 0, 0
    
    # Observed agreement (species-level)
    p_obs = full_agree / total if total > 0 else 0
    
    # Expected agreement by chance (simplified)
    # P(species agree) * P(species agree) + P(genus agree) * P(genus agree) + ...
    p_chance = ((full_agree + genus_agree) / total) ** 2 + (discordant / total) ** 2 + (ambiguous + no_signal) / total * 0.1
    
    kappa = (p_obs - p_chance) / (1 - p_chance) if p_chance < 1 else 0
    
    return kappa, full_agree + genus_agree, discordant


def main():
    parser = argparse.ArgumentParser(
        description='Compute multi-marker concordance for HÆMA blood-meal host calls'
    )
    parser.add_argument('--host-call-table', required=True,
                        help='Path to rambo_host_call_table.tsv or equivalent')
    parser.add_argument('--endpoint-manifest', required=True,
                        help='Path to bloodmeal_master_endpoint.tsv')
    parser.add_argument('--output', required=True,
                        help='Output concordance TSV path')
    parser.add_argument('--summary-output', default='',
                        help='Output concordance summary JSON path')
    args = parser.parse_args()
    
    # Read input data
    host_call_rows = read_tsv(args.host_call_table)
    endpoint_rows = read_tsv(args.endpoint_manifest)
    
    # Compute concordance
    concordance_results = compute_concordance(host_call_rows)
    
    # Compute Cohen's Kappa
    kappa, n_agree, n_disagree = compute_cohens_kappa(concordance_results)
    
    # Write concordance TSV
    if concordance_results:
        fieldnames = list(concordance_results[0].keys())
        with open(args.output, 'w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
            writer.writeheader()
            for row in concordance_results:
                writer.writerow(row)
    
    # Write summary JSON
    if args.summary_output:
        summary = {
            'total_samples': len(concordance_results),
            'full_species_agreement': sum(1 for r in concordance_results if r['concordance_status'] == 'full_species_agreement'),
            'genus_agreement': sum(1 for r in concordance_results if r['concordance_status'] == 'genus_agreement'),
            'discordant': sum(1 for r in concordance_results if r['concordance_status'] == 'discordant'),
            'single_marker_only': sum(1 for r in concordance_results if r['concordance_status'] == 'single_marker_only'),
            'no_marker_signal': sum(1 for r in concordance_results if r['concordance_status'] == 'no_marker_signal'),
            'ambiguous_lca_only': sum(1 for r in concordance_results if r['concordance_status'] == 'ambiguous_lca_only'),
            'cohens_kappa': round(kappa, 4),
            'n_agreeing': n_agree,
            'n_disagreeing': n_disagree,
            'numt_flags': sum(1 for r in concordance_results if r['possible_numt_flag']),
            'mixed_meal_flags': sum(1 for r in concordance_results if r['possible_mixed_meal_flag']),
            'notes': [
                'Concordance computed from host calls across markers for each specimen.',
                'Discordant calls may indicate NUMT co-amplification, reference gaps, or genuine mixed meals.',
                'Cohen\'s Kappa is a simplified measure; weighted kappa by marker count is recommended for production.',
                'NUMT flag is a caution indicator, not proof of NUMT contamination.',
            ]
        }
        with open(args.summary_output, 'w') as fh:
            json.dump(summary, fh, indent=2)
            fh.write('\n')
    
    # Print summary to stdout
    print(f"Multi-marker concordance computed:")
    print(f"  Total samples: {len(concordance_results)}")
    print(f"  Full species agreement: {summary['full_species_agreement'] if args.summary_output else sum(1 for r in concordance_results if r['concordance_status'] == 'full_species_agreement')}")
    print(f"  Genus agreement: {summary['genus_agreement'] if args.summary_output else sum(1 for r in concordance_results if r['concordance_status'] == 'genus_agreement')}")
    print(f"  Discordant: {summary['discordant'] if args.summary_output else sum(1 for r in concordance_results if r['concordance_status'] == 'discordant')}")
    print(f"  Cohen's Kappa: {kappa:.4f}")
    print(f"  Output: {args.output}")


if __name__ == '__main__':
    main()
