#!/usr/bin/env python3
import argparse
import json
from pathlib import Path
import hashlib


def read_json(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open() as handle:
        return json.load(handle)


def compute_sha256(path):
    """Compute SHA-256 checksum of a file, returning None if the file does not exist."""
    path = Path(path)
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()


def count_fasta_records(path):
    """Count the number of sequences in a FASTA file, returning None on error."""
    path = Path(path)
    if not path.exists():
        return None
    count = 0
    with path.open() as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def build_database_provenance(params):
    """Build comprehensive database provenance from pipeline parameters.
    
    Records for every database used:
    - database_name, database_type, database_source, database_path
    - database_version, download_date, build_date
    - number_of_sequences, number_of_taxa
    - marker, sha256_checksum, taxonomy_map_path, taxdump_path
    """
    provenance = {}
    
    # Curated panel
    ref_fasta = params.get('reference_fasta', '')
    if ref_fasta:
        checksum = compute_sha256(ref_fasta)
        seq_count = count_fasta_records(ref_fasta)
        curated_meta = params.get('curated_reference_metadata', '')
        provenance['curated_panel'] = {
            'database_name': 'HÆMA Curated Vertebrate Panel',
            'database_type': 'fasta',
            'database_source': 'bundled',
            'database_path': ref_fasta,
            'database_version': params.get('reference_version', 'unversioned'),
            'download_date': None,
            'build_date': None,
            'number_of_sequences': seq_count,
            'number_of_taxa': None,
            'marker': 'co1_short,co1_long,cyt_b',
            'sha256_checksum': checksum,
            'taxonomy_map_path': curated_meta if curated_meta else None,
            'taxdump_path': None,
            'notes': 'Ghana/West African peridomestic vertebrate hosts'
        }
    
    # Broad BLAST database
    blast_db = params.get('blast_db', '')
    if blast_db:
        db_path = blast_db
        checksum = compute_sha256(db_path)
        provenance['broad_blast_db'] = {
            'database_name': params.get('blast_db_label', 'broad_blast'),
            'database_type': 'blast_db',
            'database_source': params.get('blast_db_source', 'user_supplied'),
            'database_path': blast_db,
            'database_version': params.get('blast_db_version', 'unknown'),
            'download_date': params.get('blast_db_download_date', None),
            'build_date': params.get('blast_db_build_date', None),
            'number_of_sequences': None,
            'number_of_taxa': None,
            'marker': 'all',
            'sha256_checksum': checksum,
            'taxonomy_map_path': None,
            'taxdump_path': params.get('taxdump_dir', ''),
            'notes': 'Broad reference database for taxonomy assignment'
        }
    
    # Fallback BLAST database
    fallback_db = params.get('fallback_blast_db', '')
    if fallback_db:
        checksum = compute_sha256(fallback_db)
        provenance['fallback_blast_db'] = {
            'database_name': params.get('fallback_blast_db_label', 'ncbi_nt'),
            'database_type': 'blast_db',
            'database_source': 'ncbi_nt',
            'database_path': fallback_db,
            'database_version': params.get('fallback_blast_db_version', 'unknown'),
            'download_date': params.get('fallback_blast_db_download_date', None),
            'build_date': params.get('fallback_blast_db_build_date', None),
            'number_of_sequences': None,
            'number_of_taxa': None,
            'marker': 'all',
            'sha256_checksum': checksum,
            'taxonomy_map_path': None,
            'taxdump_path': params.get('taxdump_dir', ''),
            'notes': 'NCBI nt remote fallback for unresolved queries'
        }
    
    # BOLD-derived COI database
    bold_fasta = params.get('bold_fasta', '')
    if bold_fasta:
        checksum = compute_sha256(bold_fasta)
        bold_tax = params.get('bold_taxonomy', '')
        provenance['bold_coi'] = {
            'database_name': 'BOLD-derived COI Database',
            'database_type': 'fasta',
            'database_source': 'bold',
            'database_path': bold_fasta,
            'database_version': params.get('bold_version', 'unknown'),
            'download_date': params.get('bold_download_date', None),
            'build_date': params.get('bold_build_date', None),
            'number_of_sequences': count_fasta_records(bold_fasta),
            'number_of_taxa': None,
            'marker': 'co1_short,co1_long',
            'sha256_checksum': checksum,
            'taxonomy_map_path': bold_tax if bold_tax else None,
            'taxdump_path': None,
            'notes': 'BOLD-derived COI barcodes for Mode D (BOLD-aware COI)'
        }
    
    # Taxdump
    taxdump_dir = params.get('taxdump_dir', '')
    if taxdump_dir:
        provenance['taxdump'] = {
            'database_name': 'NCBI Taxonomy',
            'database_type': 'taxdump',
            'database_source': 'ncbi_taxonomy',
            'database_path': taxdump_dir,
            'database_version': params.get('taxdump_version', 'unknown'),
            'download_date': params.get('taxdump_download_date', None),
            'build_date': None,
            'number_of_sequences': None,
            'number_of_taxa': None,
            'marker': 'all',
            'sha256_checksum': None,
            'taxonomy_map_path': None,
            'taxdump_path': taxdump_dir,
            'notes': 'NCBI nodes.dmp + names.dmp for taxid-backed LCA'
        }
    
    # Warn about missing provenance
    warnings = []
    if not ref_fasta and not blast_db and not bold_fasta:
        warnings.append('No reference database FASTA or BLAST database configured.')
    if params.get('enable_ncbi_remote_fallback', False) and not fallback_db:
        warnings.append('NCBI remote fallback enabled but --fallback_blast_db is not set.')
    if params.get('enable_bold_fallback', False) and not bold_fasta:
        warnings.append('BOLD fallback enabled but --bold_fasta is not set.')
    
    provenance['_warnings'] = warnings
    return provenance


def build_numt_risk(params):
    """Build NUMT risk reporting from parameters or defaults."""
    marker_risk_str = params.get('marker_numt_risk', '')
    if marker_risk_str:
        try:
            return json.loads(marker_risk_str)
        except (json.JSONDecodeError, TypeError):
            pass
    # Default NUMT risk values
    return {
        'co1_short': {'risk': 'moderate', 'rationale': 'Borderline above <200 bp risk threshold (Reeves et al. 2018)'},
        'co1_long': {'risk': 'low', 'rationale': 'Longer amplicon reduces NUMT co-amplification probability'},
        'cyt_b': {'risk': 'low', 'rationale': 'Mitochondrial multicopy advantage confirmed (Hadj-Henni et al. 2015)'}
    }


def main():
    parser = argparse.ArgumentParser(description="Build a machine-readable HÆMA pipeline run manifest")
    parser.add_argument("--endpoint-manifest", required=True)
    parser.add_argument("--input-validation-report", required=True)
    parser.add_argument("--production-preflight-report", required=True)
    parser.add_argument("--pipeline-version", default="")
    parser.add_argument("--workflow-run-name", default="")
    parser.add_argument("--workflow-session-id", default="")
    parser.add_argument("--workflow-profile", default="")
    parser.add_argument("--workflow-command-line", default="")
    parser.add_argument("--parameters-json", default="")
    parser.add_argument("--outputs-json", default="")
    parser.add_argument("--parameters-json-file", default="")
    parser.add_argument("--outputs-json-file", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    parameters_text = Path(args.parameters_json_file).read_text() if args.parameters_json_file else args.parameters_json
    outputs_text = Path(args.outputs_json_file).read_text() if args.outputs_json_file else args.outputs_json
    try:
        parameters = json.loads(parameters_text)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Could not parse parameters JSON: {exc}")
    try:
        outputs = json.loads(outputs_text)
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Could not parse outputs JSON: {exc}")

    # Build comprehensive database provenance
    database_provenance = build_database_provenance(parameters)
    
    # Build NUMT risk reporting
    numt_risk = build_numt_risk(parameters)
    
    # Build reference mode summary
    reference_mode = parameters.get('reference_mode', 'curated_panel')
    taxonomy_strategy = parameters.get('taxonomy_strategy', 'curated_then_fallback')
    
    manifest = {
        "pipeline": {
            "name": "haema/bloodmeal-metabarcoding",
            "version": args.pipeline_version,
            "run_name": args.workflow_run_name,
            "session_id": args.workflow_session_id,
            "profile": args.workflow_profile,
            "command_line": args.workflow_command_line,
        },
        "parameters": parameters,
        "reference_database": {
            "mode": reference_mode,
            "strategy": taxonomy_strategy,
            "provenance": database_provenance,
            "numt_risk": numt_risk,
            "markers": {
                "co1_short": {"target_size_bp": 234, "species_identity_threshold": parameters.get('coi_species_identity_threshold', 98.0)},
                "co1_long": {"target_size_bp": 359, "species_identity_threshold": parameters.get('coi_species_identity_threshold', 98.0)},
                "cyt_b": {"target_size_bp": 359, "species_identity_threshold": parameters.get('cytb_species_identity_threshold', 97.0)}
            },
            "identity_thresholds": {
                "global_min_identity": parameters.get('min_blast_identity', 97.0),
                "global_min_coverage": parameters.get('min_blast_coverage', 80.0),
                "top_bitscore_delta": parameters.get('top_bitscore_delta', 2.0),
                "evalue_threshold": parameters.get('blast_evalue', '1e-20'),
                "assignment_method": parameters.get('taxonomy_assignment_method', 'conservative_lca')
            },
            "fallback_chain": {
                "curated_panel": parameters.get('enable_curated_panel_check', True),
                "ncbi_remote": parameters.get('enable_ncbi_remote_fallback', True),
                "bold_coi": parameters.get('enable_bold_fallback', False),
                "bold_mode": parameters.get('bold_mode', 'local_fasta')
            }
        },
        "databases": {
            "reference_fasta": parameters.get("reference_fasta", ""),
            "curated_reference_metadata": parameters.get("curated_reference_metadata", ""),
            "blast_db": parameters.get("blast_db", ""),
            "fallback_blast_db": parameters.get("fallback_blast_db", ""),
            "taxdump_dir": parameters.get("taxdump_dir", ""),
            "bold_fasta": parameters.get("bold_fasta", ""),
            "bold_taxonomy": parameters.get("bold_taxonomy", ""),
            # Content hashes of the resolved reference assets (computed host-side at run time), so a
            # run is reproducible/auditable even if the panel or sidecar is edited or overridden.
            "reference_checksums_sha256": parameters.get("reference_checksums_sha256", {}),
        },
        "outputs": outputs,
        "input_validation": read_json(args.input_validation_report),
        "production_preflight": read_json(args.production_preflight_report),
        "endpoint_manifest": read_json(args.endpoint_manifest),
        "notes": [
            "Rows are flagged rather than silently removed.",
            "Advanced demux, Medaka, phyloseq, decontam, and curated taxid LCA are parameter-controlled optional stages.",
            "Per-process versions.yml files are published with their respective result directories.",
            "Database provenance is recorded in the 'reference_database' section above.",
            "NUMT risk is a per-marker caution flag, not a hard filter.",
            "Read fractions are supporting evidence only; they do not represent proportional ingested blood volumes."
        ],
    }
    with Path(args.output).open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
