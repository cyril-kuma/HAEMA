#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def read_json(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open() as handle:
        return json.load(handle)


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
        "databases": {
            "reference_fasta": parameters.get("reference_fasta", ""),
            "curated_reference_metadata": parameters.get("curated_reference_metadata", ""),
            "blast_db": parameters.get("blast_db", ""),
            "fallback_blast_db": parameters.get("fallback_blast_db", ""),
            "taxdump_dir": parameters.get("taxdump_dir", ""),
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
        ],
    }
    with Path(args.output).open("w") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")


if __name__ == "__main__":
    main()
