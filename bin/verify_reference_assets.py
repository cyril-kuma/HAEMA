#!/usr/bin/env python3
"""Verify bundled reference assets: file checksums and curated taxonomy sidecar fields.

Used by the pre-release validation script and CI, and by anyone who edits the curated panel.
Exits non-zero (with clear messages) on any failure. Does NOT fabricate data — it only checks
that what is shipped is internally consistent and complete.

Usage:
    verify_reference_assets.py --assets-dir assets/references
"""
import argparse
import csv
import hashlib
import sys
from pathlib import Path

REQUIRED_SIDECAR_FIELDS = ["seqid", "taxid", "scientific_name", "rank", "source_accession"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_checksums(assets_dir, errors):
    manifest = assets_dir / "CHECKSUMS.sha256"
    if not manifest.exists():
        errors.append(f"Missing checksum manifest: {manifest}")
        return
    for line in manifest.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        expected, _, name = line.partition(" *") if " *" in line else line.partition("  ")
        name = name.strip()
        target = assets_dir / name
        if not target.exists():
            errors.append(f"Checksum manifest lists a missing file: {name}")
            continue
        actual = sha256(target)
        if actual != expected.strip():
            errors.append(f"Checksum mismatch for {name}: expected {expected.strip()[:12]}…, got {actual[:12]}…")


def verify_sidecar(assets_dir, errors):
    sidecar = assets_dir / "vertebrate_dna_ref_panel.taxonomy.tsv"
    fasta = assets_dir / "vertebrate_dna_ref_panel.fasta"
    if not sidecar.exists():
        errors.append(f"Missing taxonomy sidecar: {sidecar}")
        return
    fasta_ids = set()
    if fasta.exists():
        fasta_ids = {ln[1:].split()[0] for ln in fasta.read_text().splitlines() if ln.startswith(">")}
    with sidecar.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing_cols = [c for c in REQUIRED_SIDECAR_FIELDS if c not in (reader.fieldnames or [])]
        if missing_cols:
            errors.append(f"Sidecar missing required column(s): {', '.join(missing_cols)}")
            return
        n = 0
        for i, row in enumerate(reader, start=2):
            n += 1
            for field in REQUIRED_SIDECAR_FIELDS:
                if not (row.get(field) or "").strip():
                    errors.append(f"Sidecar line {i}: empty required field '{field}' (seqid={row.get('seqid','?')})")
            taxid = (row.get("taxid") or "").strip()
            if taxid and (not taxid.isdigit() or taxid == "0"):
                errors.append(f"Sidecar line {i}: taxid '{taxid}' is not a positive integer")
            seqid = (row.get("seqid") or "").strip()
            if fasta_ids and seqid and seqid not in fasta_ids:
                errors.append(f"Sidecar line {i}: seqid not found in panel FASTA: {seqid}")
    # coverage: every FASTA record should have a sidecar taxonomy row (re-read; small file)
    with sidecar.open(newline="") as fh:
        sidecar_ids = {(r.get("seqid") or "").strip() for r in csv.DictReader(fh, delimiter="\t")}
    for fid in sorted(fasta_ids - sidecar_ids):
        errors.append(f"Panel FASTA record has no sidecar taxonomy row: {fid}")


def main():
    ap = argparse.ArgumentParser(description="Verify bundled reference assets (checksums + sidecar fields).")
    ap.add_argument("--assets-dir", default="assets/references")
    args = ap.parse_args()
    assets_dir = Path(args.assets_dir)
    errors = []
    if not assets_dir.is_dir():
        print(f"ERROR: assets dir not found: {assets_dir}", file=sys.stderr)
        sys.exit(2)
    verify_checksums(assets_dir, errors)
    verify_sidecar(assets_dir, errors)
    if errors:
        print(f"Reference-asset verification FAILED ({len(errors)} issue(s)):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)
    print("Reference-asset verification PASSED (checksums + sidecar fields OK).")


if __name__ == "__main__":
    main()
