#!/usr/bin/env python3
"""Build a broad, reproducible blood-meal reference FASTA from a public source.

Normalises a RefSeq-mitochondrion-style FASTA (or any GenBank-deflined FASTA) into
a species-deflined FASTA suitable for the HÆMA `broad_blast` / `bold_aware` reference
modes, and writes a provenance JSON (source, version, build date, command, SHA-256,
sequence and taxon counts). The normalised deflines are `>Genus_species|accession`
(no whitespace) so that `parse_blast_assignments.py` recovers the host from the BLAST
`sseqid` (first `|`-delimited token) without needing `-parse_seqids`/taxids — mirroring
the curated panel convention. This keeps the broad database reproducible and offline.

Rationale: reference-database completeness is the primary limiting factor for amplicon
blood-meal identification (Channumsin 2021; Kipp 2023; Santos 2019). RefSeq
mitochondrion gives broad, versioned vertebrate (and other) mitogenome coverage of
both COI and CytB, catching hosts absent from the small curated panel, while remaining
a fixed, checksummed, offline database rather than a live API.

Example:
  python3 build_reference_db.py \
      --source-fasta mitochondrion.*.genomic.fna.gz \
      --source-name "NCBI RefSeq mitochondrion" --source-version "release 235" \
      --db-name refseq_mito_r235 \
      --output-fasta refseq_mito_r235.fasta \
      --output-provenance refseq_mito_r235.provenance.json
"""
import argparse
import datetime as _dt
import gzip
import hashlib
import json
import re
import sys
from pathlib import Path

# Words that mark the end of the organism name in a GenBank/RefSeq mitochondrion defline.
_STOP = re.compile(
    r"\b(mitochondrion|mitochondrial|complete|genome|isolate|voucher|strain|"
    r"haplotype|clone|DNA|chromosome|unplaced|cds|gene)\b",
    re.IGNORECASE,
)


def _open(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def normalise_organism(title):
    """Extract 'Genus_species' from a GenBank/RefSeq defline title (after the accession)."""
    # Drop the leading accession token.
    parts = title.split(None, 1)
    rest = parts[1] if len(parts) > 1 else parts[0]
    # Cut at the first stop-word or comma.
    cut = _STOP.search(rest)
    name = rest[: cut.start()] if cut else rest
    name = name.split(",")[0].strip()
    tokens = [t for t in re.split(r"\s+", name) if t]
    if not tokens:
        return ""
    # Keep genus + species epithet (+ subspecies epithet when the 3rd token is lowercase alpha).
    keep = tokens[:2]
    if len(tokens) >= 3 and re.fullmatch(r"[a-z-]+", tokens[2]):
        keep.append(tokens[2])
    genus_species = "_".join(keep)
    # Guard against non-binomial junk (e.g. 'UNVERIFIED').
    if not re.match(r"^[A-Z][A-Za-z.'-]+_[a-z]", genus_species):
        return ""
    return genus_species


def iter_records(paths):
    for path in paths:
        header = None
        seq = []
        with _open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    if header is not None:
                        yield header, "".join(seq)
                    header = line[1:]
                    seq = []
                else:
                    seq.append(line.strip())
            if header is not None:
                yield header, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-fasta", nargs="+", required=True, help="One or more source FASTA(.gz) files.")
    ap.add_argument("--output-fasta", required=True)
    ap.add_argument("--output-provenance", required=True)
    ap.add_argument("--db-name", required=True)
    ap.add_argument("--source-name", default="unknown")
    ap.add_argument("--source-version", default="unknown")
    ap.add_argument("--source-url", default="")
    ap.add_argument("--marker-coverage", default="mitogenome (COI + CytB)")
    ap.add_argument("--min-length", type=int, default=200, help="Drop sequences shorter than this.")
    args = ap.parse_args()

    out = Path(args.output_fasta)
    n_in = n_out = n_skipped_name = n_skipped_len = 0
    taxa = set()
    sha = hashlib.sha256()
    with out.open("w") as ofh:
        for header, seq in iter_records(args.source_fasta):
            n_in += 1
            accession = header.split()[0]
            organism = normalise_organism(header)
            if not organism:
                n_skipped_name += 1
                continue
            if len(seq) < args.min_length:
                n_skipped_len += 1
                continue
            defline = f">{organism}|{accession}"
            block = f"{defline}\n{seq}\n"
            ofh.write(block)
            sha.update(block.encode())
            taxa.add(organism)
            n_out += 1

    provenance = {
        "database_name": args.db_name,
        "database_type": "fasta",
        "database_source": args.source_name,
        "database_version": args.source_version,
        "source_url": args.source_url,
        "build_date": _dt.date.today().isoformat(),
        "build_command": "build_reference_db.py " + " ".join(sys.argv[1:]),
        "marker_coverage": args.marker_coverage,
        "number_of_sequences": n_out,
        "number_of_taxa": len(taxa),
        "min_length": args.min_length,
        "sequences_in": n_in,
        "sequences_dropped_no_binomial": n_skipped_name,
        "sequences_dropped_short": n_skipped_len,
        "output_fasta": str(out),
        "output_fasta_sha256": sha.hexdigest(),
        "defline_format": "Genus_species|accession",
        "notes": [
            "Deflines normalised to species|accession so host is recovered from BLAST sseqid.",
            "Reproducible, offline, versioned database; not a live API.",
        ],
    }
    Path(args.output_provenance).write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"[build_reference_db] {n_out} sequences / {len(taxa)} taxa written to {out}")
    print(f"[build_reference_db] dropped: no-binomial={n_skipped_name} short(<{args.min_length})={n_skipped_len}")
    print(f"[build_reference_db] sha256={sha.hexdigest()}")
    print(f"[build_reference_db] provenance -> {args.output_provenance}")


if __name__ == "__main__":
    main()
