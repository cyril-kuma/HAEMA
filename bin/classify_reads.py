#!/usr/bin/env python3
"""Reference-guided read classification (ONT-adapted Logue 2016 'classify-then-count').

Classifies each read against a marker/host reference with minimap2 (ONT preset), assigns it to a
host taxon by best hit passing identity/coverage floors, then *counts reads per taxon* and denoises
by **abundance, not geometry** — replacing the fragile unsupervised UMAP/HDBSCAN k-mer clustering.
Output is the same schema RAMBO/ecological-indices already consume, so everything downstream is
unchanged. Vector/non-host self-hits are excluded by genus (see --non-host-genera).

Inputs (either):
  * --reads reads.fastq[.gz] --reference ref.fasta  (runs minimap2 itself), or
  * --paf alignments.paf  (pre-computed minimap2 -x map-ont PAF; lets the caller run minimap2 in a
    container).

Reference deflines must start with the taxon, `Genus_species|...` (curated panel / RefSeq-mito
build convention), so the host is recovered from the PAF target name.
"""
import argparse
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

NON_HOST_DEFAULT = {
    "Anopheles", "Culex", "Aedes", "Ochlerotatus", "Culiseta", "Mansonia",
    "Culicidae", "Culicinae", "Anophelinae", "Culicoides", "Phlebotomus",
    "Lutzomyia", "Simulium", "Glossina", "Ixodes", "Rhipicephalus", "Amblyomma",
}
UNRESOLVED = {"", "unassigned", "ambiguous", "unresolved"}


def taxon_from_target(tname):
    return tname.split("|")[0].replace("_", " ").strip()


def is_host(taxon, non_host_genera):
    if taxon.lower() in UNRESOLVED:
        return False
    genus = taxon.split()[0] if taxon else ""
    return genus not in non_host_genera


def run_minimap2(reads, reference, threads, minimap2="minimap2"):
    cmd = [minimap2, "-x", "map-ont", "--secondary=no", "-t", str(threads), reference, reads]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"minimap2 failed: {proc.stderr[:500]}")
    return proc.stdout.splitlines()


def best_hits_from_paf(paf_lines, min_identity, min_coverage):
    """Return {read_id: taxon} for the best hit per read passing identity/coverage floors."""
    best = {}  # qname -> (nmatch, taxon)
    for line in paf_lines:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        qname, qlen = f[0], int(f[1])
        qstart, qend = int(f[2]), int(f[3])
        tname = f[5]
        nmatch, alnlen = int(f[9]), int(f[10])
        if alnlen == 0 or qlen == 0:
            continue
        identity = 100.0 * nmatch / alnlen
        coverage = 100.0 * (qend - qstart) / qlen
        if identity < min_identity or coverage < min_coverage:
            continue
        if qname not in best or nmatch > best[qname][0]:
            best[qname] = (nmatch, taxon_from_target(tname))
    return {q: t for q, (_, t) in best.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--reads", help="reads FASTQ(.gz); runs minimap2 itself")
    src.add_argument("--paf", help="pre-computed minimap2 -x map-ont PAF")
    ap.add_argument("--reference", help="reference FASTA (required with --reads)")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--sample-uid", default="sample")
    ap.add_argument("--marker", default="NA")
    ap.add_argument("--min-identity", type=float, default=85.0, help="min %% identity (ONT-appropriate)")
    ap.add_argument("--min-coverage", type=float, default=60.0, help="min %% query coverage")
    ap.add_argument("--min-reads-per-host", type=int, default=3)
    ap.add_argument("--min-host-fraction", type=float, default=0.01)
    ap.add_argument("--global-min-count", type=int, default=1,
                    help="drop a taxon supported by <= this many reads (Logue global filter)")
    ap.add_argument("--non-host-genera", default="")
    ap.add_argument("--output-counts", required=True)
    ap.add_argument("--output-per-read", default="")
    args = ap.parse_args()

    non_host = ({g.strip() for g in args.non_host_genera.split(",") if g.strip()}
                if args.non_host_genera else NON_HOST_DEFAULT)

    if args.paf:
        paf_lines = Path(args.paf).read_text().splitlines()
    else:
        if not args.reference:
            sys.exit("--reference required with --reads")
        paf_lines = run_minimap2(args.reads, args.reference, args.threads, args.minimap2)

    read_taxon = best_hits_from_paf(paf_lines, args.min_identity, args.min_coverage)

    counts = defaultdict(int)
    for taxon in read_taxon.values():
        counts[taxon] += 1
    total_classified = sum(counts.values())

    # Abundance denoise (Logue-style: count-based, not geometry-based).
    kept = {}
    for taxon, c in counts.items():
        if not is_host(taxon, non_host):
            continue
        if c <= args.global_min_count:
            continue
        frac = c / total_classified if total_classified else 0.0
        if c >= args.min_reads_per_host and frac >= args.min_host_fraction:
            kept[taxon] = (c, frac)

    with open(args.output_counts, "w") as fh:
        fh.write("sample_uid\tmarker\thost_assignment\tcount\tfraction\tassignment_method\n")
        for taxon, (c, frac) in sorted(kept.items(), key=lambda kv: -kv[1][0]):
            fh.write(f"{args.sample_uid}\t{args.marker}\t{taxon}\t{c}\t{frac:.4f}\tread_classification\n")

    if args.output_per_read:
        with open(args.output_per_read, "w") as fh:
            fh.write("read_id\thost_assignment\n")
            for rid, taxon in read_taxon.items():
                fh.write(f"{rid}\t{taxon}\n")

    mixed = len(kept) >= 2
    print(f"[classify_reads] {args.sample_uid}/{args.marker}: {total_classified} classified reads; "
          f"kept hosts={sorted(kept)}; call={'mixed' if mixed else ('single' if kept else 'none')}")


if __name__ == "__main__":
    main()
