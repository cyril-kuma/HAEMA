#!/usr/bin/env python3
"""End-to-end test of taxid-LCA assignment for the curated panel.

Proves that taxid-backed LCA works WITHOUT native BLAST staxids: the curated taxonomy sidecar
backfills taxids, and a synthetic taxdump resolves the lowest common ancestor. A query that hits two
divergent reference taxa is conservatively escalated to their common ancestor.

Run from the repo root:  python3 tests/test_taxid_assignment.py
"""
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PARSE = REPO / "bin" / "parse_blast_assignments.py"

# Synthetic taxdump: Homo sapiens (9606) and Bos taurus (9913) under Mammalia (40674).
NODES = [
    ("1", "1", "no rank"), ("40674", "1", "class"),
    ("9604", "40674", "family"), ("9605", "9604", "genus"), ("9606", "9605", "species"),
    ("9895", "40674", "family"), ("9903", "9895", "genus"), ("9913", "9903", "species"),
]
NAMES = {"1": "root", "40674": "Mammalia", "9606": "Homo sapiens", "9913": "Bos taurus",
         "9605": "Homo", "9604": "Hominidae", "9903": "Bos", "9895": "Bovidae"}

HUMAN = "Human_test_host|test_reference|synthetic"
CATTLE = "Cattle_test_host|test_reference|synthetic"

COUNTS = (
    "sample_uid\trun_id\tsample_id\tbarcode_id\tmarker\tcluster_id\tasv_id\tsequence\tcount\tfraction\tretained\n"
    "S1\tR\tS1\tbc1\tcyt_b\tall\tASV0001\tACGT\t50\t1.0\ttrue\n"
    "S1\tR\tS1\tbc1\tcyt_b\tall\tASV0002\tACGT\t40\t1.0\ttrue\n"
)
# outfmt 6: qseqid sseqid pident length qlen slen evalue bitscore stitle staxids  (staxids empty -> backfilled)
BLAST = (
    f"ASV0001\t{HUMAN}\t99.0\t400\t400\t400\t1e-50\t700\tN/A\t\n"          # single taxon -> Homo sapiens
    f"ASV0002\t{HUMAN}\t98.0\t400\t400\t400\t1e-48\t690\tN/A\t\n"          # two divergent taxa within delta
    f"ASV0002\t{CATTLE}\t98.0\t400\t400\t400\t1e-48\t689\tN/A\t\n"         # -> escalate to Mammalia
)
SIDECAR = (
    "seqid\ttaxid\tscientific_name\trank\n"
    f"{HUMAN}\t9606\tHomo sapiens\tspecies\n"
    f"{CATTLE}\t9913\tBos taurus\tspecies\n"
)


def main():
    failures = []
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\n" for t, p, r in NODES))
        (d / "names.dmp").write_text("".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in NAMES.items()))
        (d / "counts.tsv").write_text(COUNTS)
        (d / "blast.tsv").write_text(BLAST)
        (d / "sidecar.tsv").write_text(SIDECAR)
        out = d / "tax.tsv"
        r = subprocess.run(
            [sys.executable, str(PARSE), "--counts", str(d / "counts.tsv"), "--blast", str(d / "blast.tsv"),
             "--assignment-method", "taxid_lca", "--taxdump-dir", str(d),
             "--reference-taxonomy", str(d / "sidecar.tsv"),
             "--min-identity", "90", "--min-coverage", "80", "--top-bitscore-delta", "2.0",
             "--output", str(out)],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            failures.append(f"parse_blast_assignments exited {r.returncode}")
        else:
            import csv
            rows = {row["asv_id"]: row for row in csv.DictReader(out.open(), delimiter="\t")}
            a1, a2 = rows.get("ASV0001", {}), rows.get("ASV0002", {})
            if a1.get("host_assignment") != "Homo sapiens":
                failures.append(f"ASV0001 host {a1.get('host_assignment')!r}, want 'Homo sapiens' (sidecar-backfilled taxid LCA)")
            if a1.get("lca_taxid") != "9606":
                failures.append(f"ASV0001 lca_taxid {a1.get('lca_taxid')!r}, want '9606'")
            if a2.get("host_assignment") != "Mammalia" or a2.get("lca_taxid") != "40674":
                failures.append(f"ASV0002 should escalate to Mammalia/40674, got {a2.get('host_assignment')!r}/{a2.get('lca_taxid')!r}")
            if a2.get("taxon_rank") != "class":
                failures.append(f"ASV0002 rank {a2.get('taxon_rank')!r}, want 'class'")

    if failures:
        print("taxid-assignment tests FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        sys.exit(1)
    print("taxid-assignment tests PASSED (curated-panel taxid LCA via sidecar; conservative escalation)")


if __name__ == "__main__":
    main()
