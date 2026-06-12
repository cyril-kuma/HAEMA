#!/usr/bin/env python3
"""Unit test for taxdump-backed LCA in bin/parse_blast_assignments.py using a tiny synthetic tree.

Run from the repo root:  python3 tests/test_lca.py

Synthetic tree (real NCBI taxids, minimal Hominidae subtree):
    1 root
    └─ 9604 Hominidae (family)
       ├─ 9605 Homo (genus) ── 9606 Homo sapiens (species)
       └─ 9593 Gorilla (genus) ── 9595 Gorilla gorilla (species)
"""
import importlib.util
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pba", REPO / "bin" / "parse_blast_assignments.py")
pba = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pba)

NODES = [
    ("1", "1", "no rank"),
    ("9604", "1", "family"),
    ("9605", "9604", "genus"),
    ("9606", "9605", "species"),
    ("9593", "9604", "genus"),
    ("9595", "9593", "species"),
]
NAMES = {
    "1": "root", "9604": "Hominidae", "9605": "Homo",
    "9606": "Homo sapiens", "9593": "Gorilla", "9595": "Gorilla gorilla",
}


def main():
    failures = []
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "nodes.dmp").write_text("".join(f"{t}\t|\t{p}\t|\t{r}\t|\n" for t, p, r in NODES))
        (d / "names.dmp").write_text("".join(f"{t}\t|\t{n}\t|\t\t|\tscientific name\t|\n" for t, n in NAMES.items()))

        tax = pba.load_taxdump(str(d))
        if tax is None:
            failures.append("load_taxdump returned None for a valid fixture")
        else:
            checks = [
                (["9606", "9595"], "9604", "LCA(Homo sapiens, Gorilla gorilla) should escalate to Hominidae"),
                (["9606", "9605"], "9605", "LCA(Homo sapiens, Homo) should be Homo"),
                (["9606"], "9606", "LCA of a single taxon is itself"),
                (["9606", "9606"], "9606", "LCA of duplicates is the taxon"),
            ]
            for taxids, want, msg in checks:
                got = pba.lca_taxid(taxids, tax)
                if got != want:
                    failures.append(f"{msg}: got {got!r}, want {want!r}")
            if tax["names"].get("9604") != "Hominidae":
                failures.append("names.dmp not parsed (expected 9604 -> Hominidae)")

    if failures:
        print("LCA tests FAILED:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        sys.exit(1)
    print("LCA tests PASSED (synthetic taxdump: escalation, exact, single, duplicate)")


if __name__ == "__main__":
    main()
