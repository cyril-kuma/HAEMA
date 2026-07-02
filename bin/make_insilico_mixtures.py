#!/usr/bin/env python3
"""Build synthetic known-ratio mixed-host controls for in-silico denoising calibration.

The "bridge" for when wet-lab known-ratio mixed controls are not available: pool reads from
*pure* single-host sources at defined ratios so the true minor-host read fraction is known, then
run them through classification (`classify_reads.py`) to measure detection limit and calibrate
thresholds — without any wet-lab work. Adapts Logue et al. (2016) classify-then-count.

Two read sources, same output contract:

  * `--pool NAME=reads.fastq[.gz]` (repeatable) — REAL pure single-host read pools (preferred once
    the full run yields cattle/goat/sheep-only samples).
  * `--simulate` with `--reference panel.fasta --primers primers.csv --hosts A,B,...` — simulate
    ONT-like reads by in-silico PCR of each host's marker amplicon(s) + an ONT error model. Fully
    reproducible; use when real non-target-host pools are not yet available.

Honest scope: simulated/pooled reads capture *sequencing + bioinformatic* behaviour (detection
floor, threshold choice) but NOT extraction/PCR/mtDNA-copy/digestion bias — so this calibrates the
denoising thresholds, it does NOT license quantitative blood-fraction claims
(`host_fractions_benchmarked` stays false). See docs/denoising_redesign_plan.md.
"""
import argparse
import gzip
import random
import re
from pathlib import Path

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMP = str.maketrans("ACGTNRYSWKMBDHV", "TGCANYRSWMKVHDB")


def revcomp(s):
    return s.translate(COMP)[::-1]


def iupac_regex(primer):
    return re.compile("".join(f"[{IUPAC.get(b, 'ACGT')}]" for b in primer.upper()))


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            else:
                seq.append(line.strip().upper())
        if name:
            yield name, "".join(seq)


def _open_out(path):
    return gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w")


def find_primer(seq, primer, max_mismatch=3):
    """Return start index of the best (fewest-mismatch) match of an IUPAC primer, else None."""
    plen = len(primer)
    allowed = [set(IUPAC.get(b, "ACGT")) for b in primer.upper()]
    best_i, best_mm = None, max_mismatch + 1
    for i in range(0, len(seq) - plen + 1):
        mm = 0
        for j in range(plen):
            if seq[i + j] not in allowed[j]:
                mm += 1
                if mm >= best_mm:
                    break
        if mm < best_mm:
            best_mm, best_i = mm, i
            if mm == 0:
                break
    return best_i


def insilico_pcr(ref_seq, fwd, rev, max_mismatch=3):
    """Extract the amplicon between fwd (sense) and revcomp(rev) (antisense) on ref_seq."""
    f = find_primer(ref_seq, fwd, max_mismatch)
    rc = revcomp(rev)
    if f is None:
        return None
    r = find_primer(ref_seq[f + len(fwd):], rc, max_mismatch)
    if r is None:
        return None
    end = f + len(fwd) + r + len(rc)
    amp = ref_seq[f:end]
    return amp if 80 <= len(amp) <= 2000 else None


def simulate_read(amplicon, rng, sub=0.015, ins=0.008, dele=0.010, min_frac=0.6):
    """Emit one ONT-like read: random sub-span of the amplicon with sub/indel errors."""
    L = len(amplicon)
    start = rng.randint(0, int(L * (1 - min_frac)))
    end = rng.randint(start + int(L * min_frac), L)
    span = amplicon[start:end]
    out = []
    for base in span:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(rng.choice("ACGT"))
            out.append(base)
        elif r < dele + ins + sub:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def write_fastq(records, path):
    with _open_out(path) as fh:
        for rid, seq in records:
            fh.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")


def load_pool_fastq(path):
    reads = []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip()
            fh.readline()
            fh.readline()
            if s:
                reads.append(s)
    return reads


def build_pools_simulated(reference, primers_csv, hosts, reads_per_pool, rng, args):
    """Return {host: [read_seq,...]} simulated across all markers in primers.csv."""
    refs = {name.split("|")[0]: seq for name, seq in read_fasta(reference)}
    markers = []
    with open(primers_csv) as fh:
        header = fh.readline().rstrip().split(",")
        for line in fh:
            parts = dict(zip(header, line.rstrip().split(",")))
            markers.append((parts["Gene"], parts["Forward_Primer"], parts["Reverse_Primer"]))
    pools = {}
    for host in hosts:
        if host not in refs:
            raise SystemExit(f"host '{host}' not found in reference {reference}")
        amps = []
        for gene, fwd, rev in markers:
            amp = insilico_pcr(refs[host], fwd, rev)
            if amp:
                amps.append((gene, amp))
        if not amps:
            raise SystemExit(f"no amplicon recovered by in-silico PCR for {host}")
        reads = []
        per = max(1, reads_per_pool // len(amps))
        for gene, amp in amps:
            for k in range(per):
                reads.append(simulate_read(amp, rng, args.sub, args.ins, args.dele))
        pools[host] = reads
        print(f"[simulate] {host}: {len(reads)} reads across {len(amps)} markers "
              f"(amplens={[len(a) for _, a in amps]})")
    return pools


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pool", action="append", default=[], metavar="NAME=reads.fastq",
                    help="Real pure single-host read pool (repeatable).")
    ap.add_argument("--simulate", action="store_true", help="Simulate pools from a reference FASTA.")
    ap.add_argument("--reference", help="Reference FASTA (deflines Genus_species|...) for --simulate.")
    ap.add_argument("--primers", help="primers.csv (Gene,Forward_Primer,Reverse_Primer,Size).")
    ap.add_argument("--hosts", help="Comma-separated host names (must match reference) for --simulate.")
    ap.add_argument("--reads-per-pool", type=int, default=2000)
    ap.add_argument("--sub", type=float, default=0.015)
    ap.add_argument("--ins", type=float, default=0.008)
    ap.add_argument("--dele", type=float, default=0.010)
    ap.add_argument("--ratios", default="99:1,95:5,90:10,75:25,50:50",
                    help="Major:minor read ratios; each is built both directions.")
    ap.add_argument("--depth", type=int, default=2000, help="Total reads per synthetic mixture.")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    if args.simulate:
        hosts = [h.strip() for h in args.hosts.split(",") if h.strip()]
        pools = build_pools_simulated(args.reference, args.primers, hosts, args.reads_per_pool, rng, args)
    else:
        pools = {}
        for spec in args.pool:
            name, _, path = spec.partition("=")
            pools[name] = load_pool_fastq(path)
            print(f"[pool] {name}: {len(pools[name])} reads from {path}")
    hosts = list(pools)
    if len(hosts) < 2:
        raise SystemExit("need >=2 host pools to build mixtures")

    ratios = []
    for r in args.ratios.split(","):
        a, b = r.split(":")
        ratios.append((int(a), int(b)))

    manifest = out / "mixture_manifest.tsv"
    with manifest.open("w") as mh:
        mh.write("mixture_id\tmajor_host\tminor_host\tmajor_ratio\tminor_ratio\t"
                 "true_minor_fraction\tdepth\treplicate\n")
        # every ordered host pair, every ratio, every replicate
        for i, major in enumerate(hosts):
            for minor in hosts:
                if major == minor:
                    continue
                for maj, mino in ratios:
                    for rep in range(args.replicates):
                        n_major = round(args.depth * maj / (maj + mino))
                        n_minor = args.depth - n_major
                        reads = []
                        for src, n, tag in ((major, n_major, "MAJ"), (minor, n_minor, "MIN")):
                            pool = pools[src]
                            pick = [rng.choice(pool) for _ in range(n)]  # with replacement
                            reads += [(f"{tag}_{src}_{k}", s) for k, s in enumerate(pick)]
                        rng.shuffle(reads)
                        mid = f"{major}__{minor}__{maj}-{mino}__r{rep}"
                        write_fastq(reads, out / f"{mid}.fastq")
                        true_minor = n_minor / (n_major + n_minor)
                        mh.write(f"{mid}\t{major}\t{minor}\t{maj}\t{mino}\t"
                                 f"{true_minor:.4f}\t{args.depth}\t{rep}\n")
    n_mix = sum(1 for _ in manifest.open()) - 1
    print(f"[mixtures] wrote {n_mix} synthetic mixtures + {manifest}")


if __name__ == "__main__":
    main()
