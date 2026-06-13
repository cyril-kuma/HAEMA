#!/usr/bin/env python3
import gzip
from pathlib import Path


def open_text(path, mode="rt"):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode)


def read_fastq(path):
    with open_text(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Truncated FASTQ record in {path}")
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def write_fastq_record(handle, header, seq, qual):
    handle.write(f"{header}\n{seq}\n+\n{qual}\n")


def mean_q(qual):
    if not qual:
        return 0.0
    return sum(max(0, ord(ch) - 33) for ch in qual) / len(qual)


IUPAC = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T", "N"},
}

COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1].upper()


def primer_mismatches(query, primer):
    mismatches = 0
    query = query.upper()
    primer = primer.upper()
    for base, code in zip(query, primer):
        allowed = IUPAC.get(code, {code})
        if base not in allowed:
            mismatches += 1
    return mismatches


def best_primer_match(seq, primer, region, max_error_rate):
    seq = seq.upper()
    primer = primer.upper()
    plen = len(primer)
    if plen == 0 or len(seq) < plen:
        return None
    max_mismatches = int(plen * max_error_rate + 1e-9)
    # region is always a (start, stop) window; callers (find_terminal_pair) pass terminal windows.
    if isinstance(region, tuple):
        start, stop = region
        start = max(0, start)
        stop = min(len(seq) - plen, stop)
        stop = max(start, stop) + 1
    else:
        raise ValueError(f"Unsupported match region: {region}")

    best = None
    for pos in range(start, stop):
        mismatches = primer_mismatches(seq[pos : pos + plen], primer)
        if mismatches <= max_mismatches:
            candidate = {
                "pos": pos,
                "end": pos + plen,
                "mismatches": mismatches,
                "error_rate": mismatches / plen,
                "primer": primer,
            }
            if best is None or (candidate["mismatches"], candidate["pos"]) < (best["mismatches"], best["pos"]):
                best = candidate
    return best


def find_terminal_pair(seq, forward_primer, reverse_primer, window, max_error_rate):
    seq = seq.upper()
    fwd = forward_primer.upper()
    rev = reverse_primer.upper()
    rc_rev = revcomp(rev)
    rc_fwd = revcomp(fwd)
    n = len(seq)
    left_region = (0, min(window, n))
    right_region = (max(0, n - window - max(len(rc_rev), len(rc_fwd))), n)

    fwd_start = best_primer_match(seq, fwd, left_region, max_error_rate)
    rev_end = best_primer_match(seq, rc_rev, right_region, max_error_rate)
    candidates = []
    if fwd_start and rev_end and fwd_start["end"] <= rev_end["pos"]:
        candidates.append(
            {
                "orientation": "forward",
                "start": fwd_start["end"],
                "end": rev_end["pos"],
                "mismatches": fwd_start["mismatches"] + rev_end["mismatches"],
                "error_rate": (fwd_start["error_rate"] + rev_end["error_rate"]) / 2,
            }
        )

    rev_start = best_primer_match(seq, rev, left_region, max_error_rate)
    fwd_end = best_primer_match(seq, rc_fwd, right_region, max_error_rate)
    if rev_start and fwd_end and rev_start["end"] <= fwd_end["pos"]:
        candidates.append(
            {
                "orientation": "reverse",
                "start": rev_start["end"],
                "end": fwd_end["pos"],
                "mismatches": rev_start["mismatches"] + fwd_end["mismatches"],
                "error_rate": (rev_start["error_rate"] + fwd_end["error_rate"]) / 2,
            }
        )

    if not candidates:
        return None
    return sorted(candidates, key=lambda c: (c["mismatches"], -1 * (c["end"] - c["start"])))[0]
