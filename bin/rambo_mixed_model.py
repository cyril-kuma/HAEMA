#!/usr/bin/env python3
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def is_assigned(row):
    return row.get("host_assignment") not in {"", "unassigned", "ambiguous"}


def norm_taxon(name):
    return " ".join((name or "").strip().lower().replace("_", " ").split())


def taxon_match(expected, observed):
    """Return 'species', 'genus', or '' for how well an observed host matches the expected host."""
    e, o = norm_taxon(expected), norm_taxon(observed)
    if not e or not o:
        return ""
    if e == o:
        return "species"
    if e.split() and o.split() and e.split()[0] == o.split()[0]:
        return "genus"
    return ""


def split_expected_hosts(value):
    """Parse a samplesheet expected-host field into a set of declared host names.

    Supports a single host or a multi-host lab-prepared mixture declared as a
    ';'-, '|'- or '+'-separated list, e.g. "Homo sapiens; Capra hircus".
    """
    if not value:
        return []
    raw = value.replace("|", ";").replace("+", ";").split(";")
    return [h.strip() for h in raw if h.strip()]


def recovers(expected_host, observed_hosts):
    """Best match level ('species'/'genus'/'') of an expected host against observed calls."""
    levels = [taxon_match(expected_host, o) for o in observed_hosts]
    if "species" in levels:
        return "species"
    if "genus" in levels:
        return "genus"
    return ""


def check_positive_controls(rows, host_call_rows):
    """Compare control host calls to the host(s) declared in the samplesheet.

    Handles BOTH single-host positive controls and lab-prepared MIXED-host controls (multiple
    expected hosts). Any sample with a non-empty `expected_host_scientific_name` is evaluated, so
    lab-prepared mixed-feeding samples (sample_type=sample) are checked once their known
    composition is declared. This is a control-integrity / recovery check, NOT a calibrated
    benchmark of host fractions. Returns (check_rows, counts).
    """
    expected_by_sample = {}
    control_status_by_sample = {}
    for row in rows:
        uid = row.get("sample_uid")
        if not uid:
            continue
        control_status_by_sample[uid] = row.get("control_status", "") or row.get("sample_type", "")
        exp = (row.get("expected_host_scientific_name") or "").strip()
        if exp:
            expected_by_sample[uid] = exp

    observed_by_sample = defaultdict(list)
    markers_by_sample = defaultdict(set)
    for hc in host_call_rows:
        host = hc.get("host_assignment") or ""
        if host and host not in {"unassigned", "ambiguous"}:
            observed_by_sample[hc.get("sample_uid")].append(host)
            markers_by_sample[hc.get("sample_uid")].add(hc.get("marker"))

    # Evaluate any control with a declared expected host, plus positive controls without one.
    candidates = sorted(set(expected_by_sample) | {
        uid for uid, st in control_status_by_sample.items() if st == "positive_control"
    })
    check_rows = []
    counts = Counter()
    for uid in candidates:
        expected = split_expected_hosts(expected_by_sample.get(uid, ""))
        observed_hosts = sorted(set(observed_by_sample.get(uid, [])))
        n_expected = len(expected)
        control_kind = "mixed_host_control" if n_expected > 1 else "single_host_control"

        recovered = [e for e in expected if recovers(e, observed_hosts)]
        missing = [e for e in expected if not recovers(e, observed_hosts)]
        unexpected = [o for o in observed_hosts if not any(taxon_match(e, o) for e in expected)] if expected else []

        if not expected:
            status, note = "indeterminate_no_expected_host", "no expected_host_scientific_name declared for this control"
        elif not observed_hosts:
            status, note = "fail_no_host_signal", "control produced no retained host calls"
        elif n_expected == 1:
            lvl = recovers(expected[0], observed_hosts)
            if lvl == "species":
                status, note = "pass", "expected host recovered (species-level match)"
            elif lvl == "genus":
                status, note = "pass_genus", "expected host recovered at genus level only"
            else:
                status, note = "fail_unexpected_host", "expected host not among observed host calls"
        else:  # mixed-host control: did we recover all declared hosts?
            if len(recovered) == n_expected and not unexpected:
                status, note = "mixed_pass_all", f"all {n_expected} expected hosts recovered"
            elif len(recovered) == n_expected:
                status, note = "mixed_pass_with_extra", f"all {n_expected} expected recovered, plus unexpected host(s)"
            elif recovered:
                status, note = "mixed_partial", f"recovered {len(recovered)}/{n_expected} expected hosts"
            else:
                status, note = "mixed_fail", "no expected hosts recovered"
        counts[status] += 1
        check_rows.append(
            {
                "sample_uid": uid,
                "control_kind": control_kind,
                "control_status": control_status_by_sample.get(uid, ""),
                "expected_hosts": ";".join(expected) if expected else "(none provided)",
                "observed_hosts": ";".join(observed_hosts) if observed_hosts else "(none)",
                "n_expected": n_expected,
                "n_recovered": len(recovered),
                "missing_hosts": ";".join(missing) if missing else "",
                "unexpected_hosts": ";".join(unexpected) if unexpected else "",
                "markers_with_signal": ";".join(sorted(markers_by_sample.get(uid, []))) or "(none)",
                "status": status,
                "note": note,
            }
        )
    return check_rows, counts


def main():
    parser = argparse.ArgumentParser(description="RAMBO-style mixed host evidence modelling from endpoint feature calls")
    parser.add_argument("--master-endpoint", required=True)
    parser.add_argument("--marker-summary", required=True)
    parser.add_argument("--min-host-reads", type=int, default=3)
    parser.add_argument("--min-host-fraction", type=float, default=0.02)
    parser.add_argument("--include-contaminants", default="false")
    parser.add_argument("--output-evidence", required=True)
    parser.add_argument("--output-host-calls", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-control-check", required=True)
    args = parser.parse_args()

    include_contaminants = str(args.include_contaminants).strip().lower() in {"1", "true", "yes", "y", "on"}
    rows = read_tsv(args.master_endpoint)
    marker_summary = {
        (row.get("sample_uid"), row.get("marker")): row
        for row in read_tsv(args.marker_summary)
    }

    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.get("sample_uid"), row.get("marker"))].append(row)

    evidence_rows = []
    host_call_rows = []
    summary_counter = Counter()
    for (sample_uid, marker), feature_rows in sorted(grouped.items()):
        total_reads = sum(as_int(row.get("count")) for row in feature_rows)
        host_counts = Counter()
        host_features = defaultdict(list)
        for row in feature_rows:
            count = as_int(row.get("count"))
            host = row.get("host_assignment") or "unassigned"
            contaminant = row.get("contamination_flag") == "true"
            retained_for_model = (
                is_assigned(row)
                and (include_contaminants or not contaminant)
            )
            if retained_for_model:
                host_counts[host] += count
                host_features[host].append(row)

            evidence_rows.append(
                {
                    "sample_uid": sample_uid,
                    "run_id": row.get("run_id", ""),
                    "sample_id": row.get("sample_id", ""),
                    "barcode_id": row.get("barcode_id", ""),
                    "control_status": row.get("control_status", ""),
                    "marker": marker,
                    "cluster_id": row.get("cluster_id", ""),
                    "asv_id": row.get("asv_id", ""),
                    "host_assignment": host,
                    "taxon_rank": row.get("taxon_rank", ""),
                    "count": count,
                    "feature_fraction": row.get("fraction", ""),
                    "sample_marker_fraction": f"{(count / total_reads) if total_reads else 0.0:.6f}",
                    "confidence": row.get("confidence", ""),
                    "assignment_status": row.get("assignment_status", ""),
                    "blast_source": row.get("blast_source", ""),
                    "contamination_flag": row.get("contamination_flag", "false"),
                    "retained_for_mixed_model": str(retained_for_model).lower(),
                    "model_reason": "assigned_noncontaminant" if retained_for_model else "excluded_unassigned_or_contaminant",
                }
            )

        retained_hosts = []
        for host, count in host_counts.most_common():
            fraction = count / total_reads if total_reads else 0.0
            if count >= args.min_host_reads and fraction >= args.min_host_fraction:
                retained_hosts.append((host, count, fraction))

        status = "no_host_signal"
        if len(retained_hosts) == 1:
            status = "single_host"
        elif len(retained_hosts) > 1:
            status = "mixed_host"
        summary_counter[status] += 1

        marker_meta = marker_summary.get((sample_uid, marker), {})
        for rank, (host, count, fraction) in enumerate(retained_hosts, start=1):
            examples = host_features.get(host, [])
            best = sorted(
                examples,
                key=lambda row: (-as_float(row.get("bitscore")), -as_float(row.get("pident")), row.get("asv_id", "")),
            )[0] if examples else {}
            host_call_rows.append(
                {
                    "sample_uid": sample_uid,
                    "run_id": best.get("run_id", ""),
                    "sample_id": best.get("sample_id", ""),
                    "barcode_id": best.get("barcode_id", ""),
                    "control_status": best.get("control_status", ""),
                    "marker": marker,
                    "best_cluster_id": best.get("cluster_id", ""),
                    "host_assignment": host,
                    "host_rank": rank,
                    "host_reads": count,
                    "host_fraction": f"{fraction:.6f}",
                    "n_supporting_features": len(examples),
                    "best_feature_id": best.get("asv_id", ""),
                    "best_confidence": best.get("confidence", ""),
                    "best_assignment_status": best.get("assignment_status", ""),
                    "mixed_status": status,
                    "total_marker_reads": total_reads,
                    "mixed_template_warning": marker_meta.get("mixed_template_warning", "false"),
                    "model": "rambo_style_abundance_evidence",
                }
            )
        if not retained_hosts:
            host_call_rows.append(
                {
                    "sample_uid": sample_uid,
                    "run_id": feature_rows[0].get("run_id", "") if feature_rows else "",
                    "sample_id": feature_rows[0].get("sample_id", "") if feature_rows else "",
                    "barcode_id": feature_rows[0].get("barcode_id", "") if feature_rows else "",
                    "control_status": feature_rows[0].get("control_status", "") if feature_rows else "",
                    "marker": marker,
                    "best_cluster_id": "",
                    "host_assignment": "unassigned",
                    "host_rank": "",
                    "host_reads": 0,
                    "host_fraction": "0.000000",
                    "n_supporting_features": 0,
                    "best_feature_id": "",
                    "best_confidence": "",
                    "best_assignment_status": "no_model_retained_host_signal",
                    "mixed_status": status,
                    "total_marker_reads": total_reads,
                    "mixed_template_warning": marker_meta.get("mixed_template_warning", "false"),
                    "model": "rambo_style_abundance_evidence",
                }
            )

    write_tsv(
        args.output_evidence,
        evidence_rows,
        [
            "sample_uid", "run_id", "sample_id", "barcode_id", "control_status", "marker", "cluster_id", "asv_id",
            "host_assignment", "taxon_rank", "count", "feature_fraction", "sample_marker_fraction",
            "confidence", "assignment_status", "blast_source", "contamination_flag",
            "retained_for_mixed_model", "model_reason",
        ],
    )
    write_tsv(
        args.output_host_calls,
        host_call_rows,
        [
            "sample_uid", "run_id", "sample_id", "barcode_id", "control_status", "marker", "best_cluster_id", "host_assignment",
            "host_rank", "host_reads", "host_fraction", "n_supporting_features", "best_feature_id",
            "best_confidence", "best_assignment_status", "mixed_status", "total_marker_reads",
            "mixed_template_warning", "model",
        ],
    )
    control_check_rows, control_counts = check_positive_controls(rows, host_call_rows)
    write_tsv(
        args.output_control_check,
        control_check_rows,
        ["sample_uid", "control_kind", "control_status", "expected_hosts", "observed_hosts",
         "n_expected", "n_recovered", "missing_hosts", "unexpected_hosts", "markers_with_signal",
         "status", "note"],
    )

    # Control / benchmark roll-up derived from the check rows.
    single = [r for r in control_check_rows if r["control_kind"] == "single_host_control" and r["expected_hosts"] != "(none provided)"]
    mixed = [r for r in control_check_rows if r["control_kind"] == "mixed_host_control"]
    single_pass = sum(1 for r in single if r["status"] in {"pass", "pass_genus"})
    mixed_expected = sum(int(r["n_expected"]) for r in mixed)
    mixed_recovered = sum(int(r["n_recovered"]) for r in mixed)
    mixed_all_pass = sum(1 for r in mixed if r["status"] in {"mixed_pass_all", "mixed_pass_with_extra"})

    summary_rows = [
        {"metric": "sample_marker_groups", "value": sum(summary_counter.values())},
        {"metric": "single_host_groups", "value": summary_counter["single_host"]},
        {"metric": "mixed_host_groups", "value": summary_counter["mixed_host"]},
        {"metric": "no_host_signal_groups", "value": summary_counter["no_host_signal"]},
        {"metric": "min_host_reads", "value": args.min_host_reads},
        {"metric": "min_host_fraction", "value": args.min_host_fraction},
        {"metric": "include_contaminants", "value": str(include_contaminants).lower()},
        {"metric": "single_host_controls_total", "value": len(single)},
        {"metric": "single_host_controls_passed", "value": single_pass},
        {"metric": "mixed_host_controls_total", "value": len(mixed)},
        {"metric": "mixed_host_controls_all_recovered", "value": mixed_all_pass},
        {"metric": "mixed_expected_hosts", "value": mixed_expected},
        {"metric": "mixed_expected_hosts_recovered", "value": mixed_recovered},
        {"metric": "mixed_host_recovery_rate", "value": f"{(mixed_recovered/mixed_expected):.3f}" if mixed_expected else "NA"},
        # Scientific safeguard: host fractions are abundance EVIDENCE, not validated quantitative
        # estimates. The recovery rate above is a sensitivity benchmark only when the declared
        # mixed-host composition comes from genuine lab-prepared controls.
        # Validation-status fields (kept explicit so reports/manifests cannot overclaim):
        #  - in-silico mixtures calibrate computational DETECTION only (docs/denoising_calibration);
        #  - quantitative read-fraction<->blood-proportion calibration needs wet-lab known-ratio
        #    controls (docs/mixed_host_control_protocol.md) and is NOT done.
        {"metric": "host_fractions_benchmarked", "value": "false"},
        {"metric": "wetlab_known_ratio_validation", "value": "false"},
    ]
    write_tsv(args.output_summary, summary_rows, ["metric", "value"])


if __name__ == "__main__":
    main()
