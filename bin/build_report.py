#!/usr/bin/env python3
import argparse
import csv
import html
from collections import Counter
from pathlib import Path


def read_tsv(path):
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return []
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def positive_control_banner(pc_rows):
    """Render a pass/fail banner for control integrity (single- and mixed-host controls)."""
    if not pc_rows:
        return ('<div class="note">No control check available (no declared controls, or '
                '--enable_rambo_model was off). Set <code>expected_host_scientific_name</code> on a '
                'positive control (single host) or a lab-prepared mixed sample '
                '(semicolon-separated hosts) to enable automated control/recovery checking.</div>')
    passed = sum(1 for r in pc_rows if r.get("status", "").startswith(("pass", "mixed_pass")))
    failed = sum(1 for r in pc_rows if r.get("status", "").startswith(("fail", "mixed_fail")))
    partial = sum(1 for r in pc_rows if r.get("status", "") == "mixed_partial")
    indet = sum(1 for r in pc_rows if r.get("status", "").startswith("indeterminate"))
    cls = "ok" if failed == 0 and partial == 0 and indet == 0 else ("warn" if failed == 0 else "fail")
    mixed = [r for r in pc_rows if r.get("control_kind") == "mixed_host_control"]
    bench = ""
    if mixed:
        exp = sum(int(r.get("n_expected") or 0) for r in mixed)
        rec = sum(int(r.get("n_recovered") or 0) for r in mixed)
        bench = (f' Mixed-host controls: {rec}/{exp} declared hosts recovered '
                 f'(sensitivity {rec/exp:.0%}).' if exp else '')
    return (f'<div class="banner {cls}"><strong>Controls:</strong> {passed} passed, {partial} partial, '
            f'{failed} failed, {indet} indeterminate.{bench} A "pass" means the control recovered its '
            f'declared host(s); this is a recovery/integrity check, not a host-fraction benchmark.</div>')


def table(rows, limit=20):
    if not rows:
        return "<p>No rows.</p>"
    fields = list(rows[0].keys())
    body = ["<table><thead><tr>"]
    body.extend(f"<th>{html.escape(field)}</th>" for field in fields)
    body.append("</tr></thead><tbody>")
    for row in rows[:limit]:
        body.append("<tr>")
        body.extend(f"<td>{html.escape(str(row.get(field, '')))}</td>" for field in fields)
        body.append("</tr>")
    body.append("</tbody></table>")
    return "\n".join(body)


def main():
    parser = argparse.ArgumentParser(description="Build lightweight HÆMA HTML report")
    parser.add_argument("--master-endpoint", required=True)
    parser.add_argument("--sample-summary", required=True)
    parser.add_argument("--marker-summary", required=True)
    parser.add_argument("--contamination-flags", required=True)
    parser.add_argument("--control-check", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    master = read_tsv(args.master_endpoint)
    samples = read_tsv(args.sample_summary)
    markers = read_tsv(args.marker_summary)
    contamination = read_tsv(args.contamination_flags)
    pc_rows = read_tsv(args.control_check)

    host_counts = Counter(row.get("host_assignment", "unassigned") for row in master)
    marker_counts = Counter(row.get("marker", "") for row in master)
    control_counts = Counter(row.get("sample_type", row.get("control_status", "")) for row in master)

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>HÆMA Blood-Meal Pipeline Report</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 2rem; line-height: 1.4; color: #1f2933; }}
    h1, h2 {{ color: #102a43; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0 2rem; font-size: 0.9rem; }}
    th, td {{ border: 1px solid #d9e2ec; padding: 0.35rem 0.5rem; text-align: left; vertical-align: top; }}
    th {{ background: #f0f4f8; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; }}
    .metric {{ border: 1px solid #d9e2ec; padding: 1rem; border-radius: 6px; }}
    .metric strong {{ display: block; font-size: 1.6rem; color: #334e68; }}
    .banner {{ padding: 0.75rem 1rem; border-radius: 6px; margin: 0.5rem 0; }}
    .banner.ok {{ background: #e3f9e5; border: 1px solid #51cf66; }}
    .banner.warn {{ background: #fff9db; border: 1px solid #f0c000; }}
    .banner.fail {{ background: #ffe3e3; border: 1px solid #ff6b6b; }}
    .note {{ background: #f0f4f8; border-left: 4px solid #829ab1; padding: 0.5rem 1rem; margin: 0.5rem 0; font-size: 0.9rem; }}
    .caveat {{ background: #fff4e6; border-left: 4px solid #f08c00; padding: 0.5rem 1rem; margin: 0.5rem 0; font-size: 0.9rem; }}
  </style>
</head>
<body>
  <h1>HÆMA Blood-Meal Pipeline Report</h1>
  <div class="caveat"><strong>Scientific caveats — read before interpreting:</strong> Host
  fractions are <em>abundance evidence summaries, not validated quantitative estimates</em>.
  Mixed-host denoising thresholds are not yet benchmarked against known mixtures. Taxonomy is
  bounded by the curated reference panel; hosts absent from the panel may be unassigned or
  assigned to the nearest relative (supply an <code>nt</code> fallback to explore these). See
  <code>docs/limitations.md</code>.</div>
  <div class="grid">
    <div class="metric"><span>Endpoint rows</span><strong>{len(master)}</strong></div>
    <div class="metric"><span>Samples with retained features</span><strong>{len(samples)}</strong></div>
    <div class="metric"><span>Marker summaries</span><strong>{len(markers)}</strong></div>
    <div class="metric"><span>Contamination flags</span><strong>{len(contamination)}</strong></div>
  </div>
  <h2>Positive-Control Integrity</h2>
  {positive_control_banner(pc_rows)}
  {table(pc_rows)}
  <h2>Feature status</h2>
  <div class="note">
    <strong>Implemented &amp; tested:</strong> validation, trim/QC/split, UMAP/HDBSCAN denoising
    (with greedy fallback), curated BLAST + conservative LCA, contamination flags, mixed-host
    evidence, endpoints.<br>
    <strong>Implemented, not biologically validated:</strong> mixed-host thresholds, host-fraction
    quantitation, Medaka mixed-host polishing.<br>
    <strong>Available only with external resources:</strong> <code>nt</code> fallback BLAST,
    taxdump LCA, formal phyloseq/decontam (haema-r image).<br>
    <strong>Staged / not bundled:</strong> Barbell/Deepbinner demultiplexing, POD5/Dorado basecalling.
  </div>
  <h2>Host Assignments</h2>
  {table([{"host_assignment": k, "records": v} for k, v in host_counts.most_common()])}
  <h2>Marker Records</h2>
  {table([{"marker": k, "records": v} for k, v in marker_counts.most_common()])}
  <h2>Control Status</h2>
  {table([{"sample_type": k, "records": v} for k, v in control_counts.most_common()])}
  <h2>Sample Summary</h2>
  {table(samples)}
  <h2>Marker Summary</h2>
  {table(markers)}
  <h2>Contamination Flags</h2>
  {table(contamination)}
</body>
</html>
"""
    Path(args.output).write_text(html_text)


if __name__ == "__main__":
    main()
