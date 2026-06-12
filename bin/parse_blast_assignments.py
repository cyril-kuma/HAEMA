#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_counts(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_host(sseqid, stitle):
    token = (sseqid or "").split("|")[0].strip()
    if not token:
        token = (stitle or "").split()[0] if stitle else ""
    host = token.replace("_", " ").strip()
    parts = host.split()
    genus = parts[0] if parts else ""
    return host or "unresolved", genus


def read_reference_taxonomy(path):
    mapping = {}
    if not path or not Path(path).exists() or Path(path).stat().st_size == 0:
        return mapping
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seqid = row.get("seqid") or ""
            if seqid:
                mapping[seqid] = row
    return mapping


def read_blast(path, reference_taxonomy=None):
    hits = {}
    reference_taxonomy = reference_taxonomy or {}
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return hits
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 9:
                continue
            qseqid, sseqid, pident, length, qlen, slen, evalue, bitscore, stitle = row[:9]
            staxids = row[9] if len(row) > 9 else ""
            ref_meta = reference_taxonomy.get(sseqid, {})
            if (not split_taxids(staxids)) and ref_meta.get("taxid"):
                staxids = ref_meta.get("taxid")
            if (not stitle or stitle == "N/A") and ref_meta.get("scientific_name"):
                stitle = ref_meta.get("scientific_name")
            try:
                qlen_f = float(qlen)
                length_f = float(length)
                hit = {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "pident": float(pident),
                    "length": length_f,
                    "qlen": qlen_f,
                    "slen": float(slen),
                    "coverage": 100.0 * length_f / qlen_f if qlen_f else 0.0,
                    "evalue": float(evalue),
                    "bitscore": float(bitscore),
                    "stitle": stitle,
                    "staxids": staxids,
                    "tax_rank": ref_meta.get("rank", ""),
                }
            except ValueError:
                continue
            hits.setdefault(qseqid, []).append(hit)
    return hits


def split_taxids(value):
    taxids = []
    for item in str(value or "").replace(";", ",").split(","):
        item = item.strip()
        if item.isdigit() and item != "0":
            taxids.append(item)
    return taxids


def parse_taxdump_line(line):
    return [part.strip() for part in line.split("|")]


def load_taxdump(taxdump_dir):
    taxdump = Path(taxdump_dir) if taxdump_dir else None
    if not taxdump:
        return None
    nodes_path = taxdump / "nodes.dmp"
    names_path = taxdump / "names.dmp"
    if not nodes_path.exists() or not names_path.exists():
        return None

    parent = {}
    rank = {}
    with nodes_path.open() as handle:
        for line in handle:
            parts = parse_taxdump_line(line)
            if len(parts) >= 3:
                taxid, parent_taxid, tax_rank = parts[:3]
                parent[taxid] = parent_taxid
                rank[taxid] = tax_rank

    names = {}
    with names_path.open() as handle:
        for line in handle:
            parts = parse_taxdump_line(line)
            if len(parts) >= 4 and parts[3] == "scientific name":
                names[parts[0]] = parts[1]

    return {"parent": parent, "rank": rank, "names": names}


def lineage(taxid, taxonomy):
    parent = taxonomy["parent"]
    if taxid not in parent:
        return []
    path = [taxid]
    seen = {taxid}
    current = taxid
    while current in parent:
        next_taxid = parent[current]
        if next_taxid == current or next_taxid in seen:
            break
        path.append(next_taxid)
        seen.add(next_taxid)
        current = next_taxid
    return list(reversed(path))


def lca_taxid(taxids, taxonomy):
    lineages = [lineage(taxid, taxonomy) for taxid in taxids]
    lineages = [item for item in lineages if item]
    if not lineages:
        return None
    lca = None
    for column in zip(*lineages):
        if len(set(column)) == 1:
            lca = column[0]
        else:
            break
    return lca


def lca_assignment(top_hits, taxonomy):
    taxids = []
    for hit in top_hits:
        taxids.extend(split_taxids(hit.get("staxids", "")))
    taxids = sorted(set(taxids))
    if not taxids:
        return None
    taxid = lca_taxid(taxids, taxonomy)
    if not taxid:
        return None
    name = taxonomy["names"].get(taxid, f"taxid:{taxid}")
    rank = taxonomy["rank"].get(taxid, "no_rank")
    return name, rank, taxid, len(taxids)


def assign(hits, min_identity, min_coverage, top_bitscore_delta, blast_source, assignment_method, taxonomy):
    passing = [
        hit
        for hit in hits
        if hit["pident"] >= min_identity and hit["coverage"] >= min_coverage
    ]
    if not passing:
        return {
            "host_assignment": "unassigned",
            "taxon_rank": "none",
            "assignment_status": "no_confident_blast_hit",
            "confidence": "low",
            "n_top_hits": 0,
            "top_sseqid": "",
            "top_stitle": "",
            "pident": "",
            "coverage": "",
            "evalue": "",
            "bitscore": "",
            "blast_source": blast_source,
            "fallback_used": "false",
            "primary_assignment_status": "",
            "assignment_method": "blast_top_hit_conservative",
            "lca_taxid": "",
            "top_staxids": "",
        }
    passing.sort(key=lambda h: (-h["bitscore"], h["evalue"], -h["pident"]))
    best = passing[0]
    top = [hit for hit in passing if best["bitscore"] - hit["bitscore"] <= top_bitscore_delta]
    top_taxids = sorted(set(t for hit in top for t in split_taxids(hit.get("staxids", ""))))
    if assignment_method in {"taxid_lca", "conservative_lca"}:
        if taxonomy:
            lca = lca_assignment(top, taxonomy)
            if lca:
                name, rank, taxid, n_taxids = lca
                confidence = "high" if rank in {"species", "subspecies"} and n_taxids == 1 else "medium"
                return {
                    "host_assignment": name,
                    "taxon_rank": rank,
                    "assignment_status": "assigned_taxid_lca",
                    "confidence": confidence,
                    "n_top_hits": len(top),
                    "top_sseqid": best["sseqid"],
                    "top_stitle": best["stitle"],
                    "pident": f"{best['pident']:.3f}",
                    "coverage": f"{best['coverage']:.3f}",
                    "evalue": f"{best['evalue']:.3g}",
                    "bitscore": f"{best['bitscore']:.3f}",
                    "blast_source": blast_source,
                    "fallback_used": "false",
                    "primary_assignment_status": "",
                    "assignment_method": "blast_taxid_lca",
                    "lca_taxid": taxid,
                    "top_staxids": ";".join(top_taxids),
                }
        if len(top_taxids) == 1:
            host, _ = parse_host(best["sseqid"], best["stitle"])
            return {
                "host_assignment": host,
                "taxon_rank": best.get("tax_rank") or "species_or_reference_label",
                "assignment_status": "assigned_taxid_exact_lca",
                "confidence": "high",
                "n_top_hits": len(top),
                "top_sseqid": best["sseqid"],
                "top_stitle": best["stitle"],
                "pident": f"{best['pident']:.3f}",
                "coverage": f"{best['coverage']:.3f}",
                "evalue": f"{best['evalue']:.3g}",
                "bitscore": f"{best['bitscore']:.3f}",
                "blast_source": blast_source,
                "fallback_used": "false",
                "primary_assignment_status": "",
                "assignment_method": "blast_taxid_exact_lca",
                "lca_taxid": top_taxids[0],
                "top_staxids": top_taxids[0],
            }

    hosts = []
    genera = []
    for hit in top:
        host, genus = parse_host(hit["sseqid"], hit["stitle"])
        hosts.append(host)
        genera.append(genus)
    unique_hosts = sorted(set(hosts))
    unique_genera = sorted(set(genera))
    if len(unique_hosts) == 1:
        assignment = unique_hosts[0]
        rank = "species_or_reference_label"
        status = "assigned"
        confidence = "high"
    elif len(unique_genera) == 1 and unique_genera[0]:
        assignment = f"{unique_genera[0]} sp."
        rank = "genus"
        status = "ambiguous_species_lca_genus"
        confidence = "medium"
    else:
        assignment = "ambiguous"
        rank = "unresolved"
        status = "ambiguous_top_hits"
        confidence = "low"
    return {
        "host_assignment": assignment,
        "taxon_rank": rank,
        "assignment_status": status,
        "confidence": confidence,
        "n_top_hits": len(top),
        "top_sseqid": best["sseqid"],
        "top_stitle": best["stitle"],
        "pident": f"{best['pident']:.3f}",
        "coverage": f"{best['coverage']:.3f}",
        "evalue": f"{best['evalue']:.3g}",
        "bitscore": f"{best['bitscore']:.3f}",
        "blast_source": blast_source,
        "fallback_used": "false",
        "primary_assignment_status": "",
        "assignment_method": "blast_top_hit_conservative",
        "lca_taxid": "",
        "top_staxids": ";".join(top_taxids),
    }


def main():
    parser = argparse.ArgumentParser(description="Parse BLAST output into conservative host assignments")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--blast", required=True)
    parser.add_argument("--blast-source", default="curated_panel")
    parser.add_argument("--fallback-blast", default="")
    parser.add_argument("--fallback-source", default="nt")
    parser.add_argument("--assignment-method", default="conservative_lca")
    parser.add_argument("--taxdump-dir", default="")
    parser.add_argument("--reference-taxonomy", default="")
    parser.add_argument("--min-identity", type=float, required=True)
    parser.add_argument("--min-coverage", type=float, required=True)
    parser.add_argument("--top-bitscore-delta", type=float, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    counts = read_counts(args.counts)
    reference_taxonomy = read_reference_taxonomy(args.reference_taxonomy)
    primary_hits = read_blast(args.blast, reference_taxonomy)
    fallback_hits = read_blast(args.fallback_blast) if args.fallback_blast else {}
    taxonomy = load_taxdump(args.taxdump_dir)
    fieldnames = [
        "sample_uid",
        "run_id",
        "sample_id",
        "barcode_id",
        "marker",
        "cluster_id",
        "asv_id",
        "sequence",
        "count",
        "fraction",
        "retained",
        "host_assignment",
        "taxon_rank",
        "assignment_status",
        "confidence",
        "n_top_hits",
        "top_sseqid",
        "top_stitle",
        "pident",
        "coverage",
        "evalue",
        "bitscore",
        "blast_source",
        "fallback_used",
        "primary_assignment_status",
        "assignment_method",
        "lca_taxid",
        "top_staxids",
    ]
    with Path(args.output).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in counts:
            if row.get("retained") != "true":
                continue
            result = assign(
                primary_hits.get(row["asv_id"], []),
                args.min_identity,
                args.min_coverage,
                args.top_bitscore_delta,
                args.blast_source,
                args.assignment_method,
                taxonomy,
            )
            if args.fallback_blast and result["assignment_status"] == "no_confident_blast_hit":
                primary_status = result["assignment_status"]
                fallback_result = assign(
                    fallback_hits.get(row["asv_id"], []),
                    args.min_identity,
                    args.min_coverage,
                    args.top_bitscore_delta,
                    args.fallback_source,
                    args.assignment_method,
                    taxonomy,
                )
                if fallback_result["assignment_status"] != "no_confident_blast_hit":
                    result = fallback_result
                    result["fallback_used"] = "true"
                    result["primary_assignment_status"] = primary_status
                else:
                    result["fallback_used"] = "attempted_no_confident_hit"
                    result["primary_assignment_status"] = primary_status
            out = dict(row)
            out.update(result)
            writer.writerow(out)


if __name__ == "__main__":
    main()
