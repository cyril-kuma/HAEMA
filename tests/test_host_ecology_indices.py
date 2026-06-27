#!/usr/bin/env python3
"""Unit tests for bin/compute_host_ecology_indices.py.

Run from the repo root:  python3 tests/test_host_ecology_indices.py
"""
import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "bin" / "compute_host_ecology_indices.py"
spec = importlib.util.spec_from_file_location("hosteco", SCRIPT)
hosteco = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hosteco)

SAMPLESHEET = """samplesheet_schema_version,run_id,minknow_run_folder,barcode_id,sample_id,specimen_id,sample_type,control_type,expected_host_scientific_name,expected_host_taxid,expected_marker_result,species,sibling_species,feeding_status,collection_date,collection_time,collection_location,collection_region,bioclimatic_zone,latitude,longitude,collection_context,collection_method,specimen_sex,extraction_batch,pcr_batch,library_batch,barcode_kit,flowcell,basecalling_model,notes
1,RUN,run,barcode01,M1,M1,sample,none,,,,Anopheles_gambiae_s.l,Anopheles_coluzzii,Blood_fed,2026-01-01,00:00,SiteA,RegionA,Forest,0,0,Indoor,LTC,Female,E,P,L,SQK,R10,model,
1,RUN,run,barcode02,M2,M2,sample,none,,,,Anopheles_gambiae_s.l,Anopheles_coluzzii,Blood_fed,2026-01-01,00:00,SiteA,RegionA,Forest,0,0,Indoor,LTC,Female,E,P,L,SQK,R10,model,
1,RUN,run,barcode03,M3,M3,sample,none,,,,Anopheles_gambiae_s.l,Anopheles_gambiae_s.s,Blood_fed,2026-01-01,00:00,SiteB,RegionB,Northern_Savannah,0,0,Indoor,LTC,Female,E,P,L,SQK,R10,model,
1,RUN,run,barcode04,M4,M4,sample,none,,,,Anopheles_gambiae_s.l,Anopheles_arabiensis,Blood_fed,2026-01-01,00:00,SiteB,RegionB,Northern_Savannah,0,0,Indoor,LTC,Female,E,P,L,SQK,R10,model,
1,RUN,run,barcode95,NC,NC,negative_control,pcr_negative,,,,,,,,,,,,,,,,,E,P,L,SQK,R10,model,
"""

HOST_CALLS = """sample_uid\trun_id\tsample_id\tbarcode_id\tcontrol_status\tmarker\tbest_cluster_id\thost_assignment\thost_rank\thost_reads\thost_fraction\tn_supporting_features\tbest_feature_id\tbest_confidence\tbest_assignment_status\tmixed_status\ttotal_marker_reads\tmixed_template_warning\tmodel
RUN__barcode01__M1\tRUN\tM1\tbarcode01\tsample\tcyt_b\tcluster001\tHomo sapiens\t1\t30\t1.0\t1\tASV1\thigh\tassigned\tsingle_host\t30\tfalse\tram
RUN__barcode02__M2\tRUN\tM2\tbarcode02\tsample\tcyt_b\tcluster001\tHomo sapiens\t1\t20\t0.8\t1\tASV2\thigh\tassigned\tmixed_host\t25\tfalse\tram
RUN__barcode02__M2\tRUN\tM2\tbarcode02\tsample\tco1_short\tcluster001\tCapra hircus\t2\t5\t0.2\t1\tASV3\thigh\tassigned\tmixed_host\t25\tfalse\tram
RUN__barcode03__M3\tRUN\tM3\tbarcode03\tsample\tcyt_b\tcluster001\tOvis aries\t1\t15\t1.0\t1\tASV4\thigh\tassigned\tsingle_host\t15\tfalse\tram
RUN__barcode04__M4\tRUN\tM4\tbarcode04\tsample\tcyt_b\t\tunassigned\t\t0\t0.0\t0\t\t\tno_model_retained_host_signal\tno_host_signal\t10\tfalse\tram
RUN__barcode95__NC\tRUN\tNC\tbarcode95\tnegative_control\tcyt_b\tcluster001\tHomo sapiens\t1\t99\t1.0\t1\tASV5\thigh\tassigned\tsingle_host\t99\tfalse\tram
"""


def read_tsv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    failures = []
    if hosteco.canonical_species("An_coluzzii") != "Anopheles_coluzzii":
        failures.append("species alias An_coluzzii should canonicalise to samplesheet spelling")

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        samplesheet = tmp / "samplesheet.csv"
        host_calls = tmp / "host_calls.tsv"
        outdir = tmp / "figure_data"
        samplesheet.write_text(SAMPLESHEET)
        host_calls.write_text(HOST_CALLS)

        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--host-calls", str(host_calls), "--samplesheet", str(samplesheet), "--outdir", str(outdir)],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            failures.append(f"compute_host_ecology_indices exited {result.returncode}: {result.stderr}")
        else:
            matrix = read_tsv(outdir / "vector_host_matrix.tsv")
            by_species = {row["species"]: row for row in matrix}
            if by_species["Anopheles_coluzzii"]["Homo sapiens"] != "2":
                failures.append("Anopheles_coluzzii should have two human-positive incidences")
            if by_species["Anopheles_coluzzii"]["Capra hircus"] != "1":
                failures.append("mixed meal should contribute one Capra incidence")

            indices = read_tsv(outdir / "host_ecology_indices.tsv")
            overall_zpi = next((row for row in indices if row["index"] == "zooprophylaxis_index" and row["stratum"] == "Overall"), {})
            if overall_zpi.get("value") != "0.6667":
                failures.append(f"overall zooprophylaxis should be 2/3, got {overall_zpi.get('value')!r}")
            if not any(row["index"] == "niche_breadth_BA" and row["stratum"] == "Anopheles_coluzzii" for row in indices):
                failures.append("niche breadth row for Anopheles_coluzzii is missing")

            overlap = read_tsv(outdir / "niche_overlap_matrix.tsv")
            if overlap[0][overlap[0]["species"]] != "1.0":
                failures.append("Pianka diagonal should be 1.0")
            beta = read_tsv(outdir / "beta_diversity_matrix.tsv")
            if not beta:
                failures.append("beta_diversity_matrix.tsv should contain diagonal/pairwise rows")

    if failures:
        print("host-ecology tests FAILED:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        sys.exit(1)
    print("host-ecology tests PASSED (matrix, zooprophylaxis, niche breadth, overlap, turnover)")


if __name__ == "__main__":
    main()
