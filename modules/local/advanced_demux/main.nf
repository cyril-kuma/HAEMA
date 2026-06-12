process ADVANCED_DEMUX {
    label 'process_high'
    tag "${params.advanced_demux_tool}:${params.demux_run_id}"

    publishDir "${params.outdir}/00_demultiplexing", mode: 'copy'

    input:
    path pooled_fastqs

    output:
    path 'demuxed_raw', emit: demux_root
    path 'advanced_demux_summary.tsv', emit: summary
    path 'versions.yml', emit: versions

    when:
    params.enable_advanced_demux || params.input_type == 'pooled_fastq'

    script:
    def tool = params.advanced_demux_tool as String
    def command_template = params.advanced_demux_command_template as String
    """
    mkdir -p demuxed_raw
    printf '%s\\n' ${pooled_fastqs} > pooled_fastqs.list

    if [[ '${tool}' == 'header_tag' ]]; then
        python3 ${projectDir}/bin/demux_fastq_by_header.py \\
            --fastq-list pooled_fastqs.list \\
            --run-id '${params.demux_run_id}' \\
            --output-root demuxed_raw \\
            --summary advanced_demux_summary.tsv \\
            --min-reads-per-barcode '${params.advanced_demux_min_reads_per_barcode}'
    elif [[ '${tool}' == 'command_template' || '${tool}' == 'barbell' || '${tool}' == 'deepbinner' ]]; then
        if [[ -z '${command_template}' ]]; then
            cat >&2 <<'END_ERROR'
Advanced demultiplexing was enabled with a command-wrapper tool, but
--advanced_demux_command_template is empty.

Provide a command that writes MinKNOW-like FASTQ output under {output}.
Available placeholders: {input}, {input_list}, {output}, {barcode_kit}, {threads}.
END_ERROR
            exit 2
        fi
        cat > command_template.txt <<'END_TEMPLATE'
${command_template}
END_TEMPLATE
        python3 - <<'PY'
from pathlib import Path
template = Path("command_template.txt").read_text().strip()
cmd = template.format(
    input=str(Path("pooled_fastqs.list").read_text().splitlines()[0]),
    input_list="pooled_fastqs.list",
    output="demuxed_raw",
    barcode_kit="${params.barcode_kit}",
    threads="${task.cpus}",
)
Path("run_advanced_demux.sh").write_text("#!/usr/bin/env bash\\nset -euo pipefail\\n" + cmd + "\\n")
PY
        bash run_advanced_demux.sh
        if [[ ! -d demuxed_raw ]]; then
            echo "Advanced demux command completed but did not create demuxed_raw." >&2
            exit 2
        fi
        python3 - <<'PY'
from pathlib import Path
import csv
rows = []
for fastq in sorted(Path("demuxed_raw").rglob("*.fastq*")):
    barcode = next((p.name for p in fastq.parents if p.name.startswith("barcode")), "unknown")
    rows.append({
        "run_id": "${params.demux_run_id}",
        "source_file": fastq.name,
        "barcode_id": barcode,
        "reads": "",
        "retained": "true",
        "method": "${tool}",
    })
with Path("advanced_demux_summary.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["run_id", "source_file", "barcode_id", "reads", "retained", "method"], delimiter="\\t")
    writer.writeheader()
    writer.writerows(rows)
PY
    else
        echo "Unsupported --advanced_demux_tool '${tool}'. Use header_tag, command_template, barbell, or deepbinner." >&2
        exit 2
    fi

    cat > versions.yml <<-END_VERSIONS
    "${task.process}":
      wrapper: "demux_fastq_by_header.py"
      advanced_demux_tool: "${tool}"
    END_VERSIONS
    """

    stub:
    """
    mkdir -p demuxed_raw/${params.demux_run_id}/advanced_demux/fastq_pass/barcode01
    cat > demuxed_raw/${params.demux_run_id}/advanced_demux/fastq_pass/barcode01/${params.demux_run_id}_barcode01.advanced_demux.fastq <<-END_FASTQ
    @stub barcode=barcode01
    ACGTACGTACGT
    +
    IIIIIIIIIIII
    END_FASTQ
    cat > advanced_demux_summary.tsv <<-END
    run_id	source_file	barcode_id	reads	retained	method
    ${params.demux_run_id}	stub.fastq	barcode01	1	true	stub
    END
    cat > versions.yml <<-END
    "${task.process}":
      stub: true
    END
    """
}
