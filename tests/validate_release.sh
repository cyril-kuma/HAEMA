#!/usr/bin/env bash
# Pre-release validation for the HÆMA pipeline. Runs every check that does not need external
# data, registry credentials, or a GPU. Safe to run locally and in CI.
#
#   bash tests/validate_release.sh            # static checks + fail-fast tests (fast, no Docker)
#   bash tests/validate_release.sh --run      # also runs the test profile (needs Docker)
#
# Requires: nextflow (Java 17-24 on PATH or JAVA_HOME), python3. Run from the repo root.
set -uo pipefail
cd "$(dirname "$0")/.."

NXF="${NXF:-nextflow}"
RUN_PIPELINE=0
[[ "${1:-}" == "--run" ]] && RUN_PIPELINE=1
fail=0
section() { printf '\n=== %s ===\n' "$1"; }
check()   { if "$@"; then echo "  PASS"; else echo "  FAIL: $*"; fail=1; fi; }

section "1. Python byte-compile + unit tests"
check python3 -m py_compile bin/*.py tests/*.py
check python3 tests/test_validation.py
check python3 tests/test_positive_controls.py
check python3 tests/test_lca.py
check python3 tests/test_taxid_assignment.py
check python3 tests/test_ecological_indices.py
check python3 tests/test_host_ecology_indices.py

section "2. Schema is valid JSON and covers required inputs"
check python3 -c "import json; d=json.load(open('nextflow_schema.json')); io=next(x['\$ref'] for x in d['allOf'] if x['\$ref'].endswith('/input_output')); block=d['\$defs'][io.rsplit('/',1)[-1]]; assert {'input','raw_data_dir'} <= set(block['required']); print('  schema OK, input/output required fields:', ', '.join(block['required']))"

section "3. nextflow.config parses for all profiles"
for p in test local docker singularity apptainer slurm gpu production; do
  if $NXF config . -profile "$p" >/dev/null 2>&1; then echo "  $p: OK"; else echo "  $p: FAIL"; fail=1; fi
done

section "4. Bundled reference assets (checksums + sidecar fields)"
check python3 bin/verify_reference_assets.py --assets-dir assets/references

section "5. Fail-fast on missing required inputs"
# Capture output first: Nextflow exits non-zero on the error, which would trip pipefail in a pipe.
ff_out="$($NXF run . --outdir /tmp/haema_failfast 2>&1 || true)"
if grep -Eq "Missing (--input|required parameter\\(s\\): input)" <<<"$ff_out"; then
  echo "  PASS (clear missing-input error)"
else
  echo "  FAIL (no clear missing-input error)"; fail=1
fi

section "6. Release artifacts present"
for f in LICENSE CITATION.cff CHANGELOG.md CONTRIBUTING.md CODE_OF_CONDUCT.md \
         .github/workflows/ci.yml docs/usage.md docs/output.md docs/limitations.md; do
  [[ -e "$f" ]] && echo "  $f: OK" || { echo "  $f: MISSING"; fail=1; }
done

if [[ "$RUN_PIPELINE" == "1" ]]; then
  section "7. Test-profile run (curated BLAST taxonomy)"
  if $NXF run . -profile test,docker --skip_taxonomy false \
        --outdir results/validate_release --log_dir logs/validate_release >/dev/null 2>&1; then
    echo "  run: OK"
    python3 tests/test_endpoint_columns.py results/validate_release || fail=1
  else
    echo "  run: FAIL"; fail=1
  fi
fi

section "RESULT"
if [[ "$fail" == "0" ]]; then echo "  ALL CHECKS PASSED"; else echo "  SOME CHECKS FAILED"; fi
exit "$fail"
