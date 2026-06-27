# HÆMA pipeline — common developer / release tasks.
# Usage: `make help`. Requires Nextflow, Docker, Java 17-24, python3.

NXF        ?= nextflow
PROFILE    ?= test,docker
OUTDIR     ?= results/make_test
PY_IMAGE   ?= haema-python:0.3.0
R_IMAGE    ?= haema-r:0.3.0
FIG_IMAGE  ?= haema-figures:0.4.0
MEDAKA_IMG ?= haema-medaka:0.3.0
# Set REGISTRY to your registry namespace to push, e.g. REGISTRY=ghcr.io/USER
REGISTRY   ?=

.PHONY: help
help:  ## Show this help
	@grep -E '^[a-zA-Z0-9_-]+:.*?## ' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-22s\033[0m %s\n", $$1, $$2}'

## ---- validation / tests (no containers needed except `test`) ----
.PHONY: lint
lint: pycompile schema config verify-refs unit  ## Run all static checks (python, schema, config, refs, unit tests)

.PHONY: pycompile
pycompile:  ## Byte-compile all bin/*.py and tests/*.py
	python3 -m py_compile bin/*.py tests/*.py && echo "python: OK"

.PHONY: unit
unit:  ## Run python unit tests (validation, positive-control, LCA, taxid, ecological indices)
	python3 tests/test_validation.py && python3 tests/test_positive_controls.py \
	  && python3 tests/test_lca.py && python3 tests/test_taxid_assignment.py \
	  && python3 tests/test_ecological_indices.py && python3 tests/test_host_ecology_indices.py

.PHONY: schema
schema:  ## Validate nextflow_schema.json is valid JSON
	python3 -c "import json; json.load(open('nextflow_schema.json')); print('schema: OK')"

.PHONY: config
config:  ## Parse nextflow.config for all profiles
	@for p in test local docker singularity apptainer slurm gpu production; do \
	  $(NXF) config . -profile $$p >/dev/null 2>&1 && echo "config $$p: OK" || { echo "config $$p: FAIL"; exit 1; }; \
	done

.PHONY: verify-refs
verify-refs:  ## Verify bundled reference checksums + sidecar fields
	python3 bin/verify_reference_assets.py --assets-dir assets/references

.PHONY: test
test:  ## Run the bundled test profile (curated BLAST taxonomy)
	$(NXF) run . -profile $(PROFILE) --skip_taxonomy false --outdir $(OUTDIR) --log_dir logs/make_test

.PHONY: stub
stub:  ## Whole-DAG stub run (all optional features on)
	$(NXF) run . -profile test,docker -stub-run \
	  --skip_taxonomy false --enable_medaka true --medaka_container python:3.11 \
	  --enable_r_outputs true --enable_multiqc true --multiqc_container python:3.11 \
	  --enable_figures true --enable_publication_figures true --figures_container python:3.11 \
	  --blast_container python:3.11 --outdir results/make_stub --log_dir logs/make_stub

.PHONY: prerelease
prerelease:  ## Run the full pre-release validation script
	bash tests/validate_release.sh

## ---- custom container images (production) ----
.PHONY: images
images: image-python image-r image-figures  ## Build the custom images (denoising, R, figures)

image-python:  ## Build haema-python
	docker build -t $(PY_IMAGE) -f containers/haema-python/Dockerfile .
image-r:  ## Build haema-r
	docker build -t $(R_IMAGE) -f containers/haema-r/Dockerfile .
image-figures:  ## Build haema-figures (publication figure stack)
	docker build -t $(FIG_IMAGE) -f containers/haema-figures/Dockerfile .
image-medaka:  ## Build the OPTIONAL model-asserted Medaka image
	docker build -t $(MEDAKA_IMG) -f containers/haema-medaka/Dockerfile .

.PHONY: push
push:  ## Tag+push custom images to $REGISTRY (set REGISTRY=...). Prints digests to pin.
	@test -n "$(REGISTRY)" || { echo "Set REGISTRY=<namespace>, e.g. REGISTRY=ghcr.io/USER"; exit 1; }
	docker tag $(PY_IMAGE)  $(REGISTRY)/$(PY_IMAGE)  && docker push $(REGISTRY)/$(PY_IMAGE)
	docker tag $(R_IMAGE)   $(REGISTRY)/$(R_IMAGE)   && docker push $(REGISTRY)/$(R_IMAGE)
	docker tag $(FIG_IMAGE) $(REGISTRY)/$(FIG_IMAGE) && docker push $(REGISTRY)/$(FIG_IMAGE)
	@echo "Now pin the printed @sha256: digests in nextflow.config (params.python_container / r_container / figures_container)."
