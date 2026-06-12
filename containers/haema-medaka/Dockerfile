# OPTIONAL hardened Medaka image for HÆMA consensus polishing.
#
# The production profile uses the pinned upstream ONT image directly
# (ontresearch/medaka:sha<...>), and the model is verified at RUNTIME by
# MEDAKA_MODEL_PREFLIGHT, so this image is not required. It exists only for users
# who additionally want the required model asserted at BUILD time, or who want to
# bake the model into a frozen image for an air-gapped HPC.
#
# The base defaults to a pinned immutable upstream sha tag (NOT :latest) so the
# build is reproducible. Override MEDAKA_BASE_IMAGE to pin a different tag/digest.
ARG MEDAKA_BASE_IMAGE=ontresearch/medaka:shaf39439188c323053e66cf79618fae1ab33d1b38a
FROM ${MEDAKA_BASE_IMAGE}
ARG MEDAKA_BASE_IMAGE

LABEL org.opencontainers.image.title="haema-medaka"
LABEL org.opencontainers.image.description="Medaka runtime for HÆMA ONT consensus polishing"
LABEL org.opencontainers.image.version="0.3.0"
LABEL org.opencontainers.image.base.name="${MEDAKA_BASE_IMAGE}"
LABEL haema.medaka.required_model="r1041_e82_400bps_sup_v4.3.0"

# Assert the required SUP/R10 model is present at build time.
RUN medaka tools list_models | grep -F "r1041_e82_400bps_sup_v4.3.0"
