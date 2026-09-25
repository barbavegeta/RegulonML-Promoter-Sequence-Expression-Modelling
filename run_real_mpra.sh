#!/usr/bin/env bash
# Full RegulonML run: download the Kircher et al. 2019 MPRA data, build sequence and
# JASPAR motif features, evaluate with leave-one-locus-out CV, and write figures.
set -euo pipefail
export PYTHONPATH="$(cd "$(dirname "$0")" && pwd)/src${PYTHONPATH:+:$PYTHONPATH}"
python3 -m regulonml.pipeline --outdir results "$@"
