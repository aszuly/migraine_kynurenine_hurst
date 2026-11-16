#!/bin/bash
# Runs the analysis and writes the results, tables and figures.
#
#   ./run_all.sh                 steps 2-5, from the deposited parcel data
#   ./run_all.sh --from-images   steps 0-5; needs the original workspace file,
#                                the image folders and private/subject_key.csv,
#                                which are not shared
#
# Steps 3-5 also need data/participants.csv (age, aura status and plasma
# metabolite concentrations), which is not part of this repository; see
# README.md. Without it, check_hurst.m recomputes the Hurst exponents from
# the output of step 2, compares them with results/hurst_*.csv, and the
# script stops. Steps 1 and 2 call FSL, found through FSLDIR. The
# interpreters can be set through PYTHON, MATLAB and RSCRIPT, the toolbox
# folder through HURST_TOOLBOXES.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python3}
MATLAB=${MATLAB:-matlab}
RSCRIPT=${RSCRIPT:-Rscript}

if [ "${1:-}" = "--from-images" ]; then
    "$MATLAB" -batch step0_export_inputs
    "$PYTHON" -I step1_extract.py
fi
"$PYTHON" -I step2_denoise.py
if [ ! -f data/participants.csv ]; then
    "$MATLAB" -batch check_hurst
    echo "run_all: step 2 and the Hurst check finished. Steps 3-5 need"
    echo "data/participants.csv (see README.md, section 'Clinical table'); stopping here."
    exit 0
fi
"$MATLAB" -batch step3_main_analysis
"$MATLAB" -batch step4_sensitivity
"$RSCRIPT" step5_figures.R
echo "run_all: finished"
