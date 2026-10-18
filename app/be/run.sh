#!/usr/bin/env bash
# Run the beryllium HHG pipeline (MCTDHF vs TD-2RDM, Lackner et al. PRA 95, 033414).
#
# usage: run.sh <inputDir> <runDir> [step ...]
#   steps: relaxA relaxB average mctdhf_<pulse> td2rdm_<pulse>
#          relaxA_hf relaxB_hf tdhf_<pulse>                     (default: all)
#   pulses: weak4c strong2c
#
# Each step runs in <runDir>/<step> and copies its input files from the step it
# depends on:
#   relaxA                    imaginary-time MCTDHF at small lmax
#   relaxB    <- relaxA       imaginary-time MCTDHF on the propagation grid
#   average   <- relaxB       TD-2RDM initial state (field-free time average)
#   mctdhf_*  <- relaxB       MCTDHF in the laser pulse
#   td2rdm_*  <- average      TD-2RDM in the laser pulse
#   relaxA_hf, relaxB_hf, tdhf_*   the same for TDHF (MCTDHF with 1s, 2s only)
#
# Requires the executables (cody build app) and `source setBuildVars.sh`.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
inputs="$(realpath "$1")"
runs="$(realpath -m "$2")"
shift 2
steps=("$@")
if [ ${#steps[@]} -eq 0 ]; then
  steps=(relaxA relaxB average mctdhf_weak4c td2rdm_weak4c mctdhf_strong2c td2rdm_strong2c
         relaxA_hf relaxB_hf tdhf_weak4c tdhf_strong2c)
fi

export GFORTRAN_UNBUFFERED_PRECONNECTED=y
export CODY_PRETTY_PRINT=OFF

run_step() {
  local step=$1 exe=$2 json=$3 from=$4
  shift 4
  local dir="$runs/$step"
  mkdir -p "$dir"
  if [ -n "$from" ]; then
    for pattern in "$@"; do
      cp "$runs/$from"/$pattern "$dir"/
    done
  fi
  echo "[$(date '+%F %T')] $step: $exe $json"
  if ! (cd "$dir" && { time "$here/$exe" "$inputs/$json"; } > log.txt 2>&1); then
    echo "[$(date '+%F %T')] $step FAILED (see $dir/log.txt)"
    exit 1
  fi
  echo "[$(date '+%F %T')] $step done"
}

for step in "${steps[@]}"; do
  case $step in
    relaxA)   run_step relaxA relax/P_ImagTimeMctdhf.exe stageA.json "" ;;
    relaxB)   run_step relaxB relax/P_ImagTimeMctdhf.exe stageB.json relaxA "orb*.in" "coeffs.in" ;;
    average)  run_step average td2rdmInit/P_Td2rdmTimeAverage.exe timeAverage.json relaxB "orb*.in" "rdm2Bt*.in" ;;
    mctdhf_*) run_step "$step" laser/P_LaserMctdhf.exe "$step.json" relaxB "orb*.in" "coeffs.in" ;;
    td2rdm_*) run_step "$step" laser/P_LaserTd2rdm.exe "$step.json" average "orb*.in" "rdm2Bt*.in" ;;
    relaxA_hf) run_step relaxA_hf relax/P_ImagTimeMctdhf.exe stageA_hf.json "" ;;
    relaxB_hf) run_step relaxB_hf relax/P_ImagTimeMctdhf.exe stageB_hf.json relaxA_hf "orb*.in" "coeffs.in" ;;
    tdhf_*)   run_step "$step" laser/P_LaserMctdhf.exe "$step.json" relaxB_hf "orb*.in" "coeffs.in" ;;
    *) echo "unknown step: $step" >&2; exit 1 ;;
  esac
done
