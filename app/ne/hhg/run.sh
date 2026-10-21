#!/usr/bin/env bash
# Run the neon strong-field MCTDHF pipeline with the beryllium programs.
#
# usage: run.sh <inputDir> <runDir> [step ...]
#   steps: relaxA relaxB mctdhf_<pulse>          (default: all)
#
#   relaxA                    imaginary-time MCTDHF at lmax = 1
#   relaxB    <- relaxA       embedding into the propagation grid (one check)
#   mctdhf_*  <- relaxB       MCTDHF in the laser pulse
#
# Requires the executables (cody build app) and `source setBuildVars.sh`.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
programs="$here/../../be"
inputs="$(realpath "$1")"
runs="$(realpath -m "$2")"
shift 2
steps=("$@")
if [ ${#steps[@]} -eq 0 ]; then
  steps=(relaxA relaxB mctdhf_strong2c)
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
  if ! (cd "$dir" && { time "$programs/$exe" "$inputs/$json"; } > log.txt 2>&1); then
    echo "[$(date '+%F %T')] $step FAILED (see $dir/log.txt)"
    exit 1
  fi
  echo "[$(date '+%F %T')] $step done"
}

for step in "${steps[@]}"; do
  case $step in
    relaxA)   run_step relaxA relax/P_ImagTimeMctdhf.exe stageA.json "" ;;
    relaxB)   run_step relaxB relax/P_ImagTimeMctdhf.exe stageB.json relaxA "orb*.in" "coeffs.in" ;;
    mctdhf_*) run_step "$step" laser/P_LaserMctdhf.exe "$step.json" relaxB "orb*.in" "coeffs.in" ;;
    *) echo "unknown step: $step" >&2; exit 1 ;;
  esac
done
