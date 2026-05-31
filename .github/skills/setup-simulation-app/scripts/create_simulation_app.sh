#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 <appGroup> <scenario> <programBase>"
  echo "Example: $0 hubbard1d sweepA P_Hubbard1d_TdciPropagate"
  exit 1
fi

app_group="$1"
scenario="$2"
program_base="$3"

if [[ ! "$program_base" =~ ^P_[A-Za-z0-9_]+$ ]]; then
  echo "Error: programBase must start with 'P_' and contain only [A-Za-z0-9_]"
  exit 1
fi

if [[ ! "$app_group" =~ ^[A-Za-z0-9_]+$ ]]; then
  echo "Error: appGroup must contain only [A-Za-z0-9_]"
  exit 1
fi

if [[ ! "$scenario" =~ ^[A-Za-z0-9_]+$ ]]; then
  echo "Error: scenario must contain only [A-Za-z0-9_]"
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../../../.." && pwd)"

target_dir="$repo_root/app/$app_group/$scenario"
mkdir -p "$target_dir"

f90_target="$target_dir/${program_base}.f90"
json_target="$target_dir/${program_base}.json"
json_rel_path="app/$app_group/$scenario/${program_base}.json"

cp "$script_dir/../templates/P_Simulation_Scaffold.f90" "$f90_target"
cp "$script_dir/../templates/P_Simulation_Scaffold.json" "$json_target"

sed -i "s/__PROGRAM_NAME__/${program_base}/g" "$f90_target"
sed -i "s|__JSON_REL_PATH__|${json_rel_path}|g" "$f90_target"

echo "Created scaffold:"
echo "  $f90_target"
echo "  $json_target"
