#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <category> <system> <groundstate|dynamics> <method> <programBase>"
  echo "Example: $0 hubbard chain6 dynamics tdci P_TdciPropagate"
  exit 1
fi

category="$1"
system="$2"
stage="$3"
method="$4"
program_base="$5"

if [[ ! "$program_base" =~ ^P_[A-Za-z0-9_]+$ ]]; then
  echo "Error: programBase must start with 'P_' and contain only [A-Za-z0-9_]"
  exit 1
fi

for name in "$category" "$system" "$method"; do
  if [[ ! "$name" =~ ^[A-Za-z0-9_]+$ ]]; then
    echo "Error: category, system and method must contain only [A-Za-z0-9_]"
    exit 1
  fi
done

if [[ "$stage" != "groundstate" && "$stage" != "dynamics" ]]; then
  echo "Error: stage must be 'groundstate' or 'dynamics'"
  exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../../../.." && pwd)"

rel_dir="app/$category/$system/$stage/$method"
target_dir="$repo_root/$rel_dir"
mkdir -p "$target_dir"

f90_target="$target_dir/${program_base}.f90"
json_target="$target_dir/${program_base}.json"
json_rel_path="$rel_dir/${program_base}.json"

cp "$script_dir/../templates/P_Simulation_Scaffold.f90" "$f90_target"
cp "$script_dir/../templates/P_Simulation_Scaffold.json" "$json_target"

sed -i "s/__PROGRAM_NAME__/${program_base}/g" "$f90_target"
sed -i "s|__JSON_REL_PATH__|${json_rel_path}|g" "$f90_target"

echo "Created scaffold:"
echo "  $f90_target"
echo "  $json_target"
