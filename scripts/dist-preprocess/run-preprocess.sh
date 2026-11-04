#!/usr/bin/env bash
# =============================================================================
# run-preprocess.sh — Distance-matrix preprocessing for RazorNet experiments
#
# Edit the "Settings" block below, then run:
#   ./run-preprocess.sh
#
# Processing order (when both MDS and rounding are enabled):
#   ROUND_WHEN=after  -> input -> MDS -> round -> output  (default, v1 behaviour)
#   ROUND_WHEN=before -> input -> round -> MDS -> output
# With only one of MDS or rounding enabled, order is not configurable.
# The original input file is never modified.
#
# Requirements: Java 17+, Maven (mvn), network on first run (downloads ojalgo).
# =============================================================================

set -euo pipefail

# --- Paths (this script lives in scripts/dist-preprocess/) ---------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# =============================================================================
# Quick presets (copy one block into Settings below)
# =============================================================================
#
# MDS only (dimension 2, default):
#   USE_MDS=true   MDS_DIMENSION=2
#   USE_ROUND=false
#   OUTPUT_FILE="/tmp/example.mds2.dist"
#
# Rounding only (2 significant digits):
#   USE_MDS=false
#   USE_ROUND=true   ROUND_DIGITS=2
#   OUTPUT_FILE="/tmp/example.r2.dist"
#
# MDS then rounding (dimension 2, 3 significant digits):
#   USE_MDS=true   MDS_DIMENSION=2
#   USE_ROUND=true   ROUND_DIGITS=3   ROUND_WHEN=after
#   OUTPUT_FILE="/tmp/example.mds2.r3.dist"
#
# Rounding then MDS (dimension 2, 3 significant digits):
#   USE_MDS=true   MDS_DIMENSION=2
#   USE_ROUND=true   ROUND_DIGITS=3   ROUND_WHEN=before
#   OUTPUT_FILE="/tmp/example.r3.mds2.dist"
#
# =============================================================================
# Settings — change these for each experiment
# =============================================================================

# Input distance matrix (.dist). Relative paths are from REPO_ROOT.
INPUT_FILE="tests/example-distances/AmericanIndustry_sqrt1-r.dist"

# Output file. Parent directories are created automatically.
# Naming tips when both MDS and rounding are used:
#   MDS -> round:  foo.mds2.r3.dist
#   round -> MDS:  foo.r3.mds2.dist
OUTPUT_FILE="tests/example-distances/AmericanIndustry_sqrt1-r.mds3.noround.dist"

# --- Classical MDS (optional) ------------------------------------------------
# Set USE_MDS=true to embed distances in d-dimensional Euclidean space and
# recompute pairwise Euclidean (L2) distances.
USE_MDS=true

# Number of MDS dimensions (only used when USE_MDS=true). Default in the Java tool is 2.
MDS_DIMENSION=3

# --- Significant-digit rounding (optional) -----------------------------------
# Set USE_ROUND=true to round every distance to ROUND_DIGITS significant figures.
USE_ROUND=false

# Significant digits (only used when USE_ROUND=true). Example: 3.333333 -> 3.3 at 2 digits.
ROUND_DIGITS=3

# When both USE_MDS and USE_ROUND are true, choose processing order:
#   after  = input -> MDS -> round -> output  (default, same as v1)
#   before = input -> round -> MDS -> output
ROUND_WHEN=after

# =============================================================================
# End of settings — usually no need to edit below this line
# =============================================================================

resolve_path() {
  local p="$1"
  if [[ "$p" = /* ]]; then
    echo "$p"
  else
    echo "$REPO_ROOT/$p"
  fi
}

INPUT="$(resolve_path "$INPUT_FILE")"
OUTPUT="$(resolve_path "$OUTPUT_FILE")"

if [[ ! -f "$INPUT" ]]; then
  echo "error: input file not found: $INPUT" >&2
  echo "  (INPUT_FILE=$INPUT_FILE)" >&2
  exit 1
fi

if [[ "$INPUT" == "$OUTPUT" ]]; then
  echo "error: OUTPUT_FILE must differ from INPUT_FILE (input is never overwritten)" >&2
  exit 1
fi

if [[ "$USE_MDS" != true && "$USE_MDS" != false ]]; then
  echo "error: USE_MDS must be true or false (got: $USE_MDS)" >&2
  exit 1
fi
if [[ "$USE_ROUND" != true && "$USE_ROUND" != false ]]; then
  echo "error: USE_ROUND must be true or false (got: $USE_ROUND)" >&2
  exit 1
fi
if [[ "$USE_MDS" != true && "$USE_ROUND" != true ]]; then
  echo "error: enable at least one of USE_MDS or USE_ROUND (a no-op copy is refused)" >&2
  exit 1
fi
if [[ "$ROUND_WHEN" != after && "$ROUND_WHEN" != before ]]; then
  echo "error: ROUND_WHEN must be 'after' or 'before' (got: $ROUND_WHEN)" >&2
  exit 1
fi

ARGS=(--input "$INPUT" --output "$OUTPUT")

if [[ "$USE_MDS" == true ]]; then
  ARGS+=(--mds --dimension "$MDS_DIMENSION")
fi
if [[ "$USE_ROUND" == true ]]; then
  ARGS+=(--round "$ROUND_DIGITS")
fi
if [[ "$USE_MDS" == true && "$USE_ROUND" == true ]]; then
  ARGS+=(--round-when "$ROUND_WHEN")
fi

echo "=== dist-preprocess run ==="
echo "input:       $INPUT"
echo "output:      $OUTPUT"
echo "USE_MDS:     $USE_MDS"
if [[ "$USE_MDS" == true ]]; then
  echo "MDS_DIMENSION: $MDS_DIMENSION"
fi
echo "USE_ROUND:   $USE_ROUND"
if [[ "$USE_ROUND" == true ]]; then
  echo "ROUND_DIGITS:  $ROUND_DIGITS"
fi
if [[ "$USE_MDS" == true && "$USE_ROUND" == true ]]; then
  echo "ROUND_WHEN:    $ROUND_WHEN"
fi
echo "command:     mvn -q exec:java -Dexec.args=\"${ARGS[*]}\""
echo "---------------------------"

cd "$SCRIPT_DIR"
mvn -q exec:java -Dexec.args="${ARGS[*]}"

echo "Done. Output: $OUTPUT"
