#!/bin/bash
# Collects the reports of the two simulation studies and writes the tables and figures of the manuscript.
#   make-results.sh <directory of study 1> <directory of study 2>
# Needs python3 with matplotlib (set PYTHON to choose the interpreter).
set -e
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python3}
$PYTHON collect.py "$1" results/simulation-study-1.tsv
$PYTHON collect.py "$2" results/simulation-study-2.tsv
$PYTHON figures.py results/simulation-study-1.tsv results/simulation-study-2.tsv figures results
$PYTHON tables.py results tables
