#!/bin/bash
# Simulation study 1: RazorNet and Neighbor-Net on sequences that evolved down gene trees which each differ from a
# species tree by one subtree-prune-and-regraft move. Generates the data sets, then runs the comparison, one
# process per data set, three at a time; a data set that has a report is skipped, so the script can be restarted.
# REPLICATES and PARALLEL can be set in the environment.
#   run-simulation-study-1.sh <output directory>
# JAVA_CP must hold the class path of the app module with its dependencies; iqtree2 must be on the PATH.
set -u
OUT=${1:?output directory}
REPLICATES=${REPLICATES:-10}   # data sets per setting
PARALLEL=${PARALLEL:-3}       # comparison processes at a time, each with four threads
LAST=$(printf rep%03d "$REPLICATES")
run() { java -Xmx6g -Djava.awt.headless=true -cp "$JAVA_CP" "$@"; }
mkdir -p "$OUT"
s=0
for setting in inf low mid high; do
  s=$((s+1))
  case $setting in
    inf)  model="--height 0.01 --length 6000 --infiniteSites";;
    low)  model="--height 0.002 --length 30000";;
    mid)  model="--height 0.01 --length 6000";;
    high) model="--height 0.05 --length 1200";;
  esac
  for n in 10 20 50 100; do
    for k in 1 2 3 4 5 6; do
      d="$OUT/$setting/n${n}k$k"
      [ -f "$d/$LAST/hamming.dist" ] || run razornet.testsuite.SimulateTreeMixtures --output "$d" --taxa $n --trees $k \
        --replicates $REPLICATES --seed $((s*100000+n*100+k)) $model
    done
  done
done
one() {
  d=$(dirname "$1")
  [ -s "$d/report.tsv" ] && return
  run razornet.testsuite.CompareDistanceMethods --truth genes.mutations.nwk --restarts 1000 --maxSeconds 20 --threads 4 \
    --timeLimit 120 --silent --networks "$d/nets" --report "$d/report.tmp" "$1" > "$d/run.log" 2>&1 && mv "$d/report.tmp" "$d/report.tsv"
}
export -f one run
export JAVA_CP
(cd "$OUT" && find . -name hamming.dist | sort -t/ -k4,4 -k2,3 | xargs -P "$PARALLEL" -n 1 bash -c 'one "$0"')
