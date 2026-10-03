#!/bin/bash
# Simulation study 2: RazorHaploNet, median joining and TCS on samples of sequences from a coalescent genealogy.
# Generates the samples, then runs the comparison, one process per sample; a sample that has a report is skipped.
# REPLICATES and PARALLEL can be set in the environment.
#   run-simulation-study-2.sh <output directory>
# JAVA_CP must hold the class path of the app module with its dependencies; iqtree2 must be on the PATH.
set -u
OUT=${1:?output directory}
REPLICATES=${REPLICATES:-10}   # data sets per setting
PARALLEL=${PARALLEL:-1}       # comparison processes at a time, each with four threads
LAST=$(printf rep%03d "$REPLICATES")
run() { java -Xmx6g -Djava.awt.headless=true -cp "$JAVA_CP" "$@"; }
mkdir -p "$OUT"
s=0
for setting in inf finite; do
  s=$((s+1))
  model=""; [ $setting = inf ] && model="--infiniteSites"
  for n in 50 100 200; do
    for theta in 5 10 20 50; do
      d="$OUT/$setting/n${n}t$theta"
      [ -f "$d/$LAST/samples.fasta" ] || run razornet.testsuite.SimulateHaplotypes --output "$d" --samples $n --sites 1000 \
        --theta $theta --replicates $REPLICATES --seed $((s*100000+n*100+theta)) $model
    done
  done
done
one() {
  d=$(dirname "$1")
  [ -s "$d/report.tsv" ] && return
  run razornet.testsuite.CompareHaplotypeMethods --methods Distances,Truth,TCS,RazorHaploNet,RazorHaploNet+Simplify,MedianJoining \
    --truth genealogy.nwk nodes.fasta --restarts 1000 --maxSeconds 20 --threads 4 --timeLimit 120 --silent \
    --networks "$d/nets" --report "$d/report.tmp" "$1" > "$d/run.log" 2>&1 && mv "$d/report.tmp" "$d/report.tsv"
}
export -f one run
export JAVA_CP
(cd "$OUT" && find . -name samples.fasta | sort -t/ -k4,4 -k2,3 | xargs -P "$PARALLEL" -n 1 bash -c 'one "$0"')
