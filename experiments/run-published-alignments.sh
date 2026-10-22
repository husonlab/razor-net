#!/bin/bash
# The haplotype network methods on the published alignments of examples/, with their default settings, one process
# per alignment. Run from the root of the repository:
#   manuscript/experiments/run-published-alignments.sh <output directory>
# JAVA_CP must hold the class path of the app module with its dependencies. Writes results/published-alignments.tsv.
set -u
OUT=${1:?output directory}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT/reports" "$OUT/networks"
S=examples/publications/SalugaEtAl2022/doi_10_5061_dryad_h18931zp3__v20220812
i=0
{
  echo examples/rhabdoviruses/CCapSV.fas
  echo examples/rhabdoviruses/DImmSV.fas
  echo examples/rhabdoviruses/PAegRV.fas
  echo $S/ALIGNMENT_ITS_Haplotype_Network_116_SEQ.fasta
  echo $S/ALIGNMENT_cpDNA_Haplotype_Network_111_SEQ.fasta
  echo examples/publications/XuEtAl2026/Diamond_Princess.phy
  echo $S/ALIGNMENT_ITS_Population_Genetics_120_SEQ.fasta
  echo $S/ALIGNMENT_rps4_Population_Genetics_123_SEQ.fasta
  echo $S/ALIGNMENT_tRNAGly_Population_Genetics_111_SEQ.fas
  echo $S/ALIGNMENT_trnLF_Population_Genetics_117_SEQ.fasta
  echo $S/ALIGNMENT_Concatenated_Population_Genetics_120_SEQ.fasta
  echo $S/ALIGNMENT_Continental_Antarctic_Population_Genetics_58_SEQ.fas
  echo $S/ALIGNMENT_Maritime_Antarctic_Population_Genetics_49_SEQ.fas
  echo $S/ALIGNEMNT_Patagonia_Population_Genetics_10_SEQ.fas
  echo $S/ALIGNMENT_ITS_Phylogeography_127_SEQ.fasta
  echo $S/ALIGNMENT_cpDNA_Phylogeography_115_SEQ.fasta
  ls examples/Issakka_ITSClustering/Final_alignments/*.fasta
} | while IFS= read -r alignment; do
  i=$((i+1))
  report=$OUT/reports/$(printf %02d $i).tsv
  [ -s "$report" ] && continue
  java -Xmx6g -Djava.awt.headless=true -cp "$JAVA_CP" razornet.testsuite.CompareHaplotypeMethods \
    --methods TCS,RazorHaploNet,RazorHaploNet+Shave,MedianJoining --restarts 1000 --maxSeconds 60 --threads 0 \
    --timeLimit 300 --silent --networks "$OUT/networks" --report "$report.tmp" "$alignment" > "${report%.tsv}.log" 2>&1 \
    && mv "$report.tmp" "$report"
done
# one table of all reports; an alignment with a single haplotype has no report
{ head -1 "$OUT/reports/01.tsv"; for f in "$OUT"/reports/[0-9][0-9].tsv; do tail -n +2 "$f"; done; } > "$HERE/results/published-alignments.tsv"
