#!/bin/bash
# Neighbor-Net, RazorNet and RazorNet with the Simplify filter on published distance data (Table 2 of the paper and
# the candidates listed in its TODO). Alignments and character matrices enter as the Hamming distances of their
# distinct sequences or profiles, prepared as for the haplotype comparison (CharactersToDistances); the distance
# matrices of SplitsTree's examples enter as they are. Run from the root of the repository:
#   manuscript/experiments/run-published-distances.sh <output directory> [<splitstree6 examples directory>]
# JAVA_CP must hold the class path of the app module with its dependencies. Writes results/published-distances.tsv.
set -u
OUT=${1:?output directory}
ST=${2:-../splitstree6/examples}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT/distances" "$OUT/networks"
S=examples/publications/SalugaEtAl2022/doi_10_5061_dryad_h18931zp3__v20220812
convert() { # name, alignment
  [ -s "$OUT/distances/$1.nex" ] || java -Djava.awt.headless=true -cp "$JAVA_CP" razornet.testsuite.CharactersToDistances "$2" "$OUT/distances/$1.nex"
}
convert diamond-princess examples/publications/XuEtAl2026/Diamond_Princess.phy
convert DImmSV examples/rhabdoviruses/DImmSV.fas
convert PAegRV examples/rhabdoviruses/PAegRV.fas
convert CCapSV examples/rhabdoviruses/CCapSV.fas
convert ITS-haplotype-network "$S/ALIGNMENT_ITS_Haplotype_Network_116_SEQ.fasta"
convert ITS-phylogeography "$S/ALIGNMENT_ITS_Phylogeography_127_SEQ.fasta"
convert cpDNA-phylogeography "$S/ALIGNMENT_cpDNA_Phylogeography_115_SEQ.fasta"
convert agalactiae-virulence-genes "$ST/publications/MorachEtAl2018/streptococcus-agalactiae.nex"
convert huon-languages "$ST/publications/Greenhill2015/Huon-languages.nex"
cp -n "$ST/publications/FitchMargoliash1967/table3.nex" "$OUT/distances/cytochrome-c.nex" 2>/dev/null
cp -n "$ST/publications/Additional/SaitouNei1987-Table1.nex" "$OUT/distances/saitou-nei-1987.nex" 2>/dev/null
cp -n "$ST/publications/WolfEtAl2021/Cochlearia-distances.nex" "$OUT/distances/cochlearia.nex" 2>/dev/null
cp -n "$ST/publications/McVayEtAl2017/gen-2016-0206supplb.txt" "$OUT/distances/oaks.nex" 2>/dev/null
cp -n examples/publications/DotsevEtAl2026/snow-sheep-fst.nex "$OUT/distances/snow-sheep.nex" 2>/dev/null
java -Xmx8g -Djava.awt.headless=true -cp "$JAVA_CP" razornet.testsuite.CompareDistanceMethods \
  --methods NeighborNet,RazorNet,RazorNet+Simplify --restarts 100 --maxSeconds 60 --threads 4 --seed 0 --timeLimit 120 \
  --networks "$OUT/networks" --report "$OUT/report.tsv" "$OUT"/distances/*.nex > "$OUT/log.txt" 2>&1
sed "s#$OUT/distances/##" "$OUT/report.tsv" > "$HERE/results/published-distances.tsv"
