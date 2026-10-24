#!/bin/bash
# Assembles the Dryad deposit of the paper: the data sets as they entered the comparisons, the results, the
# workflow files and exports of the figures, the PopART inputs, and the README. Run from the root of the
# repository:
#   manuscript/experiments/make-dryad.sh <output directory> [<splitstree6 examples directory>]
# JAVA_CP must hold the class path of the app module with its dependencies.
set -eu
OUT=${1:?output directory}
ST=${2:-../splitstree6/examples}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT/alignments" "$OUT/distances" "$OUT/results" "$OUT/figures" "$OUT/popart"
S=examples/publications/SalugaEtAl2022/doi_10_5061_dryad_h18931zp3__v20220812

# the 26 alignments as the haplotype methods saw them: the distinct haplotypes on the kept columns
# (<name>.haplotypes.fasta) and the sequences of each haplotype (<name>.haplotypes.tsv)
alignments=(
  examples/rhabdoviruses/CCapSV.fas
  examples/rhabdoviruses/DImmSV.fas
  examples/rhabdoviruses/PAegRV.fas
  "$S/ALIGNMENT_ITS_Haplotype_Network_116_SEQ.fasta"
  "$S/ALIGNMENT_cpDNA_Haplotype_Network_111_SEQ.fasta"
  examples/publications/XuEtAl2026/Diamond_Princess.phy
  "$S/ALIGNMENT_ITS_Population_Genetics_120_SEQ.fasta"
  "$S/ALIGNMENT_rps4_Population_Genetics_123_SEQ.fasta"
  "$S/ALIGNMENT_tRNAGly_Population_Genetics_111_SEQ.fas"
  "$S/ALIGNMENT_trnLF_Population_Genetics_117_SEQ.fasta"
  "$S/ALIGNMENT_Concatenated_Population_Genetics_120_SEQ.fasta"
  "$S/ALIGNMENT_Continental_Antarctic_Population_Genetics_58_SEQ.fas"
  "$S/ALIGNMENT_Maritime_Antarctic_Population_Genetics_49_SEQ.fas"
  "$S/ALIGNEMNT_Patagonia_Population_Genetics_10_SEQ.fas"
  "$S/ALIGNMENT_ITS_Phylogeography_127_SEQ.fasta"
  "$S/ALIGNMENT_cpDNA_Phylogeography_115_SEQ.fasta"
  # the ten fungal alignments of the paper; the eleventh of the deposit, Amylap, has a single haplotype
  "examples/Issakka_ITSClustering/Final_alignments/Antser_h_seq_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Fompin_h_seq_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Fomros_h_asv_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Glosep_h_asv_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Phecen_h_asv_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Phefer_h_seq_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Phenig_h_seq_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Phevit_h_seq_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Poscae_h_seq_Nucleotide alignment.fasta"
  "examples/Issakka_ITSClustering/Final_alignments/Triabi_h_seq_Nucleotide alignment.fasta"
)
for alignment in "${alignments[@]}"; do   # one call per alignment, as run-published-alignments.sh does
  java -Djava.awt.headless=true -cp "$JAVA_CP" razornet.testsuite.CompareHaplotypeMethods \
    --methods Distances --networks "$OUT/alignments" --silent "$alignment" > /dev/null
done
for f in "$OUT"/alignments/*\ *; do [ -e "$f" ] && mv "$f" "${f// /_}"; done

# the seven distance matrices of Table 2 as RazorNet and Neighbor-Net saw them
convert() { # name, alignment
  java -Djava.awt.headless=true -cp "$JAVA_CP" razornet.testsuite.CharactersToDistances "$2" "$OUT/distances/$1.nex"
}
convert diamond-princess examples/publications/XuEtAl2026/Diamond_Princess.phy
convert DImmSV examples/rhabdoviruses/DImmSV.fas
convert PAegRV examples/rhabdoviruses/PAegRV.fas
convert agalactiae-virulence-genes "$ST/publications/MorachEtAl2018/streptococcus-agalactiae.nex"
convert huon-languages "$ST/publications/Greenhill2015/Huon-languages.nex"
cp "$ST/publications/FitchMargoliash1967/table3.nex" "$OUT/distances/cytochrome-c.nex"
cp examples/publications/DotsevEtAl2026/snow-sheep-fst.nex "$OUT/distances/snow-sheep-fst.nex"

# the results of the comparisons and the simulation studies
cp manuscript/experiments/results/*.tsv manuscript/experiments/results/simulation-numbers.tex "$OUT/results/"

# the figures: the workflow files that compute and draw them, and the exports
cp examples/publications/SalugaEtAl2022/ALIGNMENT_ITS_Haplotype_Network_116_SEQ-Figure1a-PopART-TCS.stree6 \
   examples/publications/SalugaEtAl2022/ALIGNMENT_ITS_Haplotype_Network_116_SEQ-haplotypes-fig1b.stree6 \
   examples/publications/SalugaEtAl2022/ALIGNMENT_ITS_Haplotype_Network_116_SEQ-haplotypes-fig1c.stree6 \
   examples/publications/XuEtAl2026/Diamond_Princess-four-razornet.stree6 \
   examples/publications/DotsevEtAl2026/snow-sheep-fst.stree6 \
   manuscript/figures/*.pdf "$OUT/figures/"

# the PopART inputs of the check against PopART
cp examples/popart-check/*.nex "$OUT/popart/"

cp "$HERE/dryad-README.md" "$OUT/README.md"
echo "written to $OUT:"; du -sh "$OUT"/*
