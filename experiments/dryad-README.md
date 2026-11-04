# Data for: RazorNet and RazorHaploNet: haplotype and phylogenetic networks that display distances exactly

Momoko Hayamizu (Waseda University) and Daniel H. Huson (University of Tübingen).

This deposit holds the data sets of the paper as they entered its comparisons, the results of the comparisons
and of the two simulation studies, the workflow files that compute and draw the figures, and the inputs of a
check against PopART. The software, RazorNet, with its source code, the generators of the simulated data and
the scripts that produced everything here, is at <https://github.com/husonlab/razor-net>. The published
data sets that the comparisons start from are available from their original sources, listed below; what is
deposited here is derived from them.

## Directories

| directory | contents |
|---|---|
| `alignments/` | the 26 published alignments as the haplotype network methods saw them: `<name>.haplotypes.fasta` holds the distinct haplotypes on the columns kept, the alignment without every column that contains a gap, a missing or an ambiguous character; `<name>.haplotypes.tsv` lists for each haplotype its number of sequences and their names in the original alignment |
| `distances/` | the seven distance matrices of Table 2 of the paper, in NEXUS format: the numbers of differences between the distinct sequences or profiles (`diamond-princess`, `DImmSV`, `PAegRV`, `agalactiae-virulence-genes`, `huon-languages`), the cytochrome c mutation distances (`cytochrome-c`) and the snow sheep F_ST values (`snow-sheep-fst`) |
| `results/` | one row per data set and method: `published-alignments.tsv` (Table 1 and Supplementary Table S5), `published-distances.tsv` (Table 2), `simulation-study-1.tsv` and its summary (RazorNet against Neighbor-Net on conflicting gene trees), `simulation-study-2.tsv` and its summary (the haplotype methods on samples from a population); `simulation-numbers.tex` holds the numbers quoted in the text |
| `figures/` | the SplitsTree workflow files (`.stree6`) behind Figures 1 to 3, which hold the input, the computed networks and their drawings, and the exported figures (`.pdf`), including those of the simulation studies |
| `popart/` | the seven alignments of Table 1 in the NEXUS dialect that PopART reads, for the check of our TCS construction against PopART |

Column names of the result tables: `file`, `sequences`, `haplotypes`, `sites` (or `ntax`), `method`, `status`
(`ok` or `timeout`), `seconds`, `nodes`, `auxiliary` (nodes added), `edges`, `cycles`, `mutations` (or
`length`), `excess`, `exactPairs`, `maxRatio` (and `minRatio`, `fit`, `splits` for distance networks); the
simulation tables add the measures against the truth that the Methods of the paper define.

## Sources of the published data sets

| source                                                                 | data                                                                                                                   | obtained from                                                                                                            |
|------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| Longdon et al. 2017, *Proc. R. Soc. B* 284:20162381                    | three alignments of insect rhabdoviruses, CCapSV, DImmSV and PAegRV                                                    | the alignments of the paper, Figshare <https://doi.org/10.6084/m9.figshare.3409420>; sequences GenBank KX352837–KX353310 |
| Xu et al. 2026, *Mol. Biol. Evol.* 43:msag052                          | SARS-CoV-2 genomes of the Diamond Princess outbreak and contemporaneous genomes, 482 sequences over 444 variable sites | the supplementary materials of HaploThread, Zenodo <https://doi.org/10.5281/zenodo.18194187>                             |
| Saługa et al. 2022, *J. Biogeogr.* 49:1991–2004                        | twelve alignments of nuclear and plastid markers of the mosses *Syntrichia sarconeurum* and *S. lithophila*            | Dryad <https://doi.org/10.5061/dryad.h18931zp3>                                                                          |
| Estensmo et al. 2021, *Mol. Ecol. Resour.* 21:1141–1152                | ten ITS2 alignments of fungal species                                                                                  | Dryad <https://doi.org/10.5061/dryad.h18931zjq>, `ITSClustering.zip`, `Final_alignments`                                 |
| Morach et al. 2018, *Eur. J. Clin. Microbiol. Infect. Dis.* 37:527–536 | virulence gene profiles of 256 *Streptococcus agalactiae* isolates over 40 genes                                       | as distributed with SplitsTree                                                                                           |
| Greenhill 2015, *PLoS ONE* 10:e0141563                                 | cognate sets of 14 languages of the Huon Peninsula                                                                     | TransNewGuinea.org, as distributed with SplitsTree                                                                       |
| Fitch and Margoliash 1967, *Science* 155:279–284                       | cytochrome c mutation distances between 20 species                                                                     | Table 3 of the paper, as distributed with SplitsTree                                                                     |
| Dotsev et al. 2026, *Ecol. Evol.* 16:e74182                            | pairwise F_ST between 14 populations of the snow sheep                                                                 | Table S1 of the paper's supplementary archive                                                                            |

## How the data were prepared

Every alignment was reduced to the columns without a gap, a missing or an ambiguous character, and identical
sequences were merged into one haplotype; the haplotype network methods received this common input, and the
distance matrices derived from alignments are the numbers of differences between the distinct haplotypes
(`CharactersToDistances` of RazorNet). The matrices of SplitsTree's examples were used as they are. The
protocol and the settings of every method are in the Methods and Supplementary Text S6 of the paper; the
scripts `experiments/run-published-alignments.sh`, `run-published-distances.sh`,
`run-simulation-study-1.sh` and `run-simulation-study-2.sh` of the repository reproduce the results, and
`make-dryad.sh` assembled this deposit. The simulated data sets are not deposited; they are generated again from
the seeds in the scripts.

## Formats and licence

FASTA and NEXUS are read by RazorNet, SplitsTree and PopART; the `.stree6` workflow files are read by RazorNet
and SplitsTree 6; the tables are tab-separated text. The deposited files are released under CC0; the original
data sets remain under the terms of their sources.
