# The example data

Every file under `examples/`, what it is and where it comes from. "Ours" means written by the authors or by the
program for this repository. Files that the program writes from an input (`.hamming.dist`, `.gml`, `.stree6`
workflows) are derived from the input next to them and carry its licence.

The published data sets are under `publications/`, `rhabdoviruses/` and `Issakka_ITSClustering/`; each
directory has a `README` or a `Citation.txt` naming the paper and the deposit. Please cite the source paper when
you use its data.

## Top level

| file | what it is | source |
|---|---|---|
| `ImrichEtAl1984.dist` | a metric on 24 points | the example of Imrich, Simões-Pereira and Zamfirescu 1984, *On optimal embeddings of metrics in graphs*, J. Combin. Theory B 36:1–15, transcribed from the paper |
| `ImrichEtAl1984.stree6` | the matrix above with RazorNet and the Shave filter | ours, derived |
| `cube.dist` | the 8 vertices of the unit 3-cube under the Hamming metric; the cube itself is the optimal realization (length 12), as the header says | ours, generated |
| `five.dist`, `five.gml`, `five.stree6` | the five-point metric of the paper's concept figure (Supplementary Figure S7, the idea behind RazorNet), RazorNet's network as GML, and the workflow | ours, made up |
| `handwritten11x11.dist`, `handwritten11x11.gml` | a metric on 11 points a–k and RazorNet's network | MH's hand-worked example of the Strict strategy |
| `manhattan_distances_20x20-points.tsv`, `manhattan_distances_20x20.dist`, `manhattan_distances_20x20.stree6` | the grid example of the paper: 20 random points on an integer grid, their Manhattan distances, and the workflow with Neighbor-Net, RazorNet and the Shave filter (Supplementary Figure S6) | ours, generated |
| `popart-example.stree6` | PopART's example NEXUS file (7 sequences × 56 sites with traits, modified from a DnaSP sample; the source URL is in the file), with Median-Joining and RazorHaploNet | PopART's example data, <http://popart.otago.ac.nz/examplenex.shtml> |

## `bundled/`

The six files of the program's File → Examples menu, small and quick to compute: the three distance matrices
open with RazorNet, the three alignments with RazorHaploNet. `scripts/make-examples.sh examples/bundled
app/src/main/resources/razornet/examples` copies them into the program's resources, as `.dat` files with an
`index.dat` that keeps the menu order and the display names.

| file | what it is | source |
|---|---|---|
| `snow-sheep-fst.nex` | `publications/DotsevEtAl2026/snow-sheep-fst.nex` with the display labels of `snow-sheep-fst-labels.txt` in its TAXA block | Dotsev et al. 2026 (CC BY) |
| `ImrichEtAl1984.dist` | a copy of the file at the top level | Imrich et al. 1984 |
| `manhattan_distances_20x20.dist` | a copy of the file at the top level | ours |
| `dusky_dolphins.nex` | the TAXA and CHARACTERS blocks of SplitsTree's example `dusky_dolphins.nex` | Cassens et al. 2003 |
| `CCapSV.nex` | `rhabdoviruses/CCapSV.fas` as a NEXUS file, with the countries of `CCapSV_traits.txt` as a TRAITS block with coordinates | Longdon et al. 2017 (CC BY) |
| `syntrichia-its-haplotypes.nex` | the input blocks (16 haplotypes, their sample counts as traits) of the Figure 1 workflow `publications/SalugaEtAl2022/ALIGNMENT_ITS_Haplotype_Network_116_SEQ-haplotypes-fig1b.stree6` | Saługa et al. 2022 (CC0) |

## `hiv-recombination-data/`

The HIV-1 reference genome alignment of 2023 of the Los Alamos HIV Sequence Database, <https://www.hiv.lanl.gov/>;
`Citation.txt` has the reference.

| file | what it is | source |
|---|---|---|
| `HIV1_REF_2023_genome_DNA.fasta` | the alignment as downloaded, 518 reference genomes | the Los Alamos HIV Sequence Database; the sequences are GenBank entries |
| `HIV1_REF_2023_genome_DNA.hamming.dist` | its Hamming distances as RazorHaploNet computes them (the header says so); the "HIV-1 reference genomes" row of Supplementary Table S2 | ours, derived |
| `19-taxa.txt`, `19-hamming.csv`, `19-pdistance.csv` | 14 of those sequences (three A1, four D, four G, three CRF19_cpx), their Hamming and p-distances: a recombination example, CRF19 being an A1/D/G recombinant | ours, derived |

## `manhattan3D/`

| file | what it is | source |
|---|---|---|
| `manhattan3D-151coordinates.csv`, `manhattan3D-151.dist` | 151 points of a 3-D lattice and their Manhattan distances, picked with the Lattice Picker (`webapp/`) | ours, generated |
| `manhattan3D-151.html`, `manhattan3D-151.png` | the picker page with that selection (a saved copy of `webapp/3dManhattanDataGenerator.html` with the points), and a picture of it | ours |

## `other/`

| file | what it is | source |
|---|---|---|
| `bees.stree6` | six honey bee species (*Apis andreniformis*, *A. mellifera*, *A. dorsata*, *A. cerana*, *A. florea*, *A. koschevnikovi*), 677 sites, with RazorHaploNet | SplitsTree's example `bees.nex`; the origin of the alignment is not recorded there |
| `colors.dist` | perceived differences between 10 colours | Helm 1959, *A multidimensional ratio scaling analysis of color relations*, Princeton University and Educational Testing Service, technical report |
| `six.dist`, `six.-int.dist`, `six.stree6` | six primates (*Lemur catta*, *Homo sapiens*, *Pan*, *Gorilla*, *Pongo*, *Hylobates*): the Hamming distances of their mitochondrial sequences, the same doubled with the taxa renamed t0–t5, and the workflow with Neighbor-Net and RazorNet; the "six points" row of Supplementary Table S2 | the first six taxa of SplitsTree's example `primates.nex` (12 primate mtDNA sequences, 898 sites) |

## `toydata/`

| file | what it is | source |
|---|---|---|
| `toydata1-noncactus-K23-24x24.dist` | a 24-point metric whose realization is not a cactus (it contains a K2,3) | MH, constructed |
| `toydata2-cactus-15x15.dist` | a 15-point cactus metric | MH, constructed |
| `toydata3-cactus-36x36.dist` | a 36-point cactus metric | MH, constructed |

## `popart-check/` (ours)

`README.md` and nine NEXUS files in the dialect PopART reads: the common input of the haplotype comparison for
the seven alignments of Table 1 (`CCapSV-popart.nex`, `DImmSV-popart.nex`, `PAegRV-popart.nex`,
`Diamond_Princess-popart.nex`, `ALIGNMENT_ITS_Haplotype_Network_116_SEQ-popart.nex`,
`ALIGNMENT_ITS_Phylogeography_127_SEQ-popart.nex`,
`ALIGNMENT_Maritime_Antarctic_Population_Genetics_49_SEQ-popart.nex`) and the two ITS alignments without their
0/1 indel columns (`…-dna-only-popart.nex`). Written by `PopartCheck` from the data under `publications/` and
`rhabdoviruses/`; the README has the expected numbers and the result of the check.

## `rhabdoviruses/`

Longdon et al. 2017, *Vertically transmitted rhabdoviruses are found across three insect families and have
dynamic interactions with their hosts*, Proc. R. Soc. B 284:20162381. The alignments are the authors' Figshare
deposit <https://doi.org/10.6084/m9.figshare.3409420> (CC BY); the sequences are GenBank KX352837–KX353310.
`README.txt` has the links.

| file | what it is | source |
|---|---|---|
| `CCapSV.fas`, `DImmSV.fas`, `PAegRV.fas` | the three alignments of the paper: Ceratitis capitata sigmavirus (19 sequences), Drosophila immigrans sigmavirus (87), Pararge aegeria rhabdovirus (131), 1277 to 1281 sites | the authors' deposit |
| `CCapSV_traits.txt`, `DImmSV_traits.txt`, `PAegRV_traits.txt` | one line per sequence: country and coordinates (country centroids), for the trait pie charts and the map | ours, from the paper's sampling information |
| `CCapSV.stree6`, `DImmSV.stree6`, `PAegRV.stree6`, `PAegRV-haplotypes.stree6` | workflows: Median-Joining and RazorHaploNet on each alignment; `PAegRV-haplotypes` with the Shave filter and the traits | ours, derived |

## `Issakka_ITSClustering/`

Estensmo et al. 2021, *The influence of intraspecific sequence variation during DNA metabarcoding: a case study
of eleven fungal species*, Molecular Ecology Resources 21:1141–1152. The files are the `Final_alignments`
folder of `ITSClustering.zip` in the Dryad deposit <https://doi.org/10.5061/dryad.h18931zjq> (CC0);
`Dryad_ReadMe_Issakka_ITSClustering.rtf` is the deposit's README.

| file | what it is | source |
|---|---|---|
| `Final_alignments/<species>_h_seq_Nucleotide alignment.fasta` (Amylap, Antser, Fompin, Phefer, Phenig, Phevit, Poscae, Triabi) and `<species>_h_asv_Nucleotide alignment.fasta` (Fomros, Glosep, Phecen) | the authors' ITS2 haplotype alignments of eleven fungal species, from Sanger sequences (`seq`) or amplicon variants (`asv`); Amylap has a single haplotype after complete deletion and is not in the paper | the deposit |
| `Final_alignments/<species>…alignment.hamming.dist` (eleven files) | their Hamming distances, written by `CharactersToDistances` | ours, derived |
| `Final_alignments/Haplotype_network/<species>_…_varsites.nex` (eleven files) | the authors' variable-site matrices for their haplotype networks, written by EMBOSS (`Haplotype_network.zip` of the deposit) | the deposit |
| `Final_alignments/Haplotype_network/Antser_h_seq_Nucleotide_alignment_varsites.dist`, `.gml` | Hamming distances of the Antser variable sites and RazorNet's network | ours, derived |

## `publications/`

Each directory holds a `Citation.txt` or `README.md` with the paper and the deposit.

| directory and files | what it is | source, licence |
|---|---|---|
| `AmishaEtAl2026/Citation.txt`, `indian-camels-fst.nex` | pairwise F_ST among 12 camel populations, from Additional file 1 of the paper | Amisha et al. 2026, BMC Genomics 27:784 (open access, CC BY); the table transcribed by us |
| `CassensEtAl2003/Citation.txt`, `dusky_dolphins.stree6` | 36 mitochondrial haplotypes of the dusky dolphin over 60 variable sites, with p-distances, Median-Joining, Neighbor-Net and RazorHaploNet | Cassens et al. 2003, Molecular Ecology 12:1781–1792, as distributed with SplitsTree (`dusky_dolphins.nex`) |
| `DotsevEtAl2025/Citation.txt`, `kerman-wild-sheep-fst.nex` | pairwise F_ST among 8 wild sheep populations, Supplementary Table S1 | Dotsev et al. 2025, Animals 15:238 (open access, CC BY) |
| `DotsevEtAl2026/Citation.txt`, `snow-sheep-fst.nex`, `snow-sheep-fst-labels.txt`, `snow-sheep-fst.stree6` | pairwise F_ST among 14 snow sheep populations (Table S1 of the paper's supplement), display labels with the paper's colours, and the workflow of Figure 3 (Neighbor-Net, RazorNet, Shave filter) | Dotsev et al. 2026, Ecology and Evolution 16:e74182 (open access, CC BY) |
| `Greenhill2015/Citation.txt`, `Huon-languages.stree6`, `Greenhill2015.stree6` | 391 cognate sets of 14 Huon Peninsula languages with Hamming distances and Neighbor-Net; the distances with Neighbor-Net, RazorNet and the Shave filter (Table 2) | Greenhill 2015, PLoS ONE 10:e0141563; the data from TransNewGuinea.org as distributed with SplitsTree |
| `PadillaIglesiasEtAl2024/Citation.txt`, `hunter-gatherers-genetics-fst.nex`, `-music-jaccard.nex`, `-tools-jaccard.nex`, and the three `.stree6` | F_ST and Jaccard distances among ten Central African hunter-gatherer groups, computed by us from the authors' repository as `Citation.txt` describes, and the workflows | Padilla-Iglesias et al. 2024, Nature Human Behaviour 8:1263–1275; data and code <https://github.com/ceciliapad/GenCultEvo>, which states no licence; the matrices are our own computations from its tables |
| `SalugaEtAl2022/README.md`, `doi_10_5061_dryad_h18931zp3__v20220812/` (20 alignments `.fas`/`.fasta`, their `README.txt`, and our 20 `.hamming.dist`), `ALIGNMENT_ITS_Haplotype_Network_116_SEQ.stree6`, `…-Figure1a-PopART-TCS.stree6`, `…-haplotypes-fig1b.stree6`, `…-haplotypes-fig1c.stree6`, `ALIGNMENT_cpDNA_Haplotype_Network_111_SEQ.stree6`, `…-popart.nex`, `make-popart.py`, `traits.txt`, `display-labels.txt` | the Dryad deposit of the paper (twelve alignments of *Syntrichia* markers, several in two copies with different names), the workflows of Figure 1 and of the cpDNA network, the PopART file and the script that writes it, the traits and display labels | Saługa et al. 2022, Journal of Biogeography 49:1991–2004; Dryad <https://doi.org/10.5061/dryad.h18931zp3> (CC0) |
| `XuEtAl2026/README.md`, `Diamond_Princess.phy`, `.tsv`, `.vcf` (as received), `Diamond_Princess.nex`, `-traits-matrix.tsv`, `-popart.nex`, `make-traits.py`, and the four `.stree6` | 482 SARS-CoV-2 genomes of the Diamond Princess outbreak and contemporaneous genomes over 444 variable sites, with sampling locations and dates; our NEXUS, traits and PopART versions, the script that writes them, and the workflows of Figure 2 | Xu et al. 2026, Molecular Biology and Evolution 43:msag052; the data from the supplementary materials of HaploThread on Zenodo, <https://doi.org/10.5281/zenodo.18194187>, which states no licence; used with the full citation |
| `deJongEtAl2023/Citation.txt`, `brown-bear-populations-d.nex`, `brown-bear-populations-d-no-outgroups.nex` | genome-wide distances among 27 bear populations, recomputed by us from the deposit as `Citation.txt` describes (not used in the paper) | de Jong et al. 2023, Communications Biology 6:153; Dryad <https://doi.org/10.5061/dryad.qbzkh18n6> (CC0) |
