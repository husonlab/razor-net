# The experiments behind the RazorNet paper

What is needed to repeat the comparisons of the paper (Momoko Hayamizu, additional authors, and Daniel H. Huson,
*RazorNet and RazorHaploNet: haplotype and phylogenetic networks that display distances exactly*, in preparation)
and to regenerate its figures and tables. The manuscript itself is in its own repository, `razornet-paper`, which
holds copies of `results/`, `figures/` and `tables/` as they were when the text was written. The tools are in
the `app` module (`razornet.testsuite`): `CompareHaplotypeMethods`, `CompareDistanceMethods`,
`CharactersToDistances`, `CompareStrategies`, `PopartCheck`, `SimulateTreeMixtures` and `SimulateHaplotypes`;
sequences are evolved with AliSim of IQ-TREE 2 (`iqtree2` on the `PATH`).

| file | what it does |
|---|---|
| `run-published-alignments.sh` | Median-Joining, TCS, RazorHaploNet and RazorHaploNet with the Shave filter on the published alignments of `examples/`; writes `results/published-alignments.tsv` |
| `run-published-distances.sh` | Neighbor-Net, RazorNet and RazorNet with the Shave filter on published distance data: alignments and character matrices are converted to Hamming distances of their distinct sequences by `CharactersToDistances`, the distance matrices of SplitsTree's `examples/publications` are used as they are; writes `results/published-distances.tsv` |
| `run-simulation-study-1.sh` | generates the data sets of simulation study 1 (sequences evolved down gene trees that each differ from a species tree by one transfer) and runs Neighbor-Net and RazorNet on them |
| `run-simulation-study-2.sh` | generates the samples of simulation study 2 (sequences from a coalescent genealogy) and runs the haplotype network methods on them |
| `make-results.sh` | collects the reports of the two studies (`collect.py`), draws the figures and writes the summaries and the numbers quoted in the text (`figures.py`, needs matplotlib), and writes the rows of the supplementary tables (`tables.py`) |
| `results/` | the results: one row per data set and method, the summaries, and `simulation-numbers.tex`, the numbers the text quotes as macros |
| `figures/`, `tables/` | the figures of the paper and the rows of its supplementary tables; the simulation figures and the table rows are written by `make-results.sh`, the networks of the examples are exported from RazorNet (workflow files in `examples/publications/`) |
| `manhattan-points.py` | draws the twenty points of the grid example with their taxon labels (Supplementary Figure S6a) from `examples/manhattan_distances_20x20-points.tsv` |
| `concept-figure.py` | draws panels (b) to (d) of Supplementary Figure S7, the idea behind RazorNet, on `examples/five.dist` |
| `make-dryad.sh`, `dryad-README.md` | assembles the Dryad deposit of the paper in a directory given on the command line: the common inputs of the 26 alignments, the seven distance matrices as used, the results, the workflow files and exports of the figures, the PopART inputs, and the README |

The scripts take the class path of the `app` module with its dependencies from `JAVA_CP`, with absolute paths
and the two `target/classes` directories first, so that the classes just built are used and not an installed
jar:

```bash
mvn -o -q -DskipTests package
export JAVA_CP="$PWD/app/target/classes:$PWD/algorithms/target/classes:$(cd app && mvn -o -q dependency:build-classpath -Dmdep.outputFile=/dev/stdout)"
experiments/run-simulation-study-1.sh /path/to/study1      # about 2.5 hours with three processes
experiments/run-simulation-study-2.sh /path/to/study2
PYTHON=python3 experiments/make-results.sh /path/to/study1 /path/to/study2
```

The simulated data sets are not in the repository; they are generated again from the seeds in the scripts, and
a run that is interrupted continues where it stopped. `REPLICATES` and `PARALLEL` can be set in the environment.
The check of our TCS construction against PopART is in `examples/popart-check/`.
