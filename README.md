# RazorNet

**RazorNet** computes and visualizes minimum-length networks that realize a metric exactly: given distances, or
sequences from which Hamming distances are computed, it builds a weighted graph, using auxiliary vertices where
they help, whose shortest-path distances between the taxa reproduce the input and whose total edge length is as
small as the search can make it. Momoko Hayamizu and Daniel H. Huson, 2026. GPL v3, see `LICENSE.txt`.

This repository holds:

1. **The algorithms** — module `algorithms/`, Maven artifact `razornet-algorithms`, Java module
   `razornet.algorithms`: `RazorNet` (distances → network), `RazorHaploNet` (sequences → haplotype
   network), and, for comparison on the same data under the same labeling, **TCS** with its haplotype variant.
   They register with SplitsTree6 through its `AlgorithmProvider` service, on the module path via
   `module-info.java` and on the class path via `META-INF/services`. When RazorNet is released, SplitsTree6 will
   offer them too; whether they move into the splitstree6 code base or SplitsTree6 uses `razornet-algorithms` as
   an external dependency is not decided yet.
2. **The stand-alone RazorNet program** — module `app/`, artifact `razornet-app`: a JavaFX application built on
   SplitsTree6 through its `IAppProfile` interface (`razornet.main.RazorNet`), the command-line tool
   `razornet.tools.RunRazorNet`, and the test suite `razornet.testsuite`.

The algorithm has one extension step with two strategies, chosen by the `Strategy` option. **Strict** compactifies
exactly the taxa that lie on no geodesic between two others, all in one pass, and is provably safe at the top
level. **Relaxed** also compactifies taxa that lie between others, one per call, and explores the many ties with
random restarts, which is what grids and Hamming distances with sampled intermediate haplotypes need. The
default, **Both**, runs both strategies, the likely winner first (the relaxed strategy when at least 10 % of the
taxa lie between two others), keeps the shorter network, and reports which strategy produced it.
[docs/algorithm.md](docs/algorithm.md) is the full description, with the guarantees that can be proved and the
measurements; [docs/index.md](docs/index.md) is the user manual.

## Layout

```
algorithms/        the algorithms (see above)
app/               the program, its tools and test suite
examples/          data sets, among them manhattan_distances_20x20.dist: twenty points of a grid under the L1 distance
tests/             generated test instances with reference lengths, five tiers x five sizes, and example matrices
docs/              algorithm.md (the algorithm, with proofs and measurements), index.md (the user manual)
manuscript/        the paper for MBE (LaTeX): main.tex is the working draft, main-oup.tex the journal layout,
                   supplement.tex the proofs; see manuscript/README.md
artwork/ release-notes/ scripts/ webapp/
```

## Build and run

Maven, Java 17 sources, run on JDK 23. `jloda3` and `splitstree6` are resolved from `~/.m2`, so install them
first (jloda3 → splitstree6 → this repository):

```bash
mvn -o -q -DskipTests install
```

The program: `razornet.main.RazorNet` (JavaFX; from IntelliJ, or `mvn javafx:run` in `app/`). Headless, with
`app/target/classes` plus the class path from `mvn dependency:build-classpath` in `app/`:

```bash
java -cp ... razornet.tools.RunRazorNet -i examples/manhattan_distances_20x20.dist -strategy Both -t 0
java -cp ... razornet.testsuite.RunTestSuite --input tests --maxSize 20 --maxSeconds 0
java -cp ... razornet.testsuite.CompareStrategies --input tests/tree-distances --maxSize 20
java -cp ... razornet.testsuite.RunHaplotypeNet examples/rhabdoviruses/PAegRV.fas --strategy Both -t 60
java -cp ... razornet.testsuite.CompareDistanceMethods examples/manhattan_distances_20x20.dist
java -cp ... razornet.testsuite.CompareHaplotypeMethods examples/rhabdoviruses/CCapSV.fas
```

`CompareDistanceMethods` compares NeighborNet with RazorNet (`CharactersToDistances` turns an alignment into the Hamming distances of its distinct sequences for it), and `CompareHaplotypeMethods` median joining and TCS
with RazorHaploNet, each also followed by the Shave filter. `SimulateTreeMixtures` and `SimulateHaplotypes`
generate data sets with a known truth for them (sequences are evolved with AliSim of IQ-TREE 2, which has to be
installed), and `--truth` then reports how the quartet supports displayed by a network compare with the true ones. `PopartCheck` writes the common input of the haplotype comparison in the NEXUS dialect that
PopART reads and prints the size of our TCS construction on it, for a check against PopART (`examples/popart-check/`).

`CompareStrategies` is restart-bound by default; point it at the generated tiers (`tests/<tier>-distances`)
rather than at `tests`, since `tests/example-distances` holds matrices on which a run without a time budget takes
very long.

Installers are built by `build-installers` (`apps/razornet/config.json`) and published as releases of this
repository, which is where the program's *Check for Updates* looks (`razornet.main.Version.HOME_URL`).

## Options of RazorNet

The algorithm pane shows the **basic** options. The **advanced** options appear after a separator when advanced
mode is on: start the program with `-Drazornet.advanced=true`, or set `RazorNetAdvancedOptions=true` in the
program's properties file once. The command-line tools always offer every option. (A workflow saved in advanced
mode with non-default advanced values loads in basic mode with those values back at their defaults.)

| basic option | default | meaning |
|---|---|---|
| `Preprocessing` | `None` | optional distance preprocessing (classical MDS, split pruning) |
| `Strategy` | `Both` | `Strict` (strict strategy), `Relaxed` (relaxed strategy), or `Both` (both, the likely winner first; shorter wins) |
| `SignificantDigits` | 6 | precision of the integer quantization; fewer digits give a coarser, sparser network |
| `MaxSeconds` | 60 | wall-clock budget (0 = none) |

| advanced option | default | meaning |
|---|---|---|
| `RandomRestarts` | 1000 | restarts per strategy |
| `MaxNonImprovingRestarts` | 3 | strict strategy only: its restarts run in batches of this size, and the phase stops when a batch brings no improvement (0 = run all) |
| `MaxPolishRounds` | 1 | extend+polish cycles per restart |
| `FinishingPolish` | on | after the search, polish each strategy's best graph with its own core until no shorter |
| `SplitFilterPercent` | 1 | for the split-pruning preprocessing |
| `RandomSeed` | 0 | seed |
| `LogFile` | empty | write an analysis log: slack tables under both definitions, a recursion trace, the restart-length distribution, what each strategy achieved |

`RazorHaploNet` shows `Strategy`, `Labeling`, `RemoveSuperfluous`, `ContractUnlabeled`, `MissingAsState`,
`MaskGapSites` (mask every column in which any sequence has a gap, as PopART does, instead of treating a gap as a
state; off by default) and `MaxSeconds`; `RandomRestarts` and `RandomSeed` are advanced. Sequences that are
identical on the sites kept are one haplotype: the network is computed on one taxon per haplotype, and the other
sequences of a haplotype hang off its node by a zero-length edge.

Both algorithms return a network that displays the input distances exactly. A sparser network at the price of
some distortion is obtained afterwards, by applying SplitsTree's `ShaveFilter` or `StretchFilter` to the
network; the program offers both.

A run reports what each strategy achieved and which one produced the network:

```
Both: 18 of 20 taxa lie between two others (90%, threshold 10%); running the Relaxed strategy first
Relaxed:     100 restarts, best 192 (2.36 s)
Strict:      3 restarts, best 232 (0.06 s; stopped: a whole batch brought no improvement)
after finishing polish: Relaxed 192 -> 192, Strict 232 -> 222
winner: Relaxed strategy, length 192
```
