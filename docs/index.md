# RazorNet — User Manual

> **Draft — work in progress.**
>
> **References note.** The reference list requires work.

## 1. What is RazorNet?

**RazorNet** computes a *minimum-length network representation of a distance matrix*. Given the pairwise
distances between a set of taxa, it constructs a network — a graph with weighted edges — whose shortest-path
distances reproduce the input distances **exactly** (up to the chosen numerical precision), using as little
total edge length as possible. To do so it may introduce additional (Steiner / median) nodes wherever they
let it represent the same distances with less total length.

This is what distinguishes RazorNet from methods such as Neighbor-Net, which fit an *approximate* split
network to the distances: a RazorNet is an exact realization of the metric. A network that is too complex to
read can be simplified afterwards, with a controlled amount of distortion (see **Simplifying a network** below).

## 2. What is RazorHaploNet?

**RazorHaploNet** is the sequence-based counterpart of RazorNet. Given a multiple sequence alignment, it
computes Hamming distances between the sequences, builds the RazorNet on those distances, and then **labels
the edges with the sequence changes (mutations)** they represent. The result is a haplotype network in which
the nodes are the input sequences together with any inferred intermediate haplotypes, and edge lengths
correspond to numbers of mutations.

## 3. Relationship to SplitsTree

The RazorNet program is **built on SplitsTree** (the [SplitsTree App](https://husonlab.github.io/splitstree6)): it uses SplitsTree's data model, its
readers and writers for the common phylogenetic file formats, and its interactive views for drawing,
formatting and exporting networks.

Both **RazorNet** and **RazorHaploNet** are also available as algorithms **in SplitsTree** itself.

One consequence is that **much of the detailed SplitsTree manual applies to the RazorNet program as well** —
in particular the sections on importing and exporting data, navigating and formatting the views, selecting
taxa, and saving results. This manual describes only what is specific to RazorNet; for everything else, please
refer to the [SplitsTree documentation](https://husonlab.github.io/splitstree6).

## 4. Getting started: the examples

**File → Examples** opens one of six documents that come with the program: three distance matrices with their
RazorNet networks and three sequence alignments with their RazorHaploNet networks, computed already, each
followed by the Shave filter (see **Simplifying a network** below), so what appears is the simplified network,
with the exact network as the step before the filter in the workflow; the grid example alone keeps the filter
switched off. An example is loaded from the program's
memory and never written to disk: it opens in the current window if that is still empty, otherwise in a new
window, and Save asks where to save it. The six documents are also in the `examples/bundled` directory of the
source repository, and each carries a comment naming its source.

| File → Examples | the data | source |
|---|---|---|
| Snow sheep FST, 14 populations (RazorNet) | pairwise F_ST between 14 populations of the snow sheep, with the colours of the paper as display labels: Figure 3 of the RazorNet paper. F_ST values do not satisfy the triangle inequality; RazorNet reports this in the message window and takes the metric closure of the distances before the search | Dotsev et al. (2026), Table S1 |
| Imrich et al. 1984, 24 points (RazorNet) | a metric on 24 points from the paper that introduced optimal realizations | Imrich, Simões-Pereira and Zamfirescu (1984) |
| Manhattan grid, 20 points (RazorNet) | 20 random points on an integer grid and their Manhattan distances; the network recovers the grid exactly, so this example has the Shave filter switched off (target None) and shows the exact network (Supplementary Figure S6 of the paper) | the authors |
| Dusky dolphins mtDNA, 36 haplotypes (RazorHaploNet) | 36 mitochondrial haplotypes over 60 variable sites, the data of a classic comparison of network methods | Cassens et al. (2003), as distributed with SplitsTree |
| Rhabdovirus PAegRV, 90 haplotypes at 8 sites (RazorHaploNet) | the 131 sequences of the *Pararge aegeria* rhabdovirus, 1281 sites, collapsed to their 90 distinct haplotypes as File → Open Haplotype Data does, with the eight sampling sites (six in Britain, Corsica, Sardinia) as traits with coordinates: try the trait pie charts and the world map | Longdon et al. (2017) |
| Syntrichia ITS haplotypes with counts (RazorHaploNet) | 16 ITS haplotypes of the moss *Syntrichia*, 951 sites, and the number of samples of each as traits: Figure 1 of the RazorNet paper | Saługa, Ochyra and Ronikier (2022) |

## 5. Related methods, for comparison

So that you can compare different representations of the same data, the program also provides several
established methods:

- **Neighbor-Net** — a widely used agglomerative method that fits a *split network* to a distance
  matrix (Bryant & Moulton, 2004).
- **Median-Joining network** — a popular method for building haplotype networks from sequences
  (Bandelt, Forster & Röhl, 1999).
- **Minimum spanning network** — a haplotype network built from minimum spanning trees
  (Excoffier & Smouse, 1994).
- **TCS** — a further method for distances and for sequences, implemented after PopART (Leigh & Bryant,
  2015): TCS (Templeton, Crandall & Sing, 1992) joins haplotypes in order of increasing distance.

## 6. The RazorNet algorithm (distances → network)

RazorNet takes a distance matrix and returns a minimum-length network that realizes it. The search uses many
randomized restarts and keeps the shortest valid network found within the time budget.

| Option | Default | Meaning |
|---|---|---|
| **Preprocessing** | None | Optionally transform the distances before the search. **ClassicalMDS2D / ClassicalMDS3D** embed the distances in the plane (or in 3-space) by classical multidimensional scaling — useful when the data are close to low-dimensional, e.g. geographic distances; pair with a low **Significant digits** and a **Max stretch percent** for a sparse, map-like network. **PruneSmallSplits** denoises the metric via split decomposition, dropping the weakest splits so that a "generic" (tie-poor) metric does not cause the search to over-expand. |
| **Strategy** | Both | How the search chooses the points at which it adds auxiliary nodes. **Strict** adds them only at taxa that lie on no shortest path between two others; it is fast and nearly deterministic. **Relaxed** also adds them at taxa that lie between others, and explores the many alternatives by randomized restarts; this is what grid-like distances and alignments with sampled intermediate haplotypes need. **Both** runs both strategies, the one that is likely to do better first, and keeps the shorter network. |
| **Significant digits** | 6 | The precision at which the distances are quantized before the search. Continuous or noisy metrics usually call for a **low** value (e.g. 2–3), which keeps the search sparse. |
| **Max seconds** | 60 | Wall-clock time budget for the search; the best network found so far is returned when it is reached (`0` = no limit). |

Further options (the number of random restarts, the random seed, the polishing steps, the split filter and a
log file) are available in advanced mode; see the README of the source code.

Each drawn network reports its **total length** and its **distortion** — the worst-case multiplicative spread
of realized-versus-input pairwise distances (a distortion of 1 means an exact representation; larger values
mean the network stretches some pairs relative to others).

### Simplifying a network

A RazorNet displays the distances exactly, and for noisy distances that can be a complex network. Two filters,
applied to the network after it has been computed, trade exactness for a sparser network:

- The **Shave filter** removes edges one at a time, always the edge whose removal costs least, and stops at
  a chosen target. The default target, *Balanced*, stops where further simplification would cost more than it
  gains; other targets are a maximal deviation, a number of edges, a number of cycles (0 gives a tree) or a
  percentage of the edges. The cost is the distortion of the distances or, for a haplotype network, the excess
  of the mutations displayed over the sequence differences.
- The **Stretch filter** removes edges, heaviest first, as long as no pairwise distance grows by more than a
  given percentage.

Both report the distortion they introduce.

## 7. The RazorHaploNet algorithm (sequences → haplotype network)

RazorHaploNet builds a RazorNet from the Hamming distances between aligned sequences and then labels its
edges with the inferred sequence changes. A gap counts as a character state; columns in which a sequence has
missing or ambiguous data are left out.

In addition to **Strategy** and **Max seconds** (which behave as for RazorNet, above), it offers:

| Option | Default | Meaning |
|---|---|---|
| **Labeling** | MinExcess | How the sequence changes are placed on the edges. **MinExcess** keeps the displayed pairwise distances as faithful as possible; **Parsimony** minimizes the total number of changes. |
| **Remove superfluous** | true | Remove any edge for which an alternative path of the same length exists. |
| **Contract unlabeled** | true | Contract any edge that carries no sequence change. |
| **Missing as state** | false | Treat `?`, `N` and ambiguity codes as ordinary character states instead of leaving out their columns. Unknown positions then count as mutations; use this only for data in which too little would remain otherwise. |

## 8. References

- H.-J. Bandelt and A. W. M. Dress. A canonical decomposition theory for metrics on a finite set.
  *Advances in Mathematics*, 92:47–105, 1992.
- H.-J. Bandelt, P. Forster, and A. Röhl. Median-joining networks for inferring intraspecific phylogenies.
  *Molecular Biology and Evolution*, 16(1):37–48, 1999.
- D. Bryant and V. Moulton. Neighbor-Net: an agglomerative method for the construction of phylogenetic
  networks. *Molecular Biology and Evolution*, 21(2):255–265, 2004.
- I. Cassens, K. Van Waerebeek, P. B. Best, E. A. Crespo, J. Reyes, and M. C. Milinkovitch. The phylogeography of dusky dolphins (*Lagenorhynchus obscurus*): a critical examination of network methods and rooting procedures. *Molecular Ecology*, 12(7):1781–1792, 2003.
- A. V. Dotsev, I. M. Okhlopkov, M. G. Bondar, D. I. Litovka, T. E. Deniskova, V. R. Kharzinova, M. N. Romanov, T. P. Sipko, N. V. Mamaev, M. N. Semerikova, O. A. Koshkina, N. F. Bakoev, V. A. Bagirov, A. A. Sitsko, D. G. Medvedev, A. N. Kuksin, M. H. Moradi, A. Esmailizadeh, D. K. Griffin, H. Reyer, K. Wimmers, G. M\'esz\'aros, J. S\"olkner, I. Medugorac, G. Brem, and N. A. Zinovieva. Unraveling the evolutionary history of snow sheep (Ovis nivicola): genome-wide single nucleotide polymorphism analysis reveals genetic diversity and clarifies taxonomy. *Ecology and Evolution*, 16(8):e74182, 2026.
- L. Excoffier and P. E. Smouse. Using allele frequencies and geographic subdivision to reconstruct gene
  trees within a species: molecular variance parsimony. *Genetics*, 136(1):343–359, 1994.
- J. C. Gower. Some distance properties of latent root and vector methods used in multivariate analysis.
  *Biometrika*, 53:325–338, 1966.
- D. H. Huson and D. Bryant. Application of phylogenetic networks in evolutionary studies.
  *Molecular Biology and Evolution*, 23(2):254–267, 2006.
- D. H. Huson and D. Bryant. The SplitsTree App: interactive analysis and visualization using phylogenetic
  trees and networks. *Nature Methods*, 21(10):1773–1774, 2024.
- W. Imrich, J. M. S. Simões-Pereira, and C. M. Zamfirescu. On optimal embeddings of metrics in graphs. *Journal of Combinatorial Theory, Series B*, 36:1–15, 1984.
- A. V. Karzanov. Metrics with finite sets of primitive extensions. *Annals of Combinatorics*, 1998.
- J. W. Leigh and D. Bryant. PopART: full-feature software for haplotype network construction.
  *Methods in Ecology and Evolution*, 6(9):1110–1116, 2015.
- B. Longdon, J. P. Day, N. Schulz, P. T. Leftwich, M. A. de Jong, C. J. Breuker, M. Gibbs, D. J. Obbard, L. Wilfert, S. C. L. Smith, J. E. McGonigle, T. M. Houslay, L. I. Wright, L. Livraghi, L. C. Evans, L. A. Friend, T. Chapman, J. Vontas, N. Kambouraki, and F. M. Jiggins. Vertically transmitted rhabdoviruses are found across three insect families and have dynamic interactions with their hosts. *Proceedings of the Royal Society B: Biological Sciences*, 284:20162381, 2017.
- M. Saługa, R. Ochyra, and M. Ronikier. Phylogeographical breaks and limited connectivity among multiple refugia in a pan-Antarctic moss species. *Journal of Biogeography*, 49(11):1991–2004, 2022.
- A. R. Templeton, K. A. Crandall, and C. F. Sing. A cladistic analysis of phenotypic associations with
  haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram estimation.
  *Genetics*, 132(2):619–633, 1992.
- M. Hayamizu and others and D.H. Huson. RazorNet: minimum-length network representations of metrics.
  *In preparation.*
