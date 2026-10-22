# The RazorNet algorithm

Momoko Hayamizu and Daniel H. Huson · draft of 2026-09-23 (written with an AI agent; unreviewed by MH)

This document describes the algorithm implemented in this repository: what it computes, the two strategies it
runs, the guarantees that can be proved for them, how the two are combined, and what has been measured. It is
meant to be read on its own.

**Contents.** §1 the problem · §2 definitions · §3 the pipeline · §4 the Strict strategy · §5 what can be
proved · §6 the Relaxed strategy · §7 the combination · §8 polishing, edge removal · §9 the haplotype
network · Appendix A measurements · Appendix B symbols ↔ code · references.

---

## 1 · The problem

A **metric space** `(X, d)`: `X` finite, `d(x,y) = d(y,x) ≥ 0`, `d(x,y) = 0 ⟺ x = y`, triangle inequality. A
**realization** of `(X, d)` is a graph `G` with positive edge weights and `X ⊆ V(G)` whose shortest-path metric
`d_G` satisfies `d_G(x,y) = d(x,y)` for all `x, y ∈ X`. Vertices of `V(G) ∖ X` are **auxiliary** points. The
**length** of `G` is the sum of its edge weights. A realization is **minimal** if no edge can be removed without
destroying the property, and **optimal** if its length is minimum among all realizations; `opt(X, d)` denotes
that minimum. Optimal realizations exist (Imrich, Simões-Pereira & Zamfirescu 1984); finding one is NP-hard
(Winkler 1988). The **fully labelled** problem — a minimal realization with `V(G) = X` — is easy (Lemma 4).

RazorNet returns an exact realization whose length is as small as it can find. It is a heuristic: the answer is
an upper bound on `opt`, never a claim to be it. What is guaranteed is exactness (Theorem 7) and, for the
Strict strategy, that every point it creates at the top level belongs to every optimal realization
(Theorem 6).

## 2 · Definitions

**Betweenness and slack.** For a subset `S` of the current point set and `x ∈ S` with `|S| ≥ 3`,

```
c_S(x) := min { ( D(x,y) + D(x,z) − D(y,z) ) / 2 : y, z ∈ S ∖ {x}, y ≠ z }         (c_S(x) := 0 if |S| < 3)
```

A pair `{y, z}` attaining the minimum is a **minimizer pair** of `x`. `x` is **between** `y` and `z` if
`D(x,y) + D(x,z) = D(y,z)`. A point with `c_S(x) > 0` is a **slack point** of `S`; `c_S(x)` is its **slack**. By
the triangle inequality `c_S(x) ≥ 0`, and `c_S(x) = 0` iff `x` is between two points of `S`.

The Relaxed strategy uses a second quantity, the **positive slack**

```
c⁺_S(x) := min { ( D(x,y) + D(x,z) − D(y,z) ) / 2 > 0 : y, z ∈ S ∖ {x}, y ≠ z }    (0 if no pair gives a positive value)
```

which ignores the pairs `x` lies between. `c⁺_S(x) ≥ c_S(x)`, with equality when `c_S(x) > 0`.

**Compactification.** For a slack point `x` with slack `c` and minimizer pair `{y, z}`, add a new point `x'` with
the distance column

```
maxAllColumn(D, x, c, y, z):
    col[x] ← c;  col[y] ← max(D(y,x) − c, 0);  col[z] ← max(D(z,x) − c, 0);  T ← {x, y, z}
    for each remaining point a of D, in increasing index:
        col[a] ← max( max{ D(a,s) − col[s] : s ∈ T }, 0 );  T ← T ∪ {a}
```

If some entry is 0 the new point would coincide with an existing one, and it is not added. Otherwise `x'` is the
median of `x, y, z` and `D(x, x') = c` (Lemma 2).

**Essential edges and realization.** For a metric `D'` on a point set `V`, the graph `G(D')` has vertex set `V`
and an edge `{u, v}` of weight `D'(u,v)` iff no third point `w` with `D'(u,w), D'(w,v) > 0` satisfies
`D'(u,w) + D'(w,v) ≤ D'(u,v)`. `G(D')` realizes `D'` (Lemma 4).

**Quantization.** The input is scaled by a power of ten so that its smallest and largest entries are integer-like,
and then doubled, so all off-diagonal entries are positive even integers; coincident taxa are merged. All
arithmetic is then exact (Lemma 3), and lengths are mapped back at the end.

## 3 · The pipeline

```
quantize the input (optionally after preprocessing)
record the taxa-only realization G(D) as the floor                         // Lemma 4, Theorem 7

Both only: count the interior taxa; the relaxed strategy runs first iff at least 10 % are interior   // §7
for each strategy to run, in that order:
    for each of N random restarts (in parallel; the strict strategy in batches with an adaptive stop):
        D ← the quantized matrix
        repeat R times (default R = 1):
            extend D                                                       // §4 or §6: the rule's core
            G ← G(D); record G if shorter than the best so far
            polish G (and extend D accordingly)                            // §8
        greedy global edge removal on G; record
    finishing polish of this strategy's own best graph with its own core; record
report: the order and why; restarts used, best length and seconds per strategy; the winner
```

"Record" means: if the current graph is shorter than the best so far, keep a copy. Every recorded graph
realizes the input, so a run interrupted by the time budget still returns a valid network.

## 4 · The Strict strategy

Package `razornet.algorithms.strict`. At every level of the recursion **all** slack points of the current
subset are compactified in one pass, against the matrix as it grows, and only the **new** points are handed to
the next level; the slack is the true minimum `c_S`.

```
procedure extend(D):  extendSubset(D, {0,…,n−1})

procedure extendSubset(D, S):
    if |D| ≥ extensionCap: return                                     // stop 0 (cap = max(5n, 500))

    slacklist ← [ (x, c_S(x), y, z) : x ∈ S, c_S(x) > 0 ]              // {y,z}: one minimizer pair, chosen
    if slacklist is empty: return                                     //   uniformly among those attaining c_S(x)

    n_before ← |D|
    for each (x, c, y, z) in slacklist, in decreasing c (ties: a random order):
        col ← maxAllColumn(D, x, c, y, z)         // against the CURRENT D, which grows in this loop
        if col has a zero entry: continue          // x' coincides with an existing point: do not add it
        append col to D as the new point x'
    if |D| = n_before: return                                         // stop 1'

    V(G) ← ( S ∖ {x : x in slacklist} ) ∪ { n_before, …, |D|−1 }      // every slack point retires; every new
    E(G) ← { (a,b) : a,b ∈ V(G), a < b, no w ∉ {a,b} with D(a,b) ≥ D(a,w) + D(w,b) }   //   point joins
    settled     ← { v ∈ V(G) : deg_G(v) ≤ 2 }
    unconfirmed ← { (u,v) ∈ E(G) : u, v ∉ settled }
    if unconfirmed is empty: return                                   // stop 2

    for each connected component C of (V(G), unconfirmed):
        extendSubset(D, V(C))
```

An equivalent formulation keeps a set `E*` instead: attach each slack `x` as a pendant vertex of its `x'`,
remove redundant edges, remove degree-1 vertices (saving their edges in `E*`), remove degree-≤2 vertices
(saving their edges), and recurse on the components that still have edges. It is the same algorithm.

**Randomness.** Only two tie-breaks: the processing order among equal slacks and the choice of the minimizer pair.
Neither changes *which* points are compactified. On instances without ties the restarts are identical.

## 5 · What can be proved

Throughout, `D` is a metric on a finite point set `P` (the taxa and the auxiliary points added so far),
`S ⊆ P`, `x ∈ S`, `c = c_S(x) > 0` with minimizer pair `{y, z}`.

**Lemma 1 (slack and betweenness).** `c_S(x) ≥ 0`; `c_S(x) = 0` iff `x` is between two points of `S ∖ {x}`;
and `c_S(x) ≤ D(x,y)` for every `y ∈ S ∖ {x}`.

*Proof.* The first two are the triangle inequality. For the bound, any pair `{y, z}` gives
`(D(x,y) + D(x,z) − D(y,z))/2 ≤ (D(x,y) + D(x,y))/2`, since `D(x,z) − D(y,z) ≤ D(x,y)`. ∎

**Lemma 2 (the compactification column).** Let `col = maxAllColumn(D, x, c, y, z)`.
(a) `col[x] = c`, `col[x] + col[y] = D(x,y)`, `col[x] + col[z] = D(x,z)`, `col[y] + col[z] = D(y,z)`: `x'` is the
**median** of `x, y, z`.
(b) If `col` has no zero entry, `D` extended by `col` is a metric on `P ∪ {x'}` whose restriction to `P` is `D`.
(c) If `S = P` (the situation at the top level), `col[a] = D(a,x) − c` for every `a ≠ x`: `x'` lies on a geodesic
from `x` to every other point, at distance `c` from `x`.

*Proof.* (a) is direct from the first three assignments and `2c = D(x,y) + D(x,z) − D(y,z)`; the values are
non-negative by Lemma 1. (b) Write `col[p] = D(x',p)`. Every placed point `s` has a witness `s* ≠ s` in `T` with
`D(x',s) = D(s,s*) − D(x',s*)` (for `x, y, z` by (a); for a later `a`, the maximizing `s`). Let `a` be the point
being placed and `s ∈ T`. `D(a,s) ≤ D(x',a) + D(x',s)` by the definition of the maximum. With `s*` the maximizer
for `a`: `D(x',a) = D(a,s*) − D(x',s*) ≤ D(a,s) + D(s,s*) − D(x',s*) ≤ D(a,s) + D(x',s)`, using the first
inequality for the pair `s, s*` (induction). With `s*` the witness of `s`:
`D(x',s) = D(s,s*) − D(x',s*) ≤ D(a,s) + D(a,s*) − D(x',s*) ≤ D(a,s) + D(x',a)`. Positivity is the hypothesis.
(c) For `a ∈ P ∖ {x,y,z}` the candidates are `D(a,x) − c` (from `s = x`), `D(a,y) − D(x,y) + c` (from `s = y`),
likewise for `z`, and `D(a,a') − D(a',x) + c` for earlier `a'` with `col[a'] = D(a',x) − c`. Each of the last
three is `≤ D(a,x) − c` iff `2c ≤ D(x,a) + D(x,·) − D(a,·)`, i.e. iff `c` is at most the slack value of the pair
`{a, ·} ⊆ S ∖ {x}` — true by minimality of `c` over `S = P`. ∎

*Remark.* Below the top level `S ⊊ P` and (c) need not hold: the column is then the largest metric-consistent
set of distances, so `x'` is pushed as far from `x` as the whole matrix allows. The maximum must run over all
placed points: a rule that uses only `x, y, z` can produce inconsistent distances.

**Lemma 3 (integrality).** Label every taxon with `π = 0` and every new point with `π(x') = c + π(x) mod 2`. With
even-integer input, `D(p,q) ≡ π(p) + π(q) (mod 2)` at all times, every slack numerator is even, and every slack
and column entry is an integer; no rounding ever occurs.

*Proof.* The numerator is `≡ 2π(x) ≡ 0`, so `c` is an integer and `c ≡ π(x') + π(x)`; then
`col[y] = D(x,y) − c ≡ π(y) + π(x')` and `col[a] = D(a,s) − col[s] ≡ π(a) + π(x')`. ∎

**Lemma 4 (realization).** For a metric `D'` on `V` without coincident points, `G(D')` realizes `D'`, every edge
of `G(D')` is essential, and `G(D')` is the unique minimal realization with vertex set exactly `V`.

*Proof.* Every path is at least as long as the distance of its endpoints, so `d_G ≥ D'`. Conversely, by
induction on `D'(u,v)`: an edge gives `d_G(u,v) ≤ D'(u,v)`; otherwise some `w` is strictly between with shorter
legs, and `d_G(u,v) ≤ D'(u,w) + D'(w,v) = D'(u,v)`. A realization on `V` lacking an essential pair connects it by
a path through a `w` that is not between, of length `> D'(u,v)`; a minimal one cannot contain a non-essential
pair. ∎

**Lemma 5 (the pendant-edge lemma).** Let `c = c_X(x) > 0`, the slack with respect to the whole set. Then in
every optimal realization `G` of `(X, d)`: (i) `x` has degree 1; (ii) after smoothing auxiliary vertices of
degree 2 (which changes neither length nor distances), its pendant edge has weight exactly `c`, and the other
endpoint `u` satisfies `d_G(u,y) = d(x,y) − c` for all `y ≠ x` — `u` realizes the compactification point `x'`.

*Proof.* (i) If `x` had edges `e_1 … e_k`, `k ≥ 2`, choose `0 < δ < min w_i`, `δ ≤ c`: shorten every `e_i` by `δ`,
add a pendant vertex `x*` at distance `δ` from `x`, move the label `x` to `x*`. The length drops by `(k−1)δ`;
`d(x, ·)` is unchanged; a path between other taxa through the old `x` lost `2δ` but was at least
`d(x,y) + d(x,z) ≥ d(y,z) + 2c` long. Contradiction. (ii) For the minimizer pair,
`2c = 2w + d_G(u,y) + d_G(u,z) − d_G(y,z) ≥ 2w`, so `w ≤ c`. If `w < c`, then `u` is on no geodesic between
taxa other than `x` and has neighbours `x, v_1 … v_k`, `k ≥ 2` after smoothing (a taxon `u` would force `w = c`
by Lemma 1); lengthening `xu` by `δ ≤ c − w` and shortening every `u v_i` by `δ` gives a shorter realization.
Hence `w = c` and `d_G(u,y) = d(x,y) − c`. ∎

**Theorem 6 (the top-level pass is lossless).** Let `D⁺` be the matrix after the top-level pass of the
Strict strategy (`S = X`). Then `opt(D⁺) = opt(D)`, and every optimal realization of `D` (with degree-2
auxiliary vertices smoothed) contains, for each slack point `x` of `X`, the pendant edge `{x, x'}` of weight
`c_X(x)`, `x'` being realized by the neighbour of `x`.

*Proof.* `opt(D⁺) ≥ opt(D)` since a realization of `D⁺` realizes `D`. Process the slack points `x_1, x_2, …` in
the algorithm's order, `c_i = c_X(x_i)`. By Lemma 2(c), `D(x_1', a) = D(x_1,a) − c_1`, so Lemma 5 applied to
`x_1` gives a realization of `D ∪ {x_1'}` of length `opt(D)`. Adding `x_1'` does not change the slack of `x_2`:
for a pair `{x_1', a}`, `(D(x_2,x_1') + D(x_2,a) − D(x_1',a))/2 = (D(x_2,x_1) + D(x_2,a) − D(x_1,a))/2 ≥ c_2`.
The algorithm's column for `x_2'` is `D(x_2,a) − c_2` on `X` (Lemma 2(c)) and
`D(x_2', x_1') = D(x_1,x_2) − c_1 − c_2` (the `s = x_2` term wins the maximum, every other term being `≤` by
minimality of `c_2`) — exactly the compactification point of Lemma 5 for `D ∪ {x_1'}`; the argument iterates,
and pairs of two new points never lower a later slack by the same computation. A coinciding `x_j'` is a point
already present, which Lemma 5(ii) identifies with the neighbour of `x_j`. ∎

**Theorem 7 (exactness, floor, termination).** For every input, RazorNet returns a graph that realizes the
input exactly on the taxa, of length at most that of the taxa-only realization `G(D)`, after finitely many
steps; each strategy's working matrix never exceeds twice its cap.

*Proof.* Extension only appends rows and columns (Lemma 2(b)); realization is exact (Lemma 4); greedy removal
and polishing keep a realization by construction, and the pipeline re-checks the result. The floor is recorded
before the restarts and only replaced by something shorter. Every call of the extension either returns or appends
at least one point and recurses only then; a call appends only if it entered below the cap, and appends at most
`|S| ≤ |D|` points. ∎

**Corollary (inheritance).** Since `Both` returns the shorter of the two strategies' results and both realize the
metric, every length statement that holds for the Strict strategy holds for `Both`: whenever Strict
returns an optimal realization, so does `Both` (Appendix A confirms `Both ≤ min` restart for restart).

**Statements not proved here.** The Strict strategy is designed to return an optimal realization on tree
metrics, cactus metrics and possibly circular metrics. The proofs are not in this document.
What is known here: by Theorem 6 the top-level pass is the classical leaf-peeling step and preserves the optimum;
the **extension alone does not always finish the tree** — six leaves in three cherries around an absent centre
give a triangle of three degree-2 vertices, all settled, of length 12 against the tree's 9 — and the polishing
step (§8) repairs it. So the tree statement is about extension plus polishing. On the test tiers all strategies
reach the known optimum on every tree instance (Appendix A).

**What is not guaranteed, for any strategy:** minimality of the output; that the recursion below the top level
preserves `opt` (the slack is local to a component and the settled rule may stop early); that polishing reaches
an optimal fixed point; independence from the tie-breaks; reproducibility under a wall-clock budget.

## 6 · The Relaxed strategy

Package `razornet.algorithms.relaxed`. The same recursion, with four changes:

- **positive slack.** `c⁺_S` replaces `c_S`: a point that lies between some pairs but not others is compactified
  too, towards the pair it is closest to lying between. Lemma 2(b) and Lemma 3 still hold (they do not use
  minimality); Lemma 5 and Theorem 6 do **not** apply to points with `c_S(x) = 0` — for those, `x'` is a
  candidate auxiliary point, not a proven one, and adding a point no optimal realization contains raises `opt` of
  the extended matrix. The heuristic therefore has to be allowed to fail, inside the best-of-N search with the
  floor of Theorem 7 as its safety net.
- **one compactification per call.** Only the point of largest `c⁺_S` whose column has no zero entry is
  compactified; the others get their turn through the recursion, each call of which draws a **fresh random
  permutation** of its subset for tie-breaking. Ties are common (on the manhattan example 40–46 % of calls have
  11–13 candidates tied at the top), so different restarts follow different paths.
- **the subset stays in play.** `V(G) ← (S ∖ {x}) ∪ {x'}`: only the compactified point leaves.
- **slack-aware settling.** `settled ← { v : deg_G(v) ≤ 2 and v is between two points of V(G) }` (the true zero
  test); a degree-2 corner that still has slack is kept for further compactification.

The extension cap is `max(20n, 2000)`.

**Why it is needed.** On many metrics almost no taxon is a slack point in the sense of §2. In the manhattan
example (twenty points of an integer grid under the L1 distance) 18 of the 20 taxa lie between some other pair;
only two have positive `c_X`. The Strict strategy therefore adds two points per restart, its recursion
stops at once, and every restart returns the same network (232 before the finishing polish, 222 after). The
optimal realization of a grid metric is a rectilinear Steiner-type network full of auxiliary corner points that
are medians of *non-pendant* taxa, and only the positive slack creates them: the Relaxed strategy reaches
191–192.

**What it costs.** Order-dependence; no optimality statement beyond Theorem 7; and on generic metrics a far
larger search (on `general-20` it is 35× slower than the Strict strategy and 4 % longer).

## 7 · Strategy `Both`: one algorithm, two strategies, the likely winner first

RazorNet is one algorithm with two strategies for the extension step: the **strict strategy** of §4 (the `Strict`
setting) and the **relaxed strategy** of §6 (the `Relaxed` setting).
`Both`, the default, runs both and keeps the shorter network; the order is decided from the data.

1. **The likely winner first.** The fraction `f` of taxa that lie between two others is counted in the quantized
   matrix (O(n³), early exit). Such *interior* taxa have zero slack under the strict strategy and keep it in every
   extension (a point between two others stays between them when points are added), so the strict strategy can never
   compactify them at the top level; only the relaxed strategy creates the medians they need. If `f ≥ 0.1` the
   relaxed strategy runs first, otherwise the strict strategy. On the real matrices measured on 2026-09-24 the relaxed
   rule, where it wins, wins by 7–47 % (PAegRV `f` = 0.37: 2266 vs 4265; cpDNA-141 `f` = 0.30: 282.5 vs 389;
   ITS-127 `f` = 0.20: 895.5 vs 1022; manhattan `f` = 0.90: 192 vs 222), whereas the strict strategy, where it wins
   on data above the threshold (general-20 `f` = 0.80, cactus-20 `f` = 0.50), wins by at most 4 % and still runs,
   second. Tree metrics, generic real-valued matrices (where the strict strategy is shorter by factors and much
   faster) and the 517 HIV-1 reference genomes have `f` = 0. Exceptions exist (colours, six points and circle-20
   have `f` = 0 and the relaxed strategy wins by 6–25 %, through its other mechanics), and the second strategy covers them
   when the run is restart-bound.
2. **Phases with their own finishing polish.** A phase is a rule's restarts followed by the finishing polish of
   that strategy's best graph with that strategy's core. The strict strategy's restarts run in parallel batches of
   `MaxNonImprovingRestarts` (default 3) and stop when a whole batch fails to shorten the best network; deciding
   on whole batches, with restart seeds fixed by index, keeps the set of executed restarts — and so the result —
   independent of thread timing. The relaxed strategy never stops early: its value is the diversity of many restarts.
3. **Budget shares.** Under a wall-clock budget the search is 95 % of `MaxSeconds`; the first phase of `Both` may
   use four fifths of it and the second gets the rest, and within a phase the restarts may use four fifths of the
   window and the finishing polish what remains (all of it when the restarts stopped early). The second strategy is
   a hedge against a wrong choice — cheap when it is the strict strategy on data of moderate size, useless when it is
   the relaxed strategy on a generic matrix — and the fifth of the budget is its price.
4. **Same seeds per phase.** Each phase draws its restart seeds from `new Random(seed)`, so a strategy inside `Both`
   executes exactly the restarts it executes alone, in either order.
5. **Both = min, restart-bound.** With its own finishing polish per strategy, `Both` is exactly the minimum of the two
   rules run alone (Appendix A: 15 of 15 instances). Under a budget it
   pays for the hedge: PAegRV at 60 s on all cores, `Both` 2034 (27 relaxed restarts) = `Relaxed` alone 2034
   (31 restarts).
6. **Provenance.** The report names the interior count and the order, the restarts used, the best length and the
   seconds per strategy, the effect of the finishing polish, and the winner (or a tie), on stderr and in the analysis
   log.

The result: the guarantees of the strict strategy are inherited (Corollary in §5), the reach of the relaxed strategy is
kept where the data call for it, and the cost of the combination is the second strategy's share of the budget.

## 8 · Polishing and edge removal

**Polishing** (`PolishGraph`, one copy per strategy, differing only in the slack they resolve): for every
vertex `v` of degree at least 3, take `N_v = {v} ∪ neighbours of v`, compactify the slack points of the metric
restricted to `N_v`, realize the extended local matrix, replace `G[N_v]` by it, and extend the global matrix
with the new points (their distances to points outside `N_v` being the minimum over `N_v` of local plus global
distance). Vertices are visited in decreasing local slack. The pipeline keeps a polished graph only if it is
shorter. Auxiliary vertices of degree 2 are smoothed away; in an optimal realization no auxiliary vertex has
degree ≤ 2.

**Greedy global edge removal**: an edge is removed iff its endpoints keep an equal-length alternative path;
every pairwise distance is preserved and the result is a minimal realization. Applied once per restart and
inside each finishing-polish cycle.

Giving up exactness for a sparser network is not part of the algorithm. It is done afterwards, by one of two
filters of SplitsTree applied to the network: the Stretch filter, which greedily removes edges, heaviest first,
while every taxon pair stays within a factor `(1 + p/100)` of its distance, and the Shave filter, which
removes edges in order of increasing cost and stops at a chosen target.

## 9 · The haplotype network

`RazorHaploNet` turns an alignment into a network: columns carrying missing or ambiguous data are
masked out (a gap is an indel and counts as a state), so the Hamming distances on the remaining columns are a
true metric; RazorNet realizes them (with the chosen strategy); the sequences of the auxiliary vertices are
inferred by parsimony, weighting each edge by its geodesic betweenness so that changes are kept off high-traffic
edges (`MinExcess`) or uniformly (`Parsimony`); from then on an edge's weight *is* its number of mutations;
superfluous edges (equal-length alternative path) are removed, unlabeled edges contracted. The TCS variant
swaps the engine and keeps everything else, so the two methods can be compared on one alignment under the same
labeling.

---

## Appendix A · Measurements (2026-09-23, this implementation)

Run through `RunRazorNet`, restart-bound (`-t 0`) and single-threaded (`-th 1`), seed 0 unless stated. Every
network realizes its input exactly.

**Manhattan example** (`examples/manhattan_distances_20x20.dist`; 20 grid points, L1):

| strategy | restarts run | length before / after finishing polish | time | winner |
|---|---|---|---|---|
| Relaxed | 100 | 192 / 192 | 2.4 s | — |
| Strict, no early stop | 100 | 232 / 222 | 1.9 s | — |
| Strict, default | 6 (stopped) | 232 / 222 | 0.3 s | — |
| Both | 6 + 100 | 192 | 2.5 s | Relaxed |
| Both, seed 3 | 6 + 100 | 201 | 2.6 s | Relaxed |

The Relaxed network has length 192, 44 nodes, 61 edges and 24 auxiliary vertices (extension counters: 1672
calls, ties in 45.6 % of them, average tie size 11.3); the Strict network has length 222, 48 nodes, 67 edges and
28 auxiliary vertices (302 calls, 200 points added per 100 restarts). Over six seeds the lengths were 191–201
(Relaxed) and 222 on every seed (Strict); the per-restart lengths over 300 restarts have 47–49 distinct values
for Relaxed and a single value for Strict.

**Slack structure of the example.** True slack points: `t11` (c = 1, 27 minimizer pairs), `t19` (c = 3, 12
pairs); the other 18 taxa are between at least one pair (`t2`: one pair, `t8`: 54). Positive slack: all 20.

**A six-leaf tree** (three cherries around an absent centre, all edges 1; tree length 9), 20 restarts:

| strategy | without polishing | with polishing |
|---|---|---|
| Strict | **12** (every restart) | 9 |
| Relaxed | 9 (12 in 2 of 20 restarts) | 9 |
| Both | 9 | 9 |

**Test tiers** (`tests/`, n ≤ 20, 100 restarts). `Both` is the minimum on all 15 instances.

| instance | reference (opt.?) | Relaxed | Strict (no early stop) | Both | winner in Both |
|---|---|---|---|---|---|
| tree-n0005 | 371 (opt.) | 371 | 371 | 371 | tie |
| tree-n0010 | 898 (opt.) | 898 | 898 | 898 | tie |
| tree-n0020 | 1856 (opt.) | 1856 | 1856 | 1856 | tie |
| cactus-n0005 | 377 | 176 | 176 | 176 | tie |
| cactus-n0010 | 441 | 313 | 313 | 313 | taxa-only floor |
| cactus-n0020 | 1199 | 865 | 856 | 856 | Strict |
| circle-n0005 | 663 | 898 | 898 | 898 | tie |
| circle-n0010 | 988 | 2330 | 2309 | 2309 | Strict |
| circle-n0020 | 1851 | 4272 | 4593 | 4272 | Relaxed |
| reticulate-n0005 | 2980 | 740 | 740 | 740 | tie |
| reticulate-n0010 | 24280 | 2325 | 2326 | 2325 | Relaxed |
| reticulate-n0020 | 127360 | 3671 | 3617 | 3617 | Strict |
| general-n0005 | 300 | 206.5 | 206.5 | 206.5 | tie |
| general-n0010 | 804 | 438.5 | 461 | 438.5 | Relaxed |
| general-n0020 | 1704 | 1198 | 1148.5 | 1148.5 | Strict |

Times on `general-20`: Strict 1.1 s, Relaxed 38 s, Both 39 s.

**Sequence data** (`PAegRV.fas`, 131 taxa, 1281 sites; `RazorHaploNet`; shortfall 0 throughout).
Unbudgeted, 20 restarts: Strict 6 restarts in 53 s, best 4265; Relaxed 20 restarts in 67 s, best
1996.5, 1855.5 after finishing polish; winner Relaxed. At the default budget of 60 s with 100 restarts: `Both` gives 2056.5 (Strict: 3 restarts in 19 s,
its share used up; Relaxed: 38 s), `Relaxed` alone 2034, `Strict` alone 4265. So on a large
alignment the default pays about 1 % for keeping the Strict strategy in the run; `Strategy = Relaxed`
recovers it.

## Appendix B · Symbols ↔ code

| here                                                                  | code                                                                                       |
|-----------------------------------------------------------------------|--------------------------------------------------------------------------------------------|
| the pipeline, strategies, restart search, finishing polish, report    | `razornet.algorithms.RazorNet`                                                             |
| `Strategy`, `MaxNonImprovingRestarts`                                 | `RazorNet.Strategy`, `RazorNet.optionMaxNonImprovingRestarts`                              |
| quantize                                                              | `utils.Quantization`, `utils.AutoScaler`, `utils.MatrixPermutation`                        |
| `c_S(x)`, random minimizer pair (Strict)                              | `strict.DistanceMatrixExtension.slackWithArgmin`                                           |
| `c⁺_S(x)` (Relaxed)                                                   | `relaxed.DistanceMatrixExtension.slackWithArgmin`                                          |
| zero-slack test for settling (Relaxed)                                | `relaxed.DistanceMatrixExtension.hasZeroSlack`                                             |
| `maxAllColumn`, coincidence test                                      | `buildAuxColumnUsingMaxAll`, `findCoincidingVertex` (both cores)                           |
| `E(G)`, the "no point between" test                                   | `isEssentialEdge` (both cores)                                                             |
| the Strict pass and recursion                                         | `strict.DistanceMatrixExtension.extendSubset`                                              |
| one point per call, fresh order                                       | `relaxed.DistanceMatrixExtension.extendSubset`, `pickBestNonCoincidingSlack`               |
| realization `G(D')`, the taxa-only floor                              | `utils.GraphFromDistances.realizeFromDistances`, `realizeTaxaOnly`                         |
| polishing                                                             | `strict.PolishGraph`, `relaxed.PolishGraph`                                                |
| greedy global edge removal                                            | `utils.GreedyGlobalEdgeRemoval`                                                            |
| the analysis log (both slack tables, trace, distribution, strategies) | `utils.AnalysisLog`, `utils.ExtensionTrace`                                                |
| haplotype network                                                     | `razornet.algorithms.RazorHaploNet`, `utils.ParsimonyLabeler`, `utils.GeodesicBetweenness` |

## References

Hakimi & Yau, *Distance matrix of a graph and its realizability*, Quart. Appl. Math. 22 (1965).
Buneman, *A note on the metric properties of trees*, J. Combin. Theory B 17 (1974).
Imrich, Simões-Pereira & Zamfirescu, *On optimal embeddings of metrics in graphs*, J. Combin. Theory B 36 (1984).
Dress, *Trees, tight extensions of metric spaces, and the cohomological dimension of certain groups*, Adv. Math. 53 (1984).
Althöfer, *On optimal realizations of finite metric spaces by graphs*, Discrete Comput. Geom. 3 (1988).
Winkler, *The complexity of metric realization*, SIAM J. Discrete Math. 1 (1988).
Leigh & Bryant, *PopART: full-feature software for haplotype network construction*, Methods Ecol. Evol. 6 (2015).
Templeton, Crandall & Sing, *A cladistic analysis of phenotypic associations with haplotypes … III. Cladogram estimation*, Genetics 132 (1992).
