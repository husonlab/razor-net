/*
 * GenerateTestInstances.java Copyright (C) 2026 Daniel H. Huson
 *
 *  (Some files contain contributions from other authors, who are then mentioned separately.)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

package razornet.testsuite;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

/**
 * Generates a suite of distance-matrix test instances for RazorNet: the four tiers of the classic
 * hierarchy tree &sub; cactus &sub; circular &sub; general, plus a reticulate (sum-of-two-trees) family:
 * <ul>
 *   <li><b>tree</b>     &mdash; additive metric: leaf-to-leaf path distances in a random weighted tree.
 *       The generating tree is the minimum-length realization, so its total edge length is the true
 *       optimum RazorNet should recover.</li>
 *   <li><b>cactus</b>   &mdash; cactus metric: shortest-path distances of a random edge-weighted cactus
 *       (a connected graph in which every edge lies in at most one cycle -- its blocks are bridges and
 *       edge-disjoint simple cycles). A cactus is outer-planar, so its metric is circular-decomposable
 *       (Kalmanson): cactus metrics are a subset of circular. The cactus is a valid realization, so its
 *       total edge length is an upper bound on the optimum.</li>
 *   <li><b>circle</b>   &mdash; circular-decomposable metric: {@code d(i,j) = } sum of the weights of the
 *       circular splits (intervals of a circular ordering) that separate i and j. Realizable by an
 *       outer-planar split network. The reference length is the total split weight.</li>
 *   <li><b>general</b>  &mdash; arbitrary metric: all-pairs shortest-path distances of a random connected
 *       weighted graph. That graph is a valid realization, so its total edge length is an upper bound on
 *       the optimum (RazorNet should realize the distances with length &le; this).</li>
 *   <li><b>reticulate</b> &mdash; sum of two tree metrics: {@code D = D(tree1) + D(tree2)} on the same taxa,
 *       where tree2 is tree1 after a small number of SPR (subtree prune-and-regraft) moves. This is the
 *       classic "two conflicting trees" signal that Neighbor-Net is used to display; there is no known
 *       minimum realization, so the reference length is the taxa-only (auxiliary-free) upper bound.</li>
 * </ul>
 * Output for each instance is a {@code .dist} file (a plain distance matrix that both RazorNet and
 * SplitsTree read, with a {@code #} metadata header and {@code t1..tn} row labels) plus a {@code .truth.txt}
 * describing the generating structure and its total length. Files land in
 * {@code <output>/{tree,cactus,circle,general,reticulate}-distances/}.
 * <p>
 * Usage: {@code GenerateTestInstances [--output DIR] [--sizes 5,10,20,40,80] [--reps 1] [--seed 42] [--maxWeight 100]}
 *
 * @author Daniel Huson, 2026
 */
public class GenerateTestInstances {

	/** One generated instance: the distance matrix, a reference length, whether that length is the proven optimum, and a ground-truth description. */
	private record Instance(int[][] distances, long referenceLength, boolean referenceIsOptimal, String truth) {
	}

	public static void main(String[] args) throws IOException {
		var outputDir = "tests";
		var sizes = new int[]{5, 10, 20, 40, 80};
		var reps = 1;
		var seed = 42L;
		var maxWeight = 100;
		var splitsPerTaxon = 1.0; // circle tier: number of non-trivial circular splits per taxon (sparse default)
		var density = 0.0;        // circle tier: if >0, dense mode = keep this fraction of all n(n-1)/2 splits instead
		var cactusCycleProb = 0.5;// cactus tier: probability of attaching a cycle (vs a bridge) at each growth step
		var numSprMoves = 2;      // reticulate tier: number of SPR moves turning tree1 into tree2
		var reticulationScale = 1.0; // reticulate tier: scale on tree2's edge weights (1.0 = equal-strength signals)

		for (var i = 0; i < args.length; i++) {
			try {
				switch (args[i]) {
					case "--output", "-o" -> outputDir = args[++i];
					case "--sizes", "-n" -> sizes = parseInts(args[++i]);
					case "--reps", "-r" -> reps = Integer.parseInt(args[++i]);
					case "--seed", "-s" -> seed = Long.parseLong(args[++i]);
					case "--maxWeight", "-w" -> maxWeight = Integer.parseInt(args[++i]);
					case "--splitsPerTaxon", "-k" -> splitsPerTaxon = Double.parseDouble(args[++i]);
					case "--density", "-d" -> density = Double.parseDouble(args[++i]);
					case "--cactusCycleProb", "-c" -> cactusCycleProb = Double.parseDouble(args[++i]);
					case "--sprMoves", "-p" -> numSprMoves = Integer.parseInt(args[++i]);
					case "--reticulationScale", "-x" -> reticulationScale = Double.parseDouble(args[++i]);
					case "--help", "-h" -> {
						printUsage();
						return;
					}
					default -> {
						System.err.println("Unknown argument: " + args[i]);
						printUsage();
						System.exit(1);
					}
				}
			} catch (ArrayIndexOutOfBoundsException | NumberFormatException ex) {
				System.err.println("Bad value for argument '" + args[i - 1 < 0 ? 0 : i] + "': " + ex.getMessage());
				printUsage();
				System.exit(1);
			}
		}

		System.err.printf("Generating test instances: sizes=%s reps=%d seed=%d maxWeight=%d circle=%s cactusCycleProb=%s reticulate=%dSPR,x%s -> %s%n",
				Arrays.toString(sizes), reps, seed, maxWeight,
				density > 0 ? "density=" + density : "splitsPerTaxon=" + splitsPerTaxon, cactusCycleProb, numSprMoves, reticulationScale, outputDir);

		// one master RNG drives per-instance seeds, so the whole suite is reproducible from a single seed
		var master = new Random(seed);
		var count = 0;
		for (var type : List.of("tree", "cactus", "circle", "general", "reticulate")) {
			var dir = Path.of(outputDir, type + "-distances");
			Files.createDirectories(dir);
			for (var n : sizes) {
				for (var r = 0; r < reps; r++) {
					var instanceSeed = master.nextLong();
					var rng = new Random(instanceSeed);
					var instance = switch (type) {
						case "tree" -> generateTree(n, rng, maxWeight);
						case "cactus" -> generateCactus(n, rng, maxWeight, cactusCycleProb);
						case "circle" -> generateCircular(n, rng, maxWeight, splitsPerTaxon, density);
						case "reticulate" -> generateReticulate(n, rng, maxWeight, numSprMoves, reticulationScale);
						default -> generateGeneral(n, rng, maxWeight);
					};
					var base = "%s-n%04d%s".formatted(type, n, reps > 1 ? "-r" + r : "");
					writeDist(dir.resolve(base + ".dist"), type, n, instanceSeed, instance);
					writeTruth(dir.resolve(base + ".truth.txt"), instanceSeed, instance);
					System.err.printf("  %-8s n=%-4d referenceLength=%-8d%s -> %s.dist%n",
							type, n, instance.referenceLength(), instance.referenceIsOptimal() ? " (optimal)" : "", base);
					count++;
				}
			}
		}
		System.err.printf("Done: wrote %d instances.%n", count);
	}

	// -----------------------------------------------------------------------------------------------
	// tree: additive metric from a random weighted tree (leaves = taxa)
	// -----------------------------------------------------------------------------------------------

	private static Instance generateTree(int n, Random rng, int maxWeight) {
		// Build a random tree topology by stepwise addition: repeatedly subdivide a random edge with a
		// new internal node and attach the next leaf there. Leaves 0..n-1 are the taxa.
		var edges = new ArrayList<int[]>(); // each entry {a, b}
		var nextInternal = n;
		if (n >= 2)
			edges.add(new int[]{0, 1});
		for (var k = 2; k < n; k++) {
			var e = edges.remove(rng.nextInt(edges.size()));
			var m = nextInternal++;
			edges.add(new int[]{e[0], m});
			edges.add(new int[]{m, e[1]});
			edges.add(new int[]{m, k});
		}
		var numNodes = Math.max(n, nextInternal);

		// random positive edge weights
		var weight = new int[edges.size()];
		var total = 0L;
		for (var i = 0; i < edges.size(); i++) {
			weight[i] = 1 + rng.nextInt(maxWeight);
			total += weight[i];
		}

		// adjacency (neighbor, weight)
		var adj = new ArrayList<List<int[]>>();
		for (var i = 0; i < numNodes; i++)
			adj.add(new ArrayList<>());
		for (var i = 0; i < edges.size(); i++) {
			var e = edges.get(i);
			adj.get(e[0]).add(new int[]{e[1], weight[i]});
			adj.get(e[1]).add(new int[]{e[0], weight[i]});
		}

		// leaf-to-leaf distances = unique tree-path lengths (one DFS per leaf)
		var d = new int[n][n];
		for (var s = 0; s < n; s++) {
			var dist = new int[numNodes];
			Arrays.fill(dist, -1);
			dist[s] = 0;
			var stack = new ArrayList<Integer>();
			stack.add(s);
			while (!stack.isEmpty()) {
				int u = stack.remove(stack.size() - 1);
				for (var nb : adj.get(u)) {
					if (dist[nb[0]] < 0) {
						dist[nb[0]] = dist[u] + nb[1];
						stack.add(nb[0]);
					}
				}
			}
			for (var t = 0; t < n; t++)
				d[s][t] = dist[t];
		}

		var sb = new StringBuilder();
		sb.append("type=tree\n").append("n=").append(n).append("\n");
		sb.append("referenceLength=").append(total).append("\n");
		sb.append("referenceIsOptimal=true\n");
		sb.append("# Additive (tree) metric; the generating tree IS the minimum-length realization.\n");
		sb.append("# nodes: leaves 0..").append(n - 1).append(" = taxa t1..t").append(n)
				.append("; internal (auxiliary) nodes ").append(n).append("..").append(numNodes - 1).append("\n");
		sb.append("# edges: 'edge u v weight'\n");
		for (var i = 0; i < edges.size(); i++) {
			var e = edges.get(i);
			sb.append("edge ").append(e[0]).append(" ").append(e[1]).append(" ").append(weight[i]).append("\n");
		}
		return new Instance(d, total, true, sb.toString());
	}

	// -----------------------------------------------------------------------------------------------
	// cactus: cactus metric from a random edge-weighted cactus (tree ⊂ cactus ⊂ circular)
	// -----------------------------------------------------------------------------------------------

	private static Instance generateCactus(int n, Random rng, int maxWeight, double cycleProb) {
		// Grow a random cactus: start with vertex 0, then repeatedly attach, at a random EXISTING vertex a,
		// either a bridge (one new pendant vertex) or a fresh simple cycle whose other vertices are all new.
		// Because every added cycle introduces only brand-new vertices (hence brand-new edges), no two cycles
		// can ever share an edge -- exactly the cactus condition "every edge lies in at most one cycle". A
		// cactus is outer-planar, so its shortest-path metric is circular-decomposable (Kalmanson): cactus
		// metrics are a subset of circular. All n vertices are taxa.
		var edges = new ArrayList<int[]>(); // {u, v, weight}
		var present = new boolean[n][n];
		final var maxCycleLen = 6;          // cap cycle size so instances vary, rather than becoming one big ring
		var created = 1;                    // vertex 0 already exists
		var cycles = 0;
		var bridges = 0;
		while (created < n) {
			var a = rng.nextInt(created);   // cut vertex: attach the new block here
			var remaining = n - created;
			if (remaining == 1 || rng.nextDouble() >= cycleProb) {
				// bridge: a -- b with one new pendant vertex b
				var b = created++;
				addEdge(edges, present, a, b, 1 + rng.nextInt(maxWeight));
				bridges++;
			} else {
				// simple cycle a - v1 - ... - vk - a using k>=2 NEW vertices (so cycle length k+1 >= 3)
				var maxNew = Math.min(maxCycleLen - 1, remaining);
				var newCount = 2 + rng.nextInt(maxNew - 1); // in [2, maxNew]
				var prev = a;
				for (var i = 0; i < newCount; i++) {
					var b = created++;
					addEdge(edges, present, prev, b, 1 + rng.nextInt(maxWeight));
					prev = b;
				}
				addEdge(edges, present, prev, a, 1 + rng.nextInt(maxWeight)); // close the cycle
				cycles++;
			}
		}

		var total = 0L;
		for (var e : edges)
			total += e[2];

		// all-pairs shortest paths (Floyd-Warshall)
		final var INF = Integer.MAX_VALUE / 4;
		var d = new int[n][n];
		for (var row : d)
			Arrays.fill(row, INF);
		for (var i = 0; i < n; i++)
			d[i][i] = 0;
		for (var e : edges) {
			d[e[0]][e[1]] = Math.min(d[e[0]][e[1]], e[2]);
			d[e[1]][e[0]] = d[e[0]][e[1]];
		}
		for (var k = 0; k < n; k++)
			for (var i = 0; i < n; i++)
				for (var j = 0; j < n; j++)
					if (d[i][k] + d[k][j] < d[i][j])
						d[i][j] = d[i][k] + d[k][j];

		var sb = new StringBuilder();
		sb.append("type=cactus\n").append("n=").append(n).append("\n");
		sb.append("referenceLength=").append(total).append("\n");
		sb.append("referenceIsOptimal=false\n");
		sb.append("# Cactus metric: all-pairs shortest-path distances of the edge-weighted cactus below\n");
		sb.append("# (a connected graph in which every edge lies in at most one cycle; here ").append(cycles)
				.append(" cycle(s) + ").append(bridges).append(" bridge(s)).\n");
		sb.append("# A cactus is outer-planar, so this metric is circular-decomposable (Kalmanson):\n");
		sb.append("# cactus metrics are a SUBSET of circular. The cactus realizes the distances, so\n");
		sb.append("# referenceLength is an UPPER BOUND on the optimum (RazorNet should realize with length <= it).\n");
		sb.append("# vertices 0..").append(n - 1).append(" = taxa t1..t").append(n).append("; edges: 'edge u v weight'\n");
		for (var e : edges)
			sb.append("edge ").append(e[0]).append(" ").append(e[1]).append(" ").append(e[2]).append("\n");
		return new Instance(d, total, false, sb.toString());
	}

	// -----------------------------------------------------------------------------------------------
	// circle: circular-decomposable metric from a random circular split system
	// -----------------------------------------------------------------------------------------------

	private static Instance generateCircular(int n, Random rng, int maxWeight, double splitsPerTaxon, double density) {
		// Random circular ordering: order[pos] = taxon at circular position pos. Every split we use is an
		// interval of this ordering (arc of positions p+1..q), so d is a circular (Kalmanson) distance --
		// the class Neighbor-Net provably reconstructs. Default is a SPARSE system: all n trivial (pendant)
		// splits plus ~splitsPerTaxon*n random non-trivial interval splits, i.e. O(n) splits, so both the
		// realization and the split network stay tractable. The full O(n^2)-split system has an impractically
		// large realization / split network (it breaks Neighbor-Net's drawing); pass --density>0 to approach
		// it on purpose.
		var order = new int[n];
		for (var i = 0; i < n; i++)
			order[i] = i;
		shuffle(order, rng);

		var d = new int[n][n];
		var total = 0L;
		var active = 0;
		var trivial = 0;
		var splitList = new StringBuilder();
		var inArc = new boolean[n];

		if (density > 0) {
			// dense/stress mode: keep each of the n(n-1)/2 interval splits with probability 'density'
			for (var p = 0; p < n; p++)
				for (var q = p + 1; q < n; q++) {
					if (rng.nextDouble() >= density)
						continue;
					fillArc(inArc, order, p, q);
					var w = 1 + rng.nextInt(maxWeight);
					applySplit(inArc, w, d, splitList);
					total += w;
					active++;
				}
		} else {
			// sparse mode: all n trivial (pendant) splits ...
			for (var t = 0; t < n; t++) {
				Arrays.fill(inArc, false);
				inArc[t] = true;
				var w = 1 + rng.nextInt(maxWeight);
				applySplit(inArc, w, d, splitList);
				total += w;
				active++;
				trivial++;
			}
			// ... plus ~splitsPerTaxon*n random non-trivial interval splits (arc size 2..n-2)
			var candidates = new ArrayList<int[]>();
			for (var p = 0; p < n; p++)
				for (var q = p + 1; q < n; q++)
					if (q - p >= 2 && q - p <= n - 2)
						candidates.add(new int[]{p, q});
			java.util.Collections.shuffle(candidates, rng);
			var take = Math.min((int) Math.round(splitsPerTaxon * n), candidates.size());
			for (var c = 0; c < take; c++) {
				var pq = candidates.get(c);
				fillArc(inArc, order, pq[0], pq[1]);
				var w = 1 + rng.nextInt(maxWeight);
				applySplit(inArc, w, d, splitList);
				total += w;
				active++;
			}
		}

		var ordering = new StringBuilder();
		for (var pos = 0; pos < n; pos++)
			ordering.append(pos == 0 ? "" : " ").append("t").append(order[pos] + 1);

		var sb = new StringBuilder();
		sb.append("type=circle\n").append("n=").append(n).append("\n");
		sb.append("referenceLength=").append(total).append("\n");
		sb.append("referenceIsOptimal=false\n");
		sb.append("# Sparse circular (Kalmanson) metric: d(i,j) = sum of weights of the circular splits separating i,j,\n");
		sb.append("# every split an interval of the circular ordering below -- the class Neighbor-Net provably\n");
		sb.append("# reconstructs, kept sparse (O(n) splits) so the realization and split network stay tractable.\n");
		sb.append("circularOrdering: ").append(ordering).append("\n");
		if (density > 0)
			sb.append("splits=").append(active).append(" (dense mode, density=").append(density).append(")\n");
		else
			sb.append("splits=").append(active).append(" (").append(trivial).append(" trivial + ").append(active - trivial)
					.append(" non-trivial; splitsPerTaxon=").append(splitsPerTaxon).append(")\n");
		sb.append("# NOTE: referenceLength (").append(total).append(") is the TOTAL SPLIT WEIGHT -- a structural measure,\n");
		sb.append("# NOT a realizing-network edge length (a realizing network is typically longer). The meaningful\n");
		sb.append("# check for circle instances is that RazorNet realizes the distances exactly.\n");
		sb.append("# splits: 'split weight=W sideA={...}' lists taxa on one side of the interval split; the rest are the other side\n");
		sb.append(splitList);
		return new Instance(d, total, false, sb.toString());
	}

	private static void fillArc(boolean[] inArc, int[] order, int p, int q) {
		Arrays.fill(inArc, false);
		for (var pos = p + 1; pos <= q; pos++)
			inArc[order[pos]] = true;
	}

	private static void applySplit(boolean[] inArc, int w, int[][] d, StringBuilder splitList) {
		splitList.append("split weight=").append(w).append(" sideA={");
		var first = true;
		for (var t = 0; t < inArc.length; t++)
			if (inArc[t]) {
				splitList.append(first ? "" : ",").append("t").append(t + 1);
				first = false;
			}
		splitList.append("}\n");
		for (var i = 0; i < inArc.length; i++)
			if (inArc[i])
				for (var j = 0; j < inArc.length; j++)
					if (!inArc[j]) {
						d[i][j] += w;
						d[j][i] += w;
					}
	}

	// -----------------------------------------------------------------------------------------------
	// general: arbitrary metric = shortest-path metric of a random connected weighted graph
	// -----------------------------------------------------------------------------------------------

	private static Instance generateGeneral(int n, Random rng, int maxWeight) {
		var edges = new ArrayList<int[]>(); // {u, v, weight}
		var present = new boolean[n][n];

		// random spanning tree (guarantees connectivity): attach each vertex to a random earlier one
		var perm = new int[n];
		for (var i = 0; i < n; i++)
			perm[i] = i;
		for (var i = n - 1; i > 0; i--) {
			var j = rng.nextInt(i + 1);
			var tmp = perm[i];
			perm[i] = perm[j];
			perm[j] = tmp;
		}
		for (var i = 1; i < n; i++) {
			var u = perm[i];
			var v = perm[rng.nextInt(i)];
			addEdge(edges, present, u, v, 1 + rng.nextInt(maxWeight));
		}
		// extra edges to make it a genuine (cyclic) network, ~n attempts
		for (var k = 0; k < n; k++) {
			var u = rng.nextInt(n);
			var v = rng.nextInt(n);
			if (u != v && !present[u][v])
				addEdge(edges, present, u, v, 1 + rng.nextInt(maxWeight));
		}

		var total = 0L;
		for (var e : edges)
			total += e[2];

		// all-pairs shortest paths (Floyd-Warshall)
		final var INF = Integer.MAX_VALUE / 4;
		var d = new int[n][n];
		for (var row : d)
			Arrays.fill(row, INF);
		for (var i = 0; i < n; i++)
			d[i][i] = 0;
		for (var e : edges) {
			d[e[0]][e[1]] = Math.min(d[e[0]][e[1]], e[2]);
			d[e[1]][e[0]] = d[e[0]][e[1]];
		}
		for (var k = 0; k < n; k++)
			for (var i = 0; i < n; i++)
				for (var j = 0; j < n; j++)
					if (d[i][k] + d[k][j] < d[i][j])
						d[i][j] = d[i][k] + d[k][j];

		var sb = new StringBuilder();
		sb.append("type=general\n").append("n=").append(n).append("\n");
		sb.append("referenceLength=").append(total).append("\n");
		sb.append("referenceIsOptimal=false\n");
		sb.append("# General metric: all-pairs shortest-path distances of the random weighted graph below.\n");
		sb.append("# That graph realizes the distances, so referenceLength is an UPPER BOUND on the optimum\n");
		sb.append("# (RazorNet should realize the distances with total length <= referenceLength).\n");
		sb.append("# vertices 0..").append(n - 1).append(" = taxa t1..t").append(n).append("; edges: 'edge u v weight'\n");
		for (var e : edges)
			sb.append("edge ").append(e[0]).append(" ").append(e[1]).append(" ").append(e[2]).append("\n");
		return new Instance(d, total, false, sb.toString());
	}

	private static void addEdge(List<int[]> edges, boolean[][] present, int u, int v, int w) {
		edges.add(new int[]{u, v, w});
		present[u][v] = present[v][u] = true;
	}

	private static void shuffle(int[] a, Random rng) {
		for (var i = a.length - 1; i > 0; i--) {
			var j = rng.nextInt(i + 1);
			var t = a[i];
			a[i] = a[j];
			a[j] = t;
		}
	}

	// -----------------------------------------------------------------------------------------------
	// reticulate: sum of two tree metrics (tree2 = tree1 after a few SPR moves) -- conflicting-trees signal
	// -----------------------------------------------------------------------------------------------

	private static Instance generateReticulate(int n, Random rng, int maxWeight, int numSpr, double reticulationScale) {
		// Tree 1: random topology + random integer edge weights.
		var edges1 = buildTreeTopology(n, rng);
		var numNodes1 = maxNodeId(edges1) + 1;
		var w1 = new int[edges1.size()];
		var len1 = 0L;
		for (var i = 0; i < w1.length; i++) {
			w1[i] = 1 + rng.nextInt(maxWeight);
			len1 += w1[i];
		}
		var d1 = treeLeafDistances(n, edges1, w1, numNodes1);

		// Tree 2: tree 1's topology after numSpr SPR moves, then fresh (scaled) integer weights.
		var adj = new TreeMap<Integer, TreeSet<Integer>>();
		for (var e : edges1)
			addUndirected(adj, e[0], e[1]);
		var nextNode = new int[]{numNodes1};
		var sprDone = 0;
		for (var s = 0; s < numSpr; s++)
			if (sprMove(adj, rng, nextNode)) sprDone++;
		var edges2 = reindexInternalNodes(edgesFromAdj(adj), n);
		var numNodes2 = maxNodeId(edges2) + 1;
		var w2 = new int[edges2.size()];
		var len2 = 0L;
		for (var i = 0; i < w2.length; i++) {
			w2[i] = Math.max(1, (int) Math.round((1 + rng.nextInt(maxWeight)) * reticulationScale));
			len2 += w2[i];
		}
		var d2 = treeLeafDistances(n, edges2, w2, numNodes2);

		// D = D(tree1) + D(tree2): the sum of two tree metrics is again a metric (symmetric, triangle-obeying).
		var d = new int[n][n];
		for (var i = 0; i < n; i++)
			for (var j = 0; j < n; j++)
				d[i][j] = d1[i][j] + d2[i][j];

		var ref = taxaOnlyLength(d); // taxa-only (auxiliary-free) essential-edge length: a valid upper bound

		var sb = new StringBuilder();
		sb.append("type=reticulate\n").append("n=").append(n).append("\n");
		sb.append("referenceLength=").append(ref).append("\n");
		sb.append("referenceIsOptimal=false\n");
		sb.append("# Reticulate metric: D = D(tree1) + D(tree2), the SUM of two tree metrics on the same taxa --\n");
		sb.append("# tree2 is tree1 after ").append(sprDone).append(" SPR move(s) (requested ").append(numSpr)
				.append("), tree2 edge weights scaled by ").append(reticulationScale).append(".\n");
		sb.append("# The classic 'two conflicting trees' signal that Neighbor-Net is used to display; there is no\n");
		sb.append("# known minimum realization. referenceLength is the taxa-only (auxiliary-free) essential-edge length,\n");
		sb.append("# a valid UPPER BOUND (RazorNet keeps this candidate, so it realizes D with length <= it).\n");
		sb.append("# len(tree1)=").append(len1).append("  len(tree2)=").append(len2).append("\n");
		sb.append("# nodes: leaves 0..").append(n - 1).append(" = taxa t1..t").append(n).append("; internal (auxiliary) nodes have higher ids.\n");
		sb.append("# tree1 edges: 'tree1 u v weight'\n");
		for (var i = 0; i < edges1.size(); i++)
			sb.append("tree1 ").append(edges1.get(i)[0]).append(" ").append(edges1.get(i)[1]).append(" ").append(w1[i]).append("\n");
		sb.append("# tree2 edges: 'tree2 u v weight'\n");
		for (var i = 0; i < edges2.size(); i++)
			sb.append("tree2 ").append(edges2.get(i)[0]).append(" ").append(edges2.get(i)[1]).append(" ").append(w2[i]).append("\n");
		return new Instance(d, ref, false, sb.toString());
	}

	/**
	 * Random binary tree topology by stepwise leaf addition; leaves 0..n-1 are taxa, internal nodes n.. .
	 */
	private static List<int[]> buildTreeTopology(int n, Random rng) {
		var edges = new ArrayList<int[]>();
		var nextInternal = n;
		if (n >= 2)
			edges.add(new int[]{0, 1});
		for (var k = 2; k < n; k++) {
			var e = edges.remove(rng.nextInt(edges.size()));
			var m = nextInternal++;
			edges.add(new int[]{e[0], m});
			edges.add(new int[]{m, e[1]});
			edges.add(new int[]{m, k});
		}
		return edges;
	}

	/**
	 * Leaf-to-leaf path distances of a weighted tree (one DFS per leaf); node ids must be 0..numNodes-1.
	 */
	private static int[][] treeLeafDistances(int n, List<int[]> edges, int[] weight, int numNodes) {
		var adj = new ArrayList<List<int[]>>();
		for (var i = 0; i < numNodes; i++)
			adj.add(new ArrayList<>());
		for (var i = 0; i < edges.size(); i++) {
			var e = edges.get(i);
			adj.get(e[0]).add(new int[]{e[1], weight[i]});
			adj.get(e[1]).add(new int[]{e[0], weight[i]});
		}
		var d = new int[n][n];
		for (var s = 0; s < n; s++) {
			var dist = new int[numNodes];
			Arrays.fill(dist, -1);
			dist[s] = 0;
			var stack = new ArrayList<Integer>();
			stack.add(s);
			while (!stack.isEmpty()) {
				int u = stack.remove(stack.size() - 1);
				for (var nb : adj.get(u))
					if (dist[nb[0]] < 0) {
						dist[nb[0]] = dist[u] + nb[1];
						stack.add(nb[0]);
					}
			}
			for (var t = 0; t < n; t++)
				d[s][t] = dist[t];
		}
		return d;
	}

	/**
	 * One SPR (subtree prune-and-regraft) move on an unrooted binary tree held as an adjacency map. Picks a
	 * prune edge whose attachment endpoint u is internal (degree 3), detaches u's neighbor v (and its subtree),
	 * suppresses u (joining its two other neighbors), then regrafts the subtree by subdividing a random edge of
	 * the remaining tree with a new internal node. Keeps the tree binary and the leaf set unchanged.
	 *
	 * @return true if a move was made
	 */
	private static boolean sprMove(TreeMap<Integer, TreeSet<Integer>> adj, Random rng, int[] nextNode) {
		var edgeList = new ArrayList<int[]>();
		for (var a : adj.keySet())
			for (var b : adj.get(a))
				if (a < b) edgeList.add(new int[]{a, b});
		Collections.shuffle(edgeList, rng);
		for (var e : edgeList) {
			int a = e[0], b = e[1];
			var cands = new ArrayList<Integer>();
			if (adj.get(a).size() == 3) cands.add(a);
			if (adj.get(b).size() == 3) cands.add(b);
			if (cands.isEmpty()) continue;                 // need an internal attachment node to suppress
			int u = cands.get(rng.nextInt(cands.size()));
			int v = (u == a) ? b : a;
			var others = new ArrayList<>(adj.get(u));       // u has degree 3
			others.remove((Integer) v);
			int w1 = others.get(0), w2 = others.get(1);
			// prune v's subtree, suppress u
			removeUndirected(adj, u, v);
			removeUndirected(adj, u, w1);
			removeUndirected(adj, u, w2);
			adj.remove(u);
			addUndirected(adj, w1, w2);
			// regraft target: an edge of the component containing w1, excluding the just-created (w1,w2)
			var main = component(adj, w1);
			var targets = new ArrayList<int[]>();
			for (var p : main)
				for (var q : adj.get(p))
					if (p < q && main.contains(q) && !(Math.min(p, q) == Math.min(w1, w2) && Math.max(p, q) == Math.max(w1, w2)))
						targets.add(new int[]{p, q});
			int p, q;
			if (targets.isEmpty()) {           // only the undo edge exists (tiny tree): reattach there
				p = Math.min(w1, w2);
				q = Math.max(w1, w2);
			} else {
				var t = targets.get(rng.nextInt(targets.size()));
				p = t[0];
				q = t[1];
			}
			var m = nextNode[0]++;
			removeUndirected(adj, p, q);
			addUndirected(adj, p, m);
			addUndirected(adj, m, q);
			addUndirected(adj, m, v);
			return true;
		}
		return false;
	}

	private static void addUndirected(TreeMap<Integer, TreeSet<Integer>> adj, int a, int b) {
		adj.computeIfAbsent(a, k -> new TreeSet<>()).add(b);
		adj.computeIfAbsent(b, k -> new TreeSet<>()).add(a);
	}

	private static void removeUndirected(TreeMap<Integer, TreeSet<Integer>> adj, int a, int b) {
		if (adj.containsKey(a)) adj.get(a).remove(b);
		if (adj.containsKey(b)) adj.get(b).remove(a);
	}

	private static Set<Integer> component(TreeMap<Integer, TreeSet<Integer>> adj, int start) {
		var seen = new TreeSet<Integer>();
		var stack = new ArrayList<Integer>();
		stack.add(start);
		seen.add(start);
		while (!stack.isEmpty()) {
			var x = stack.remove(stack.size() - 1);
			for (var y : adj.get(x))
				if (seen.add(y)) stack.add(y);
		}
		return seen;
	}

	private static List<int[]> edgesFromAdj(TreeMap<Integer, TreeSet<Integer>> adj) {
		var out = new ArrayList<int[]>();
		for (var a : adj.keySet())
			for (var b : adj.get(a))
				if (a < b) out.add(new int[]{a, b});
		return out;
	}

	/**
	 * Relabel internal nodes to n, n+1, ... (leaves 0..n-1 stay fixed) so ids are contiguous 0..numNodes-1.
	 */
	private static List<int[]> reindexInternalNodes(List<int[]> edges, int n) {
		var map = new TreeMap<Integer, Integer>();
		var internals = new TreeSet<Integer>();
		for (var e : edges) {
			if (e[0] >= n) internals.add(e[0]);
			if (e[1] >= n) internals.add(e[1]);
		}
		var next = n;
		for (var id : internals)
			map.put(id, next++);
		var out = new ArrayList<int[]>();
		for (var e : edges) {
			var a = e[0] < n ? e[0] : map.get(e[0]);
			var b = e[1] < n ? e[1] : map.get(e[1]);
			out.add(new int[]{a, b});
		}
		return out;
	}

	private static int maxNodeId(List<int[]> edges) {
		var max = 0;
		for (var e : edges)
			max = Math.max(max, Math.max(e[0], e[1]));
		return max;
	}

	/**
	 * Length of the taxa-only (auxiliary-free) realization of a metric d: keep edge (i,j) unless some third point
	 * k gives an equal-or-shorter two-hop path (which for a metric implies an equal-or-shorter path of any
	 * length), and sum the kept edge weights. A valid UPPER BOUND on the optimal realization length.
	 */
	private static long taxaOnlyLength(int[][] d) {
		var n = d.length;
		var total = 0L;
		for (var i = 0; i < n; i++)
			for (var j = i + 1; j < n; j++) {
				var redundant = false;
				for (var k = 0; k < n && !redundant; k++) {
					if (k == i || k == j) continue;
					if (d[i][k] > 0 && d[k][j] > 0 && d[i][k] + d[k][j] <= d[i][j])
						redundant = true;
				}
				if (!redundant)
					total += d[i][j];
			}
		return total;
	}

	// -----------------------------------------------------------------------------------------------
	// output
	// -----------------------------------------------------------------------------------------------

	private static void writeDist(Path path, String type, int n, long seed, Instance instance) throws IOException {
		var d = instance.distances();
		var sb = new StringBuilder();
		sb.append("# RazorNet test instance\n");
		sb.append("# type=").append(type).append(" n=").append(n).append(" seed=").append(seed)
				.append(" referenceLength=").append(instance.referenceLength())
				.append(" referenceIsOptimal=").append(instance.referenceIsOptimal()).append("\n");
		sb.append(n).append("\n");
		for (var i = 0; i < n; i++) {
			sb.append("t").append(i + 1);
			for (var j = 0; j < n; j++)
				sb.append("\t").append(d[i][j]);
			sb.append("\n");
		}
		Files.writeString(path, sb.toString());
	}

	private static void writeTruth(Path path, long seed, Instance instance) throws IOException {
		Files.writeString(path, "# Ground truth for the sibling .dist instance\nseed=" + seed + "\n" + instance.truth());
	}

	private static int[] parseInts(String s) {
		var parts = s.trim().split("[,\\s]+");
		var result = new int[parts.length];
		for (var i = 0; i < parts.length; i++)
			result[i] = Integer.parseInt(parts[i]);
		return result;
	}

	private static void printUsage() {
		System.err.println("""
				GenerateTestInstances - create RazorNet distance-matrix test instances

				Tiers (hierarchy tree ⊂ cactus ⊂ circular ⊂ general):
				  tree    additive metric of a random weighted tree (reference length = true optimum)
				  cactus  shortest-path metric of a random edge-weighted cactus (every edge in ≤1 cycle)
				  circle  circular-decomposable (Kalmanson) metric of a random circular split system
				  general shortest-path metric of an arbitrary random weighted graph
				  reticulate  D = D(tree1)+D(tree2), tree2 = tree1 after a few SPR moves (conflicting-trees signal)

				Options:
				  -o, --output DIR      output directory (default: tests); subdirs {tree,cactus,circle,general,reticulate}-distances are created
				  -n, --sizes LIST      comma/space separated taxon counts (default: 5,10,20,40,80)
				  -r, --reps N          instances per size per type (default: 1)
				  -s, --seed N          master random seed (default: 42)
				  -w, --maxWeight N       maximum edge/split weight (default: 100)
				  -k, --splitsPerTaxon F  circle tier: non-trivial circular splits per taxon (default: 1 = sparse, tractable)
				  -d, --density F         circle tier: dense STRESS mode = keep this fraction of ALL n(n-1)/2 splits instead (default: 0 = off)
				  -c, --cactusCycleProb F cactus tier: probability of a cycle (vs a bridge) at each growth step (default: 0.5)
				  -p, --sprMoves N        reticulate tier: number of SPR moves turning tree1 into tree2 (default: 2)
				  -x, --reticulationScale F reticulate tier: scale on tree2 edge weights, 1.0 = equal signals (default: 1.0)
				  -h, --help              show this help

				Each instance is written as <type>-nNNNN[.rR].dist (a plain distance matrix that both RazorNet
				and SplitsTree read) plus a .truth.txt describing the generating tree/cactus/split-system/graph
				(or tree pair) and its total length.""");
	}
}
