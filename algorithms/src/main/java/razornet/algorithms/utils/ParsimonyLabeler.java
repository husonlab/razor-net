/*
 * ParsimonyLabeler.java Copyright (C) 2026 Daniel H. Huson
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
package razornet.algorithms.utils;

import java.util.*;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

/**
 * Weighted small parsimony on a connected, undirected graph in which some nodes are labeled by sequences: assigns
 * a sequence to every node such that the total weight of the edges whose endpoints differ, summed over the sites,
 * is as small as possible. With weight 1 on every edge this minimizes the number of changes; with the geodesic
 * betweenness as weight it minimizes the excess of the displayed distances.
 * <p>
 * Each site is solved on its own. A site is either DNA ({A,C,G,T}, with IUPAC ambiguity codes) or binary ({0,1});
 * a gap is an ordinary state, and '?' and 'N' allow any state unless {@link #setMissingIsState} is set. Sites with
 * at most two states are solved exactly by a minimum s-t cut; sites with three or more states by alpha-expansion
 * (Boykov, Veksler &amp; Zabih 2001) from several starting labelings. Maximum flow is computed with Dinic's
 * algorithm. Nodes and edges are processed in the order of their string representations, so that the result is
 * deterministic.
 */
public class ParsimonyLabeler<Node, Edge> {
	/**
	 * Maximum α-expansion iterations per restart.
	 */
	private static final int MAX_ALPHA_ITERATIONS = 50;

	/**
	 * Number of random restarts (in addition to the deterministic greedy initialization).
	 */
	private static final int NUM_RANDOM_RESTARTS = 4;

	/**
	 * Seed for random restart initializations (deterministic by default).
	 */
	private static final long DEFAULT_SEED = 42L;

	/**
	 * Set to true to log convergence-cap warnings.
	 */
	public static boolean verbose = false;

	private final List<Node> nodes;          // sorted, deterministic order
	private final List<Edge> edges;          // sorted, deterministic order
	private final Function<Edge, Node> u;
	private final Function<Edge, Node> v;
	private final long seed;
	private final ToDoubleFunction<Edge> edgeWeight;

	private boolean missingIsState = false;   // see setMissingIsState

	public ParsimonyLabeler(Iterable<Node> nodes, Iterable<Edge> edges,
							Function<Edge, Node> u, Function<Edge, Node> v) {
		this(nodes, edges, u, v, DEFAULT_SEED);
	}

	/**
	 * @param seed seed for random restarts, ensuring reproducibility
	 */
	public ParsimonyLabeler(Iterable<Node> nodes, Iterable<Edge> edges,
							Function<Edge, Node> u, Function<Edge, Node> v, long seed) {
		this(nodes, edges, u, v, seed, e -> 1.0);
	}

	/**
	 * @param edgeWeight per-edge cost of a state change on that edge (see the six-arg constructor)
	 */
	public ParsimonyLabeler(Iterable<Node> nodes, Iterable<Edge> edges,
							Function<Edge, Node> u, Function<Edge, Node> v, ToDoubleFunction<Edge> edgeWeight) {
		this(nodes, edges, u, v, DEFAULT_SEED, edgeWeight);
	}

	/**
	 * @param seed       seed for the random restarts
	 * @param edgeWeight the cost of a state change on an edge: 1 for parsimony, the geodesic betweenness of the edge
	 *                   to minimize the excess
	 */
	public ParsimonyLabeler(Iterable<Node> nodes, Iterable<Edge> edges,
							Function<Edge, Node> u, Function<Edge, Node> v, long seed, ToDoubleFunction<Edge> edgeWeight) {
		this.u = u;
		this.v = v;
		this.seed = seed;
		this.edgeWeight = edgeWeight;
		// nodes and edges are sorted by their string representation, for a deterministic order
		var nodeList = new ArrayList<Node>();
		nodes.forEach(nodeList::add);
		nodeList.sort(Comparator.comparing(Object::toString));
		this.nodes = Collections.unmodifiableList(nodeList);

		var edgeList = new ArrayList<Edge>();
		edges.forEach(edgeList::add);
		// Sort edges by (min(u,v), max(u,v)) using their toString representations
		edgeList.sort(Comparator
				.comparing((Edge e) -> {
					var s1 = u.apply(e).toString();
					var s2 = v.apply(e).toString();
					return s1.compareTo(s2) <= 0 ? s1 : s2;
				})
				.thenComparing(e -> {
					var s1 = u.apply(e).toString();
					var s2 = v.apply(e).toString();
					return s1.compareTo(s2) <= 0 ? s2 : s1;
				}));
		this.edges = Collections.unmodifiableList(edgeList);
	}

	/**
	 * Whether '?', 'N' and the IUPAC ambiguity codes are ordinary states rather than wildcards (default false). This
	 * has to match the way in which the distances were computed.
	 */
	public void setMissingIsState(boolean missingIsState) {
		this.missingIsState = missingIsState;
	}

	/**
	 * Main entry: label all nodes for all sites.
	 *
	 * @param knownSeqs map of pre-labeled nodes to their sequences (all equal length).
	 *                  DNA columns may contain IUPAC codes; binary columns contain 0/1.
	 *                  '?' is unconstrained in either column type unless {@link #setMissingIsState} is set;
	 *                  '-' is an ordinary state either way.
	 * @return map of every node to its deduced sequence (one concrete character per site)
	 */
	public Map<Node, String> labelAllSites(Map<Node, String> knownSeqs) {
		if (knownSeqs.isEmpty())
			throw new IllegalArgumentException("At least one labeled node required");

		int L = knownSeqs.values().iterator().next().length();
		for (String s : knownSeqs.values())
			if (s.length() != L)
				throw new IllegalArgumentException("All sequences must have equal length");

		// Per-site universe detection and allowed-mask construction.
		// Each site is independent — DNA and binary columns can be mixed in the alignment.
		List<Universe> siteUniverses = new ArrayList<>(L);
		List<Map<Node, Integer>> siteAllowed = new ArrayList<>(L);
		for (int j = 0; j < L; j++) {
			Universe universe = detectSiteUniverse(knownSeqs, j);
			siteUniverses.add(universe);
			Map<Node, Integer> m = new HashMap<>();
			for (Map.Entry<Node, String> e : knownSeqs.entrySet()) {
				m.put(e.getKey(), expandChar(e.getValue().charAt(j), universe));
			}
			siteAllowed.add(m);
		}

		// Solve each site independently
		List<Map<Node, Character>> siteLabels = new ArrayList<>(L);
		for (int j = 0; j < L; j++) {
			Universe universe = siteUniverses.get(j);
			Map<Node, Integer> allowedMasks = siteAllowed.get(j);

			// Determine feasible alphabet (restricting to states that appear in any constraint).
			// Iterating universe.bases preserves a deterministic order in sigma.
			int unionMask = 0;
			for (int m : allowedMasks.values()) unionMask |= m;
			if (unionMask == 0)
				throw new IllegalArgumentException("No allowed states at site " + j);

			Set<Character> sigma = new LinkedHashSet<>();
			for (int i = 0; i < universe.bases.length; i++)
				if ((unionMask & (1 << i)) != 0) sigma.add(universe.bases[i]);

			// a node whose allowed set is a single state is pinned; only the free nodes have to be solved for
			Map<Node, Character> pinned = new HashMap<>();
			List<Node> free = new ArrayList<>();
			for (Node x : nodes) {
				int mask = allowedMasks.getOrDefault(x, universe.allMask());
				if (Integer.bitCount(mask) == 1)
					pinned.put(x, firstAllowed(mask, universe));
				else
					free.add(x);   // includes the empty mask, which the solvers reject
			}

			Map<Node, Character> labels;
			if (sigma.size() <= 1) {
				char a = sigma.iterator().next();
				labels = new HashMap<>();
				for (Node x : nodes) labels.put(x, a);
			} else if (free.isEmpty()) {
				// nothing to optimise: the labeling is determined
				labels = new HashMap<>(pinned);
			} else {
				var solveNodes = solveNodesFor(free);
				var solveEdges = solveEdgesFor(free);
				if (sigma.size() == 2)
					labels = solveBinarySiteWithAmbiguity(allowedMasks, sigma, universe, solveNodes, solveEdges);
				else
					labels = solveSiteAlphaExpansionMultiStart(allowedMasks, sigma, universe, j, solveNodes, solveEdges);
				// nodes left out of the reduced problem keep the one state they were allowed
				for (Map.Entry<Node, Character> e : pinned.entrySet())
					labels.putIfAbsent(e.getKey(), e.getValue());
			}

			// Final feasibility check: ensure labels respect allowed masks (defensive)
			for (Node x : nodes) {
				int mask = allowedMasks.getOrDefault(x, universe.allMask());
				char bx = labels.get(x);
				if (!allowed(mask, bx, universe)) labels.put(x, firstAllowed(mask, universe));
			}

			siteLabels.add(labels);
		}

		// Concatenate per-site
		Map<Node, StringBuilder> agg = new HashMap<>();
		for (Node x : nodes) agg.put(x, new StringBuilder(L));
		for (int j = 0; j < L; j++) {
			Map<Node, Character> lab = siteLabels.get(j);
			for (Node x : nodes) agg.get(x).append(lab.get(x));
		}
		Map<Node, String> out = new HashMap<>();
		for (Map.Entry<Node, StringBuilder> e : agg.entrySet())
			out.put(e.getKey(), e.getValue().toString());
		return out;
	}

	/**
	 * the nodes that the solvers have to consider: the free nodes, and the pinned nodes adjacent to a free node. An
	 * edge between two pinned nodes contributes a constant and is left out. The order follows {@code nodes}.
	 */
	private List<Node> solveNodesFor(List<Node> free) {
		var freeSet = new HashSet<>(free);
		var involved = new HashSet<>(free);
		for (Edge e : edges) {
			Node a = u.apply(e), b = v.apply(e);
			if (freeSet.contains(a) || freeSet.contains(b)) {
				involved.add(a);
				involved.add(b);
			}
		}
		var result = new ArrayList<Node>();
		for (Node x : nodes)
			if (involved.contains(x)) result.add(x);
		return result;
	}

	/** The edges with at least one free endpoint; the rest contribute a constant (see solveNodesFor). */
	private List<Edge> solveEdgesFor(List<Node> free) {
		var freeSet = new HashSet<>(free);
		var result = new ArrayList<Edge>();
		for (Edge e : edges)
			if (freeSet.contains(u.apply(e)) || freeSet.contains(v.apply(e))) result.add(e);
		return result;
	}

	// ------------------------ Alphabet abstraction ------------------------

	/**
	 * A character "universe" defines the alphabet for a single site. The {@code bases} array
	 * is the canonical ordering; an int "mask" encodes an allowed-set as bits where bit i is
	 * set iff {@code bases[i]} is allowed. Determinism throughout the solver depends on this
	 * fixed ordering.
	 */
	private record Universe(char[] bases) {
		int allMask() {
			int m = 0;
			for (int i = 0; i < bases.length; i++) m |= (1 << i);
			return m;
		}

		/**
		 * Mask with just {@code c} set, or 0 if {@code c} isn't in this universe.
		 */
		int singleton(char c) {
			for (int i = 0; i < bases.length; i++)
				if (bases[i] == c) return 1 << i;
			return 0;
		}
	}

	// a gap is an ordinary state; '?', 'N' and the IUPAC codes are wildcards
	private static final Universe DNA = new Universe(new char[]{'A', 'C', 'G', 'T', '-'});
	private static final Universe BINARY = new Universe(new char[]{'0', '1', '-'});

	// used when missingIsState is set: every character is a state of its own ('A','C','G','T','-' come first, as
	// in DNA)
	private static final Universe DNA_ALL_STATES = new Universe(new char[]
			{'A', 'C', 'G', 'T', '-', 'U', 'N', '?', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'});
	private static final Universe BINARY_ALL_STATES = new Universe(new char[]{'0', '1', '-', 'N', '?'});

	/**
	 * Classify one column. A column with any DNA character (ACGT[U] or IUPAC code) is DNA;
	 * a column with 0/1 characters and no DNA characters is binary; an all-missing column
	 * defaults to DNA (its labeling cost is zero regardless of the choice).
	 *
	 * @throws IllegalArgumentException if the column mixes DNA and 0/1 characters
	 */
	private Universe detectSiteUniverse(Map<Node, String> seqs, int site) {
		boolean hasDna = false, has01 = false;
		for (String s : seqs.values()) {
			char c = Character.toUpperCase(s.charAt(site));
			switch (c) {
				case 'A', 'C', 'G', 'T', 'U',
					 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V' -> hasDna = true;
				case '0', '1' -> has01 = true;
				case '?', '-', 'N' -> { /* ambiguous, contributes no info to detection */ }
				default -> throw new IllegalArgumentException(
						"Unsupported character '" + c + "' at site " + site);
			}
		}
		if (hasDna && has01) {
			throw new IllegalArgumentException(
					"Site " + site + " mixes DNA and binary (0/1) characters; not supported.");
		}
		if (has01)
			return missingIsState ? BINARY_ALL_STATES : BINARY;
		else
			return missingIsState ? DNA_ALL_STATES : DNA;
	}

	/**
	 * Expand a single character to its allowed-set mask within the given universe.
	 * '?' and 'N' expand to all bases; '-' is a state and expands to itself. IUPAC ambiguity codes (only valid in DNA)
	 * expand to their canonical sets. '0' and '1' are accepted only in a binary universe.
	 * <p>
	 * Under {@link #setMissingIsState} none of that applies: each character maps to itself alone.
	 */
	private int expandChar(char c, Universe universe) {
		char up = Character.toUpperCase(c);
		if (missingIsState) {
			var mask = universe.singleton(up);
			if (mask == 0)
				throw new IllegalArgumentException("Unsupported character '" + c + "'");
			return mask;
		}
		if (up == '-') return universe.singleton('-');   // a gap is a state, not a wildcard
		if (up == '?' || up == 'N') return universe.allMask();

		if (universe.singleton('0') != 0) {   // a binary column
			return switch (up) {
				case '0' -> universe.singleton('0');
				case '1' -> universe.singleton('1');
				default -> throw new IllegalArgumentException(
						"Unsupported character '" + c + "' in binary site");
			};
		}

		// DNA universe — bits are A=1, C=2, G=4, T=8
		return switch (up) {
			case 'A' -> 0b0001;
			case 'C' -> 0b0010;
			case 'G' -> 0b0100;
			case 'T', 'U' -> 0b1000;

			case 'R' -> 0b0001 | 0b0100;          // A or G
			case 'Y' -> 0b0010 | 0b1000;          // C or T
			case 'S' -> 0b0010 | 0b0100;          // C or G
			case 'W' -> 0b0001 | 0b1000;          // A or T
			case 'K' -> 0b0100 | 0b1000;          // G or T
			case 'M' -> 0b0001 | 0b0010;          // A or C

			case 'B' -> 0b0010 | 0b0100 | 0b1000; // not A
			case 'D' -> 0b0001 | 0b0100 | 0b1000; // not C
			case 'H' -> 0b0001 | 0b0010 | 0b1000; // not G
			case 'V' -> 0b0001 | 0b0010 | 0b0100; // not T

			default -> throw new IllegalArgumentException("Unsupported IUPAC base: '" + c + "'");
		};
	}

	private static boolean allowed(int mask, char base, Universe u) {
		int bm = u.singleton(base);
		return bm != 0 && (mask & bm) != 0;
	}

	private static char firstAllowed(int mask, Universe u) {
		for (int i = 0; i < u.bases.length; i++)
			if ((mask & (1 << i)) != 0) return u.bases[i];
		throw new IllegalArgumentException("Empty allowed mask");
	}

	private static char randomAllowed(int mask, Universe u, Random rng) {
		var list = new ArrayList<Character>(u.bases.length);
		for (int i = 0; i < u.bases.length; i++)
			if ((mask & (1 << i)) != 0) list.add(u.bases[i]);
		if (list.isEmpty())
			throw new IllegalArgumentException("Empty allowed mask");
		return list.get(rng.nextInt(list.size()));
	}

	// ------------------------ Binary (|Σ|=2) exact via min-cut ------------------------

	/**
	 * Site with feasible alphabet of size 2 (either an actual binary site, or a DNA site
	 * whose constraints restrict everything to two bases). Domain restrictions are encoded
	 * as infinite unary penalties. Convention: nodes on S-side are labeled sA; on T-side sB.
	 */
	private Map<Node, Character> solveBinarySiteWithAmbiguity(
			Map<Node, Integer> allowedMasks, Set<Character> sigma, Universe universe,
			List<Node> solveNodes, List<Edge> solveEdges) {
		Iterator<Character> it = sigma.iterator();
		char sA = it.next();
		char sB = it.next();

		FlowNetwork<Node> FN = new FlowNetwork<>();
		for (Node x : solveNodes) FN.addNode(x);
		FN.addSourceSink();

		final double INF = 1e12;

		for (Node x : solveNodes) {
			int mask = allowedMasks.getOrDefault(x, universe.allMask());
			boolean allowA = allowed(mask, sA, universe);
			boolean allowB = allowed(mask, sB, universe);
			if (!allowA && !allowB) {
				throw new IllegalArgumentException(
						"Node has no allowed state among {" + sA + "," + sB + "}");
			}
			if (!allowA) FN.addCapacity(x, FN.T, INF); // forbid S-side (sA)
			if (!allowB) FN.addCapacity(FN.S, x, INF); // forbid T-side (sB)
		}

		for (Edge e : solveEdges) {
			FN.addUndirectedCapacity(u.apply(e), v.apply(e), edgeWeight.applyAsDouble(e));
		}

		FN.maxFlow();
		Set<Node> Sside = FN.minCutSourceSide();

		Map<Node, Character> lab = new HashMap<>();
		for (Node x : solveNodes) lab.put(x, Sside.contains(x) ? sA : sB);
		return lab;
	}

	// ------------------------ Multi-label (|Σ|>=3) via α-expansion with multi-start ------------------------

	/**
	 * Run α-expansion from a deterministic greedy starting point plus several random starts;
	 * return the best result found.
	 *
	 * @param siteIndex used only for logging / deterministic seed mixing
	 */
	private Map<Node, Character> solveSiteAlphaExpansionMultiStart(
			Map<Node, Integer> allowedMasks, Set<Character> sigma, Universe universe, int siteIndex,
			List<Node> solveNodes, List<Edge> solveEdges) {

		// deterministic greedy start: each node takes its first allowed state
		Map<Node, Character> bestLabels = solveSiteAlphaExpansion(
				allowedMasks, sigma, universe, greedyInit(allowedMasks, sigma, universe), solveNodes, solveEdges);
		double bestEnergy = energyPotts(bestLabels, solveEdges);

		// Random restarts. Mix seed with siteIndex so different sites get different (but
		// deterministic) random sequences.
		var rng = new Random(seed ^ ((long) siteIndex * 0x9E3779B97F4A7C15L));
		for (int trial = 0; trial < NUM_RANDOM_RESTARTS; trial++) {
			// randomInit draws for every node, in sorted order, so that the random sequence does not depend on which
			// nodes are pinned
			var initLabels = randomInit(allowedMasks, universe, rng);
			var trialLabels = solveSiteAlphaExpansion(allowedMasks, sigma, universe, initLabels, solveNodes, solveEdges);
			double trialEnergy = energyPotts(trialLabels, solveEdges);
			if (trialEnergy < bestEnergy - 1e-9) {
				bestEnergy = trialEnergy;
				bestLabels = trialLabels;
			}
		}
		return bestLabels;
	}

	/**
	 * deterministic initialisation: the first allowed state of each node (in universe order)
	 */
	private Map<Node, Character> greedyInit(
			Map<Node, Integer> allowedMasks, Set<Character> sigma, Universe universe) {
		Map<Node, Character> label = new HashMap<>();
		char defaultBase = sigma.iterator().next();
		int allMask = universe.allMask();
		for (Node x : nodes) {
			int mask = allowedMasks.getOrDefault(x, allMask);
			label.put(x, (mask == allMask) ? defaultBase : firstAllowed(mask, universe));
		}
		return label;
	}

	/**
	 * Random feasible initialisation: each node gets a uniformly-random base from its allowed mask.
	 */
	private Map<Node, Character> randomInit(
			Map<Node, Integer> allowedMasks, Universe universe, Random rng) {
		Map<Node, Character> label = new HashMap<>();
		for (Node x : nodes) {
			int mask = allowedMasks.getOrDefault(x, universe.allMask());
			label.put(x, randomAllowed(mask, universe, rng));
		}
		return label;
	}

	/**
	 * α-expansion from a given initial labeling. Returns the converged labels.
	 *
	 * @param initLabels starting labels — must be feasible (in allowed mask) for every node
	 */
	private Map<Node, Character> solveSiteAlphaExpansion(
			Map<Node, Integer> allowedMasks, Set<Character> sigma, Universe universe,
			Map<Node, Character> initLabels, List<Node> solveNodes, List<Edge> solveEdges) {

		final double INF = 1e12;
		Map<Node, Character> label = new HashMap<>(initLabels);

		boolean improved;
		int iter = 0;
		do {
			improved = false;

			for (char alpha : sigma) {
				FlowNetwork<Node> FN = new FlowNetwork<>();
				for (Node x : solveNodes) FN.addNode(x);
				FN.addSourceSink();

				// Unary feasibility constraints for this move
				for (Node x : solveNodes) {
					int mask = allowedMasks.getOrDefault(x, universe.allMask());
					char cur = label.get(x);

					double D0 = allowed(mask, cur, universe) ? 0.0 : INF;   // keeping cur must be allowed
					double D1 = allowed(mask, alpha, universe) ? 0.0 : INF; // switching to alpha must be allowed

					if (D0 > 0) FN.addCapacity(FN.S, x, D0); // cost if y=0 (T-side)
					if (D1 > 0) FN.addCapacity(x, FN.T, D1); // cost if y=1 (S-side)
				}

				// Pairwise Potts terms for α-expansion
				for (Edge e : solveEdges) {
					Node p = u.apply(e), q = v.apply(e);
					char Lp = label.get(p);
					char Lq = label.get(q);

					double we = edgeWeight.applyAsDouble(e);
					double a = (Lp != Lq) ? we : 0;        // V(0,0)
					double b = (Lp != alpha) ? we : 0;     // V(0,1)
					double c = (Lq != alpha) ? we : 0;     // V(1,0)
					double d = 0;                          // V(1,1)

					addSubmodularPairwise(FN, p, q, a, b, c, d);
				}

				double before = energyPotts(label, solveEdges);
				FN.maxFlow();
				Set<Node> Sside = FN.minCutSourceSide();

				Map<Node, Character> candidate = new HashMap<>(label);
				for (Node x : solveNodes) {
					if (Sside.contains(x)) candidate.put(x, alpha);
				}

				double after = energyPotts(candidate, solveEdges);
				if (after + 1e-9 < before) {
					label = candidate;
					improved = true;
				}
			}

			iter++;
		} while (improved && iter < MAX_ALPHA_ITERATIONS);

		if (verbose && iter >= MAX_ALPHA_ITERATIONS && improved) {
			System.err.printf("ParsimonyLabeler: α-expansion hit iteration cap (%d); result may be suboptimal.%n",
					MAX_ALPHA_ITERATIONS);
		}

		// Defensive feasibility cleanup
		for (Node x : nodes) {
			int mask = allowedMasks.getOrDefault(x, universe.allMask());
			if (!allowed(mask, label.get(x), universe))
				label.put(x, firstAllowed(mask, universe));
		}
		return label;
	}

	/**
	 * the Potts energy over the given edges: the total weight of the edges whose endpoints have different labels
	 */
	private double energyPotts(Map<Node, Character> lab, List<Edge> edgesToScore) {
		double E = 0;
		for (Edge e : edgesToScore) {
			Node a = u.apply(e), b = v.apply(e);
			if (!Objects.equals(lab.get(a), lab.get(b))) E += edgeWeight.applyAsDouble(e);
		}
		return E;
	}

	// Add a submodular binary pairwise term V with values a=V(0,0), b=V(0,1), c=V(1,0), d=V(1,1)
	private void addSubmodularPairwise(FlowNetwork<Node> FN, Node p, Node q,
									   double a, double b, double c, double d) {
		if (a + d > b + c + 1e-9) throw new IllegalArgumentException("Pairwise term not submodular");

		// V(p,q) = const + (c-a)·p + (b-a)·q − w·p·q,   w = b+c-a-d ≥ 0
		//
		// Realise −w·p·q with a single directed edge p→q of capacity w.
		// That edge contributes w·p·(1−q) = w·p − w·p·q to the cut, so we
		// subtract w from the unary coefficient on p:
		//
		//   unary on p:  (c-a) − w = d − b
		//   unary on q:  (b-a)

		double w = b + c - a - d;
		if (w > 1e-12) FN.addCapacity(p, q, w);   // directed, NOT undirected

		double up = d - b;
		double uq = b - a;

		if (up > 0) FN.addCapacity(p, FN.T, up);
		else if (up < 0) FN.addCapacity(FN.S, p, -up);

		if (uq > 0) FN.addCapacity(q, FN.T, uq);
		else if (uq < 0) FN.addCapacity(FN.S, q, -uq);
	}

	// ------------------------ maximum flow (Dinic's algorithm) ------------------------

	private static class FlowNetwork<N> {
		static class EdgeRec<N> {
			N to;
			int rev;
			double cap;

			EdgeRec(N to, int rev, double cap) {
				this.to = to;
				this.rev = rev;
				this.cap = cap;
			}
		}

		private final Map<N, List<EdgeRec<N>>> adj = new HashMap<>();
		N S, T;

		void addNode(N x) {
			adj.computeIfAbsent(x, k -> new ArrayList<>());
		}

		@SuppressWarnings("unchecked")
		void addSourceSink() {
			this.S = (N) new Object();
			this.T = (N) new Object();
			addNode(S);
			addNode(T);
		}

		void addCapacity(N a, N b, double cap) {
			addNode(a);
			addNode(b);
			List<EdgeRec<N>> A = adj.get(a), B = adj.get(b);
			A.add(new EdgeRec<>(b, B.size(), cap));
			B.add(new EdgeRec<>(a, A.size() - 1, 0));
		}

		void addUndirectedCapacity(N a, N b, double cap) {
			addCapacity(a, b, cap);
			addCapacity(b, a, cap);
		}

		double maxFlow() {
			double flow = 0;
			Map<N, Integer> level = new HashMap<>();
			Map<N, Integer> it = new HashMap<>();

			while (bfs(level)) {
				it.clear();
				for (N x : adj.keySet()) it.put(x, 0);
				double f;
				while ((f = dfs(S, T, Double.POSITIVE_INFINITY, level, it)) > 1e-12) flow += f;
			}
			return flow;
		}

		private boolean bfs(Map<N, Integer> level) {
			level.clear();
			Deque<N> dq = new ArrayDeque<>();
			level.put(S, 0);
			dq.add(S);

			while (!dq.isEmpty()) {
				N x = dq.removeFirst();
				for (EdgeRec<N> e : adj.get(x)) {
					if (e.cap > 1e-12 && !level.containsKey(e.to)) {
						level.put(e.to, level.get(x) + 1);
						dq.addLast(e.to);
					}
				}
			}
			return level.containsKey(T);
		}

		private double dfs(N x, N t, double f, Map<N, Integer> level, Map<N, Integer> it) {
			if (x.equals(t)) return f;

			List<EdgeRec<N>> L = adj.get(x);
			for (int i = it.get(x); i < L.size(); i++, it.put(x, i)) {
				EdgeRec<N> e = L.get(i);
				if (e.cap <= 1e-12) continue;
				if (!level.containsKey(e.to) || level.get(e.to) != level.get(x) + 1) continue;

				double pushed = dfs(e.to, t, Math.min(f, e.cap), level, it);
				if (pushed > 1e-12) {
					e.cap -= pushed;
					EdgeRec<N> r = adj.get(e.to).get(e.rev);
					r.cap += pushed;
					return pushed;
				}
			}
			return 0;
		}

		Set<N> minCutSourceSide() {
			Set<N> vis = new HashSet<>();
			Deque<N> dq = new ArrayDeque<>();
			dq.add(S);
			vis.add(S);

			while (!dq.isEmpty()) {
				N x = dq.removeFirst();
				for (EdgeRec<N> e : adj.get(x)) {
					if (e.cap > 1e-12 && !vis.contains(e.to)) {
						vis.add(e.to);
						dq.addLast(e.to);
					}
				}
			}
			return vis;
		}
	}
}
