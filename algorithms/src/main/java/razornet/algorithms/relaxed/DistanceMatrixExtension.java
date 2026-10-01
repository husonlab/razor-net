/*
 * DistanceMatrixExtension.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.relaxed;

import jloda.util.CanceledException;
import jloda.util.SetUtils;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.utils.ExtensionTrace;
import razornet.algorithms.utils.IntGraph;
import razornet.algorithms.utils.MutableD;

import java.util.*;

/**
 * The distance-matrix extension of the Relaxed strategy. Each call compactifies one point of the current subset
 * S: the point of largest positive slack whose new point does not coincide with an existing one. The positive
 * slack of x in S is the smallest positive value of {@code (D(x,y)+D(x,z)-D(y,z))/2} over the pairs y,z of other
 * points of S, so a point that lies between some pairs but not others is compactified too. The other points of S
 * stay in the subset, a vertex is settled only if it has degree at most 2 and lies between two others, and every
 * call breaks ties by a new random order. See {@code docs/algorithm.md}.
 * <p>
 * Daniel Huson and Momoko Hayamizu, 2026
 */
public final class DistanceMatrixExtension {
	// the working matrix is not extended beyond max(MAX_EXTENSION_FACTOR * n, MIN_EXTENSION_CAP) points
	private static final int MAX_EXTENSION_FACTOR = 20;
	private static final int MIN_EXTENSION_CAP = 2000;

	/**
	 * extends the distance matrix by auxiliary points
	 *
	 * @param distances    the distances, diagonal values 0, all off-diagonal values positive even integers
	 * @param random       used to break ties; null for a deterministic run
	 * @param progress     report progress
	 * @param extensionCap the matrix is not extended beyond this size
	 * @throws CanceledException user canceled
	 */
	public static void extend(MutableD distances, Random random, ProgressListener progress, int extensionCap) throws CanceledException {
		var all = SetUtils.fullIndexSet(distances.size());
		extendSubset(distances, all, random, progress, extensionCap, ExtensionTrace.claimTracedRestart(), 0);
	}

	/**
	 * the cap on the size of the working matrix for a whole restart, computed once from the number of taxa
	 */
	public static int extensionCap(int groundSize) {
		return Math.max(MAX_EXTENSION_FACTOR * groundSize, MIN_EXTENSION_CAP);
	}

	/**
	 * recursive extension of the distance matrix on a subset of its points
	 */
	private static void extendSubset(MutableD distances, Set<Integer> subset, Random random, ProgressListener progress,
									 int extensionCap, boolean tracing, int depth) throws CanceledException {
		if (distances.size() >= extensionCap) {
			if (tracing) ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "cap reached");
			return; // the matrix has reached the cap
		}

		// ----- Step 1: compute the slacks and compactify the best point -----

		var order = new ArrayList<>(subset);

		if (random != null) { // use random ordering
			Collections.shuffle(order, random);
		}

		var choice = pickBestNonCoincidingSlack(distances, subset, order, progress);
		if (choice == null) {
			if (tracing)
				ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "no usable slack point");
			return;
		}

		var x = choice.x();
		var slack = choice.slack();

		progress.checkForCancel();

		var newVertices = new ArrayList<Integer>();    // the new point added in this call

		var column = buildAuxColumnUsingMaxAll(distances, x, slack); // the distances of the new point xp

		// Append xp to the matrix.
		var xp = distances.appendVertex(column);
		newVertices.add(xp);
		progress.checkForCancel();


		// ----- the temporary graph G on the new point and the other points of S -----
		var verticesG = new ArrayList<>(newVertices);
		subset.stream().filter(z -> z != x).forEach(verticesG::add);

		var graph = new IntGraph(verticesG);  // simple undirected graph
		for (var a : verticesG) {
			for (var b : verticesG) {
				if (a < b && isEssentialEdge(distances, a, b))
					graph.addEdge(a, b, distances.get(a, b));
			}
		}

		// ----- confirm edges (getDegree ≤ 2) and recurse on unconfirmed connectedComponents -----
		progress.checkForCancel();

		// a vertex is settled if it has degree at most 2 in G and lies between two points of G; an edge is confirmed
		// if one of its endpoints is settled
		var settled = new HashSet<Integer>();
		for (var v : verticesG)
			if (graph.getDegree(v) <= 2 && hasZeroSlack(distances, verticesG, v))
				settled.add(v);

		var unconfirmed = new HashSet<IntGraph.Edge>();
		for (var e : graph.edges()) {
			boolean confirmed = settled.contains(e.u()) || settled.contains(e.v());
			if (!confirmed) unconfirmed.add(e);
		}

		if (unconfirmed.isEmpty()) {
			if (tracing)
				ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "all edges confirmed");
			return; // All edges in this local region are confirmed → stop here.
		}

		progress.checkForCancel();

		var H = graph.inducedByEdges(unconfirmed);// subgraph on unconfirmed edges

		if (tracing && ExtensionTrace.isEnabled()) {
			var sizes = new StringBuilder();
			for (var c : H.connectedComponents()) sizes.append(sizes.isEmpty() ? "" : ",").append(c.size());
			var keptFromS = 0;
			for (var v : verticesG) if (v != xp) keptFromS++;
			ExtensionTrace.recordCall(depth, subset.size(), 1, "c+(x)=" + slack.c(), 1, 0, verticesG.size(),
					keptFromS, graph.getNumberOfEdges(), settled.size(), unconfirmed.size(), sizes.toString(), distances.size());
		}
		ExtensionTrace.noteMatrixSize(distances.size());

		for (var comp : H.connectedComponents()) { // Recurse on this connected set of nodes (indices into current M)
			extendSubset(distances, comp, random, progress, extensionCap, tracing, depth + 1);
		}
	}

	public record Slack(int c, int y, int z) {
	}

	// Static counters for diagnostic logging.
	public static int totalCalls = 0;
	public static int tieCount = 0;          // calls where 2+ candidates shared the max slack value
	public static int totalTieSize = 0;      // sum of tie sizes (counts only ties, i.e., ≥2)
	public static int coincidenceSkips = 0;  // calls where best-slack candidate had to be skipped due to coincidence

	public static void resetCounters() {
		totalCalls = 0;
		tieCount = 0;
		totalTieSize = 0;
		coincidenceSkips = 0;
	}

	public static String summary() {
		if (totalCalls == 0) return "(no calls)";
		return String.format(
				"calls=%d, ties=%d (%.1f%% of calls), avg tie size=%.2f, coincidence skips=%d (%.1f%% of calls)",
				totalCalls,
				tieCount, 100.0 * tieCount / totalCalls,
				tieCount > 0 ? (double) totalTieSize / tieCount : 0,
				coincidenceSkips, 100.0 * coincidenceSkips / totalCalls);
	}

	private static SlackChoice pickBestNonCoincidingSlack(MutableD distances, Set<Integer> subset, List<Integer> order, ProgressListener progress) throws CanceledException {
		totalCalls++;

		// Collect all positive slacks.
		var candidates = new ArrayList<SlackChoice>();
		for (var x : order) {
			progress.checkForCancel();
			var s = slackWithArgmin(distances, subset, x);
			if (s.c() > 0) candidates.add(new SlackChoice(x, s));
		}
		if (candidates.isEmpty()) return null;


		// sort by slack desc only; ties broken by encounter order
		// (which depends on `order`, typically randomised upstream).
		candidates.sort((a, b) -> Integer.compare(b.slack().c(), a.slack().c()));


		// Track tie size at the top.
		var topValue = candidates.get(0).slack().c();
		var tiedAtTop = 0;
		for (var c : candidates) {
			if (c.slack().c() == topValue) tiedAtTop++;
			else break;
		}
		if (tiedAtTop > 1) {
			tieCount++;
			totalTieSize += tiedAtTop;
		}
		// Find the highest-slack one whose aux node doesn't coincide.
		boolean firstSkipped = false;
		try {
			for (var c : candidates) {
				var col = buildAuxColumnUsingMaxAll(distances, c.x(), c.slack());
				if (findCoincidingVertex(distances, col) == null) {
					return c;
				}
				firstSkipped = true;
			}
		} finally {
			if (firstSkipped) coincidenceSkips++;
		}
		return null;
	}

	private record SlackChoice(int x, Slack slack) {
	}

	/**
	 * the positive slack of x within {@code subset}: the smallest positive value of
	 * {@code (D(x,y)+D(x,z)-D(y,z))/2} over all pairs y,z of other points of the subset, or 0 if no pair gives a
	 * positive value, together with the first pair that attains it
	 */
	public static Slack slackWithArgmin(MutableD distances, Set<Integer> subset, int x) {
		var bc = Integer.MAX_VALUE;
		int by = 0;
		int bz = 0;
		for (int y : subset) {
			if (y != x) {
				for (int z : subset) {
					if (z != x && z != y) {
						var sum = distances.get(x, y) + distances.get(x, z) - distances.get(y, z);
						if (sum % 2 != 0)
							throw new IllegalStateException("slackWithArgmin: Sum of distance values is odd");
						var c = sum / 2;
						if (c > 0 && c < bc) {
							bc = c;
							by = y;
							bz = z;
						}
					}
				}
			}
		}
		if (bc == Integer.MAX_VALUE)
			return new Slack(0, 0, 0);
		else return new Slack(bc, by, bz);
	}

	/**
	 * returns coinciding vertex, if there is one
	 *
	 * @param distances distances
	 * @param cand      the candidate new distances
	 * @return vertex that has distance zero from new vertex
	 */
	public static Integer findCoincidingVertex(MutableD distances, int[] cand) {
		int n = distances.size();
		for (int a = 0; a < n; a++) {
			if (cand[a] == 0)
				return a;
		}
		return null;
	}

	/**
	 * the distances from the compactification point x' of x to all points: c to x, D(x,y)-c and D(x,z)-c to the
	 * minimizer pair y,z and, to each further point in index order, the smallest value that satisfies the triangle
	 * inequality with the points already placed
	 *
	 * @param distances distance matrix
	 * @param x         the point to compactify
	 * @param slack     its slack c and minimizer pair y,z
	 * @return the distances of the new point to all points of the matrix
	 */
	public static int[] buildAuxColumnUsingMaxAll(MutableD distances, int x, Slack slack) {

		var c = slack.c();
		var y = slack.y();
		var z = slack.z();

		var n = distances.size();
		var column = new int[n];
		var inS = new boolean[n];

		column[x] = c;
		inS[x] = true;
		column[y] = Math.max(distances.get(y, x) - c, 0);
		inS[y] = true;
		column[z] = Math.max(distances.get(z, x) - c, 0);
		inS[z] = true;

		for (var a = 0; a < n; a++) {
			if (!inS[a]) {
				var best = Integer.MIN_VALUE;
				for (var s = 0; s < n; s++) {
					if (inS[s]) {
						var alt = distances.get(a, s) - column[s];
						if (alt > best)
							best = alt;
					}
				}
				column[a] = Math.max(best, 0);
				inS[a] = true;
			}
		}

		return column;
	}

	/**
	 * Does v lie between two points of {@code verts}, or does it have fewer than two other points there? Such a
	 * vertex needs no auxiliary neighbour. This is the test for slack zero, not the positive slack of
	 * {@link #slackWithArgmin}.
	 */
	private static boolean hasZeroSlack(MutableD distances, Collection<Integer> verts, int v) {
		var sawPair = false;
		for (var y : verts) {
			if (y == v) continue;
			for (var z : verts) {
				if (z == v || z <= y) continue;
				sawPair = true;
				if (distances.get(v, y) + distances.get(v, z) - distances.get(y, z) == 0)
					return true;
			}
		}
		return !sawPair;
	}

	/**
	 * is {x,y} an essential pair, i.e. is there no third point z, at positive distance from both, with
	 * d(x,y) >= d(x,z) + d(z,y)?
	 */
	public static boolean isEssentialEdge(MutableD distances, int x, int y) {
		if (x == y) return false;
		var dxy = distances.get(x, y);
		for (var z = 0; z < distances.size(); z++) {
			if (z == x || z == y) continue;
			var dxz = distances.get(x, z);
			var dzy = distances.get(z, y);
			if (dxz == 0 || dzy == 0) continue;
			if (dxy >= dxz + dzy) return false;
		}
		return true;
	}
}
