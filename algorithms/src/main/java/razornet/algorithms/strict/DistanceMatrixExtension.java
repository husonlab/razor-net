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

package razornet.algorithms.strict;

import jloda.util.CanceledException;
import jloda.util.SetUtils;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.utils.MutableD;
import razornet.algorithms.utils.ExtensionTrace;
import razornet.algorithms.utils.IntGraph;

import java.util.*;

/**
 * The distance-matrix extension of the Strict strategy. At every level of the recursion, every point of the
 * current subset S that lies on no geodesic between two other points of S (a slack point) is compactified, in
 * decreasing order of slack, against the matrix as it grows; the recursion continues on the new points and the
 * points of S that are not slack points. The slack of x in S is the minimum of
 * {@code (D(x,y)+D(x,z)-D(y,z))/2} over the pairs y,z of other points of S. See {@code docs/algorithm.md}.
 * <p>
 * Momoko Hayamizu and Daniel Huson, 2026
 */
public final class DistanceMatrixExtension {
	private static final int MAX_EXTENSION_FACTOR = 5;
	private static final int MIN_EXTENSION_CAP = 500;

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
		var tracing = ExtensionTrace.claimTracedRestart(); // exactly one restart records its recursion
		extendSubset(distances, all, random, progress, extensionCap, tracing, 0);
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
	private static void extendSubset(MutableD distances, Set<Integer> subset, Random random, ProgressListener progress, int extensionCap,
									 boolean tracing, int depth) throws CanceledException {
		if (distances.size() >= extensionCap) {
			// stopping condition 0: the matrix has reached the cap
			capHits++;
			if (tracing) ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "cap reached (stop 0)");
			return;
		}

		// ----- Step 1: slack -----

		// a random order of the subset, used to break ties among equal slacks (the sort below is stable)
		var order = new ArrayList<>(subset);

		if (random != null) { // use random ordering
			Collections.shuffle(order, random);
		}

		var slacklist = buildSlacklist(distances, subset, order, random, progress);
		if (slacklist.isEmpty()) {
			if (tracing) ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "no slack point (stop 1)");
			return; // stopping condition 1: no slack point in S
		}

		// The auxiliary points added below take the consecutive indices nBefore, nBefore+1, ..., |D|-1.
		var nBefore = distances.size();
		var coincidenceSkipsHere = 0;  // this call only; the static counter is cumulative and shared across threads

		// compactify every slack point, in decreasing order of slack
		for (var slackPoint : slacklist) {
			progress.checkForCancel();
			// the column is built against the current matrix, which grows in this loop; the slack and the minimizer
			// pair are those computed before the loop
			var column = buildAuxColumnUsingMaxAll(distances, slackPoint.x(), slackPoint.slack());
			if (findCoincidingVertex(distances, column) != null) {
				// the new point would coincide with an existing point: skip it
				coincidenceSkips++;
				coincidenceSkipsHere++;
				continue;
			}
			distances.appendVertex(column);
			pointsAdded++;
		}
		noteMaxMatrixSize(distances.size());
		ExtensionTrace.noteMatrixSize(distances.size());

		if (distances.size() == nBefore) {
			if (tracing) ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "all coincided (stop 1')");
			return; // stopping condition 1': every slack point coincided, so nothing was added
		}

		// ----- the temporary graph G on the new points and the points of S that are not slack points -----
		var verticesG = new ArrayList<Integer>();
		for (var v = nBefore; v < distances.size(); v++)
			verticesG.add(v);

		{
			var slackPoints = new HashSet<Integer>();
			for (var slackPoint : slacklist)
				slackPoints.add(slackPoint.x());
			subset.stream().filter(z -> !slackPoints.contains(z)).forEach(verticesG::add);
		}

		var graph = new IntGraph(verticesG);  // simple undirected graph
		for (var a : verticesG) {
			for (var b : verticesG) {
				if (a < b && isEssentialEdge(distances, a, b))
					graph.addEdge(a, b, distances.get(a, b));
			}
		}

		// ----- confirm edges (getDegree ≤ 2) and recurse on unconfirmed connectedComponents -----
		progress.checkForCancel();

		// a vertex of degree at most 2 in G is settled; an edge is confirmed if one of its endpoints is settled
		var settled = new HashSet<Integer>();
		for (var v : verticesG)
			if (graph.getDegree(v) <= 2)
				settled.add(v);

		var unconfirmed = new HashSet<IntGraph.Edge>();
		for (var e : graph.edges()) {
			boolean confirmed = settled.contains(e.u()) || settled.contains(e.v());
			if (!confirmed) unconfirmed.add(e);
		}

		if (unconfirmed.isEmpty()) {
			if (tracing) ExtensionTrace.recordEarlyReturn(depth, subset.size(), distances.size(), "all edges confirmed (stop 2)");
			return; // All edges in this local region are confirmed → stop here.
		}

		progress.checkForCancel();

		var H = graph.inducedByEdges(unconfirmed);// subgraph on unconfirmed edges

		if (tracing && ExtensionTrace.isEnabled()) {
			var sizes = new StringBuilder();
			for (var c : H.connectedComponents()) sizes.append(sizes.isEmpty() ? "" : ",").append(c.size());
			var cs = new StringBuilder();
			for (var i = 0; i < Math.min(10, slacklist.size()); i++)
				cs.append(i == 0 ? "" : ",").append(slacklist.get(i).slack().c());
			if (slacklist.size() > 10) cs.append(",...");
			var keptFromS = 0;
			for (var v : verticesG) if (v < nBefore) keptFromS++;
			ExtensionTrace.recordCall(depth, subset.size(), slacklist.size(), cs.toString(),
					distances.size() - nBefore, coincidenceSkipsHere, verticesG.size(), keptFromS,
					graph.getNumberOfEdges(), settled.size(), unconfirmed.size(), sizes.toString(), distances.size());
		}
		for (var comp : H.connectedComponents()) { // Recurse on this connected set of nodes (indices into current M)
			extendSubset(distances, comp, random, progress, extensionCap, tracing, depth + 1);
		}
	}

	public record Slack(int c, int y, int z) {
	}

	// counters for diagnostic logging, reported by summary(); they are not synchronized, so the counts are
	// approximate when restarts run in parallel
	public static int totalCalls = 0;        // calls that reached the slack scan, i.e. passed stopping condition 0
	public static int tieCount = 0;          // calls where 2+ slack points shared the top slack value
	public static int totalTieSize = 0;      // sum of tie sizes (counts only ties, i.e., >= 2)
	public static int coincidenceSkips = 0;  // slack points skipped because their new point coincided with an existing one
	public static int pointsAdded = 0;       // auxiliary points actually appended
	public static int maxMatrixSize = 0;     // largest |D| seen
	public static int capHits = 0;           // calls that returned at stopping condition 0

	public static void resetCounters() {
		totalCalls = 0;
		tieCount = 0;
		totalTieSize = 0;
		coincidenceSkips = 0;
		pointsAdded = 0;
		maxMatrixSize = 0;
		capHits = 0;
	}

	private static void noteMaxMatrixSize(int size) {
		if (size > maxMatrixSize)
			maxMatrixSize = size;
	}

	public static String summary() {
		if (totalCalls == 0) return "(no calls)";
		return String.format(
				"calls=%d, ties=%d (%.1f%% of calls), avg tie size=%.2f, coincidence skips=%d, points added=%d, max |D|=%d, cap hits=%d",
				totalCalls,
				tieCount, 100.0 * tieCount / totalCalls,
				tieCount > 0 ? (double) totalTieSize / tieCount : 0,
				coincidenceSkips, pointsAdded, maxMatrixSize, capHits);
	}

	/**
	 * the slack points of {@code subset}, in decreasing order of slack; equal slacks keep the order given by
	 * {@code order} (the sort is stable)
	 *
	 * @param random used to choose among the minimizer pairs of a point, see
	 *               {@link #slackWithArgmin(MutableD, Set, int, Random)}
	 * @return the slack points; empty if no point of the subset has positive slack
	 */
	private static List<SlackPoint> buildSlacklist(MutableD distances, Set<Integer> subset, List<Integer> order, Random random, ProgressListener progress) throws CanceledException {
		totalCalls++;

		var slacklist = new ArrayList<SlackPoint>();
		for (var x : order) {
			progress.checkForCancel();
			var slack = slackWithArgmin(distances, subset, x, random);
			if (slack.c() > 0)
				slacklist.add(new SlackPoint(x, slack));
		}
		if (slacklist.isEmpty())
			return slacklist;

		slacklist.sort((a, b) -> Integer.compare(b.slack().c(), a.slack().c()));
		countTieAtTop(slacklist);
		return slacklist;
	}

	/**
	 * diagnostic: counts the slack points that share the largest slack value
	 */
	private static void countTieAtTop(List<SlackPoint> slacklist) {
		var topValue = slacklist.get(0).slack().c();
		var tiedAtTop = 0;
		for (var slackPoint : slacklist) {
			if (slackPoint.slack().c() == topValue) tiedAtTop++;
			else break;
		}
		if (tiedAtTop > 1) {
			tieCount++;
			totalTieSize += tiedAtTop;
		}
	}

	/**
	 * a slack point x with its slack and a minimizer pair
	 */
	private record SlackPoint(int x, Slack slack) {
	}

	/**
	 * the slack of x within {@code subset}, with the first minimizer pair in index order
	 */
	public static Slack slackWithArgmin(MutableD distances, Set<Integer> subset, int x) {
		return slackWithArgmin(distances, subset, x, null);
	}

	/**
	 * the slack of x within {@code subset}: the minimum of {@code (D(x,y)+D(x,z)-D(y,z))/2} over all pairs y,z of
	 * other points of the subset, or 0 if there are fewer than two other points, together with a pair that attains
	 * the minimum. A point with positive slack is a slack point; slack 0 means that x lies on a geodesic between
	 * two points of the subset.
	 *
	 * @param random if not null, the minimizer pair is chosen uniformly at random among all pairs that attain the
	 *               minimum (by reservoir sampling); if null, the first such pair in index order is taken
	 */
	public static Slack slackWithArgmin(MutableD distances, Set<Integer> subset, int x, Random random) {
		var bc = Integer.MAX_VALUE;
		int by = 0;
		int bz = 0;
		var attaining = 0; // how many pairs have attained bc, for uniform reservoir sampling
		var sawPair = false;
		for (int y : subset) {
			if (y != x) {
				for (int z : subset) {
					if (z != x && z != y) {
						sawPair = true;
						var sum = distances.get(x, y) + distances.get(x, z) - distances.get(y, z);
						if (sum % 2 != 0)
							throw new IllegalStateException("slackWithArgmin: Sum of distance values is odd");
						var c = sum / 2;
						if (c < 0)
							throw new IllegalStateException("slackWithArgmin: Triangle inequality violated at x=" + x + ", y=" + y + ", z=" + z);
						if (c < bc) {
							bc = c;
							by = y;
							bz = z;
							attaining = 1;
						} else if (c == bc && random != null) {
							// choose uniformly among the pairs that attain the minimum (reservoir sampling)
							attaining++;
							if (random.nextInt(attaining) == 0) {
								by = y;
								bz = z;
							}
						}
					}
				}
			}
		}
		if (!sawPair)
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
