/*
 * TCSConstruction.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.tcs;

import jloda.util.CanceledException;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.utils.IntGraph;

import java.util.*;

/**
 * The TCS construction (Templeton, Crandall &amp; Sing 1992): builds a haplotype network by joining components in
 * order of increasing distance, laying paths of unit steps through intermediate vertices. Candidate attachment
 * points are scored by how many distances between the two components the new path reproduces.
 * <p>
 * Based on the implementation in PopART (Leigh &amp; Bryant 2015, {@code src/networks/TCS.cpp}, GPL-3), modified
 * so that a join that adds nothing is skipped instead of raising an error. The input is an integer distance
 * matrix; the network need not realize it. The 95% parsimony connection limit is not implemented.
 *
 * @author Daniel Huson, 8.2026
 */
public class TCSConstruction {
	// the score of a candidate path: a bonus for each pair of taxa across the two components whose distance it
	// reproduces, a penalty for each pair that it puts too far apart, and a larger one for each pair too close
	private static final int BONUS = 20;
	private static final int SHORTCUT_PENALTY = 10;
	private static final int LONG_PENALTY = 5;

	private static final int NO_COMPONENT = -1;   // an intermediate vertex: in no component, reusable by both sides

	/**
	 * the result of a run
	 *
	 * @param graph         the network, taxa at node ids 0..n-1, edge weights in units of the input distance
	 * @param intermediates the number of intermediate vertices created (before smoothing)
	 * @param joinsSkipped  the number of joins that added nothing, because the two components were already at
	 *                      least as close as the distance being merged
	 */
	public record Result(IntGraph graph, int intermediates, int joinsSkipped) {
	}

	/**
	 * run the construction
	 *
	 * @param distances integer distance matrix, no two rows equal; not modified
	 */
	public static Result apply(int[][] distances, ProgressListener progress) throws CanceledException {
		var ntax = distances.length;
		var graph = new IntGraph();
		for (var i = 0; i < ntax; i++)
			graph.ensureNode(i);

		// component[v] for every vertex; intermediate vertices are in no component (NO_COMPONENT)
		var component = new ArrayList<Integer>();
		for (var i = 0; i < ntax; i++)
			component.add(i);

		// all pairs, grouped by distance; the groups are processed in ascending order
		var buckets = new TreeMap<Integer, List<int[]>>();
		for (var i = 0; i < ntax; i++)
			for (var j = 0; j < i; j++)
				buckets.computeIfAbsent(distances[i][j], k -> new ArrayList<>()).add(new int[]{j, i});

		progress.setMaximum(buckets.size());
		progress.setProgress(0);

		var intermediates = 0;
		var skipped = 0;

		for (var entry : buckets.entrySet()) {
			var distance = entry.getKey();
			var pending = entry.getValue();
			// one pass joins every pair between the same two components; the other pairs wait for the next pass
			while (!pending.isEmpty()) {
				progress.checkForCancel();
				var componentA = NO_COMPONENT;
				var componentB = NO_COMPONENT;
				var deferred = new ArrayList<int[]>();

				for (var pair : pending) {
					var u = pair[0];
					var v = pair[1];
					var cu = component.get(u);
					var cv = component.get(v);
					if (cu.equals(cv))
						continue;
					if (cu > cv) {
						var t = u; u = v; v = t;
						var c = cu; cu = cv; cv = c;
					}
					if (componentA == NO_COMPONENT) {
						componentA = cu;
						componentB = cv;
					}
					if (cu == componentA && cv == componentB) {
						var added = join(graph, component, distances, u, v, distance, componentA, componentB);
						if (added < 0)
							skipped++;
						else
							intermediates += added;
					} else
						deferred.add(pair);
				}

				if (componentA == NO_COMPONENT)
					break;   // nothing left to merge in this bucket

				for (var x = 0; x < component.size(); x++) {
					if (component.get(x).equals(componentB))
						component.set(x, componentA);
				}
				pending = deferred;
			}
			progress.incrementProgress();
		}
		return new Result(graph, intermediates, skipped);
	}

	/**
	 * join the component of {@code u} to the component of {@code v} with a path of {@code distance} unit steps
	 *
	 * @return the number of intermediate vertices created, or -1 if the join added nothing
	 */
	private static int join(IntGraph graph, List<Integer> component, int[][] distances, int u, int v,
							int distance, int componentA, int componentB) {
		if (distance <= 1) {
			if (!graph.hasEdge(u, v))
				graph.addEdge(u, v, 1);
			return 0;
		}

		var ntax = distances.length;
		// the distance of each vertex from u and from v (breadth-first, as all edges have weight 1)
		var fromU = breadthFirst(graph, u);
		var fromV = breadthFirst(graph, v);

		var bestScore = Integer.MIN_VALUE;
		var bestGap = Integer.MAX_VALUE;
		var bestI = u;
		var bestJ = v;

		// distances from each candidate to the taxa of its own component, cached per candidate
		var fromCandidate = new HashMap<Integer, int[]>();

		for (var i = 0; i < component.size(); i++) {
			var ci = component.get(i);
			if (ci != componentA && ci != NO_COMPONENT)
				continue;
			var du = fromU[i];
			if (du < 0 || du >= distance)
				continue;
			for (var j = 0; j < component.size(); j++) {
				var cj = component.get(j);
				if (cj != componentB && cj != NO_COMPONENT)
					continue;
				var dv = fromV[j];
				if (dv < 0 || du + dv >= distance)
					continue;
				var gap = distance - du - dv;

				var pathFromI = fromCandidate.computeIfAbsent(i, x -> breadthFirst(graph, x));
				var pathFromJ = fromCandidate.computeIfAbsent(j, x -> breadthFirst(graph, x));
				var score = score(component, distances, ntax, pathFromI, pathFromJ, componentA, componentB, gap, distance);
				if (score == Integer.MIN_VALUE)
					continue;   // inadmissible: it would put two taxa closer than the distance being merged
				if (score > bestScore || (score == bestScore && gap < bestGap)) {
					bestScore = score;
					bestGap = gap;
					bestI = i;
					bestJ = j;
				}
			}
		}

		if (bestScore == Integer.MIN_VALUE)
			bestGap = distance;   // no admissible attachment: fall back to joining u and v directly

		if (bestI == bestJ)
			return -1;

		// the two points are already at least as close as the path to be added: nothing to add
		var existing = breadthFirst(graph, bestI)[bestJ];
		if (existing >= 0 && existing <= bestGap)
			return -1;

		var previous = bestI;
		var created = 0;
		for (var step = 1; step < bestGap; step++) {
			// node ids are consecutive: component.size() is the next free id
			var w = component.size();
			graph.ensureNode(w);
			component.add(NO_COMPONENT);
			graph.addEdge(previous, w, 1);
			previous = w;
			created++;
		}
		if (!graph.hasEdge(previous, bestJ))
			graph.addEdge(previous, bestJ, 1);
		return created;
	}

	/**
	 * the score of a path of {@code gap} steps between two attachment points: how well it reproduces the distances
	 * between the taxa of the two components
	 *
	 * @return the score, or {@code Integer.MIN_VALUE} if the path would put two taxa closer than the distance
	 * being merged
	 */
	private static int score(List<Integer> component, int[][] distances, int ntax, int[] pathFromI, int[] pathFromJ,
							 int componentA, int componentB, int gap, int mergeDistance) {
		var score = 0;
		for (var a = 0; a < ntax; a++) {
			if (component.get(a) != componentA || pathFromI[a] < 0)
				continue;
			for (var b = 0; b < ntax; b++) {
				if (component.get(b) != componentB || pathFromJ[b] < 0)
					continue;
				var total = gap + pathFromI[a] + pathFromJ[b];
				if (total == distances[a][b])
					score += BONUS;
				else if (total > distances[a][b])
					score -= LONG_PENALTY;
				else if (total < mergeDistance)
					return Integer.MIN_VALUE;
				else
					score -= SHORTCUT_PENALTY;
			}
		}
		return score;
	}

	/**
	 * breadth-first distances from one vertex, in unit steps; -1 where unreachable
	 */
	private static int[] breadthFirst(IntGraph graph, int source) {
		var maximum = nextFreeNode(graph);
		var distance = new int[maximum];
		Arrays.fill(distance, -1);
		if (source >= maximum)
			return distance;
		distance[source] = 0;
		var queue = new ArrayDeque<Integer>();
		queue.add(source);
		while (!queue.isEmpty()) {
			var u = queue.poll();
			for (var v : graph.getAdjacentNodes(u)) {
				if (v < maximum && distance[v] < 0) {
					distance[v] = distance[u] + 1;
					queue.add(v);
				}
			}
		}
		return distance;
	}

	/**
	 * one past the highest node id in use, so ids stay dense and can index an array
	 */
	private static int nextFreeNode(IntGraph graph) {
		var max = -1;
		for (var v : graph.nodes())
			max = Math.max(max, v);
		return max + 1;
	}
}
