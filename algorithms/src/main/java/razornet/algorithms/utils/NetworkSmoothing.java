/*
 * NetworkSmoothing.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.CanceledException;
import jloda.util.progress.ProgressListener;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.function.IntPredicate;

/**
 * Smoothing of a network that realizes a metric: unlabeled nodes of degree at most two are suppressed (a node of
 * degree two is replaced by an edge between its two neighbours), and superfluous edges, whose endpoints are
 * joined by another path that is no longer than the edge, are removed. The two passes alternate until nothing
 * changes. Both preserve all distances between labeled nodes. If the time budget runs out, the graph as smoothed
 * so far is returned.
 *
 * @author Daniel Huson, 8.2026
 */
public class NetworkSmoothing {
	private static final int MAX_ROUNDS = 100;   // safety cap; the loop also stops as soon as nothing changes

	// the search for an alternative path expands at most this many nodes; beyond that, the edge is kept
	private static final int MAX_EXPANDED_NODES = 10000;

	/**
	 * what smoothing removed, for the run's report
	 */
	public record Result(int nodesRemoved, int edgesRemoved, boolean stoppedEarly) {
		public boolean isEmpty() {
			return nodesRemoved == 0 && edgesRemoved == 0;
		}
	}

	/**
	 * smooth the graph in place
	 *
	 * @param isLabeled       which node ids carry taxa and so may never be suppressed
	 * @param suppressNodes   run pass 1
	 * @param removeSuperfluousEdges run pass 2
	 */
	public static Result apply(IntGraph graph, IntPredicate isLabeled, boolean suppressNodes,
							   boolean removeSuperfluousEdges, ProgressListener progress) {
		var nodesRemoved = 0;
		var edgesRemoved = 0;
		var stoppedEarly = false;
		// every operation below preserves all distances, so on cancellation the graph as smoothed so far is valid
		try {
			for (var round = 0; round < MAX_ROUNDS; round++) {
				var nodes = (suppressNodes ? suppressUnlabeledNodes(graph, isLabeled, progress) : 0);
				// the two passes alternate, as each can enable the other; node suppression first, as it is the cheaper
				var edges = (removeSuperfluousEdges ? removeSuperfluousEdges(graph, progress) : 0);
				nodesRemoved += nodes;
				edgesRemoved += edges;
				if (nodes == 0 && edges == 0)
					break;
			}
		} catch (CanceledException ex) {
			stoppedEarly = true;
		}
		return new Result(nodesRemoved, edgesRemoved, stoppedEarly);
	}

	/**
	 * suppress unlabeled nodes of degree &le; 2, repeatedly, since suppressing one can drop another to degree 2
	 *
	 * @return the number of nodes removed
	 */
	private static int suppressUnlabeledNodes(IntGraph graph, IntPredicate isLabeled, ProgressListener progress) throws CanceledException {
		var removed = 0;
		for (var changed = true; changed; ) {
			changed = false;
			for (var u : new ArrayList<>(graph.nodes())) {
				progress.checkForCancel();
				if (isLabeled.test(u) || !graph.containsNode(u))
					continue;
				var degree = graph.getDegree(u);
				if (degree > 2)
					continue;
				if (degree == 2) {
					var edges = graph.getAdjacentEdges(u);
					var e1 = edges.get(0);
					var e2 = edges.get(1);
					var v1 = e1.other(u);
					var v2 = e2.other(u);
					var weight = graph.getWeight(e1) + graph.getWeight(e2);
					// if v1-v2 is already an edge, it keeps the smaller of the two weights
					if (v1 != v2) {
						if (!graph.hasEdge(v1, v2))
							graph.addEdge(v1, v2, weight);
						else if (graph.getWeight(v1, v2) > weight)
							graph.setWeight(v1, v2, weight);
					}
				}
				graph.removeNode(u);
				removed++;
				changed = true;
			}
		}
		return removed;
	}

	/**
	 * remove every edge whose endpoints are joined by another path of no greater length
	 *
	 * @return the number of edges removed
	 */
	private static int removeSuperfluousEdges(IntGraph graph, ProgressListener progress) throws CanceledException {
		var toRemove = new ArrayList<IntGraph.Edge>();
		try {
			for (var e : graph.edges()) {
				progress.checkForCancel();
				// an edge can only be superfluous if both of its endpoints have another edge
				if (graph.getDegree(e.u()) > 1 && graph.getDegree(e.v()) > 1
					&& hasAlternativePath(graph, e.u(), e.v(), graph.getWeight(e)))
					toRemove.add(e);
			}
		} finally {
			// on cancellation, the edges found so far are still removed
			if (!toRemove.isEmpty())
				graph.removeEdges(toRemove);
		}
		return toRemove.size();
	}

	/**
	 * is there a path from {@code from} to {@code to} of length at most {@code limit} that does not use the edge
	 * between them? Uses a Dijkstra search bounded by {@code limit} and by MAX_EXPANDED_NODES.
	 */
	private static boolean hasAlternativePath(IntGraph graph, int from, int to, int limit) {
		var distance = new HashMap<Integer, Integer>();
		var queue = new PriorityQueue<int[]>(Comparator.comparingInt(a -> a[1]));
		queue.add(new int[]{from, 0});
		while (!queue.isEmpty() && distance.size() < MAX_EXPANDED_NODES) {
			var top = queue.poll();
			var u = top[0];
			if (distance.containsKey(u))
				continue;
			distance.put(u, top[1]);
			if (u == to)
				return true;
			for (var v : graph.getAdjacentNodes(u)) {
				if ((u == from && v == to) || (u == to && v == from))
					continue;   // the edge under test
				if (distance.containsKey(v))
					continue;
				var d = top[1] + graph.getWeight(u, v);
				if (d <= limit)
					queue.add(new int[]{v, d});
			}
		}
		return false;
	}
}
