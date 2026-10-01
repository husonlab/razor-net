/*
 * GreedyStretchEdgeRemoval.java Copyright (C) 2026 Daniel H. Huson
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
import java.util.HashSet;
import java.util.Set;
import java.util.function.BiFunction;
import java.util.function.IntPredicate;

/**
 * Stretch filter: lossy post-processing that greedily removes edges, heaviest first, as long as no distance
 * between two labeled nodes exceeds {@code (1 + maxStretchFraction)} times its input value, so that the graph
 * stays a (1 + maxStretchFraction)-spanner of the input metric. This is the greedy spanner algorithm of Althöfer
 * et al. (1993), applied to a realization. With maxStretchFraction = 0 it is {@link GreedyGlobalEdgeRemoval}.
 */
public class GreedyStretchEdgeRemoval {

	/**
	 * Greedily remove edges while keeping every labeled-pair distance within {@code (1 + maxStretchFraction)}
	 * of the input distance.
	 *
	 * @param graph              undirected weighted graph that realizes {@code distance} on the labeled nodes
	 * @param isLabeled          is the node a taxon (an anchor for the distance constraints)?
	 * @param distance           target symmetric input distances over labeled nodes (same units as edge weights)
	 * @param maxStretchFraction maximum allowed relative inflation; e.g. 0.10 for 10%. 0 = exact (lossless).
	 * @param progress           progress listener (honours the wall-clock deadline)
	 * @return true, if the graph changed
	 */
	public static boolean apply(IntGraph graph, IntPredicate isLabeled, BiFunction<Integer, Integer, Integer> distance,
								double maxStretchFraction, ProgressListener progress) throws CanceledException {
		if (maxStretchFraction < 0)
			throw new IllegalArgumentException("maxStretchFraction must be >= 0");
		if (maxStretchFraction == 0) // exact removal
			return GreedyGlobalEdgeRemoval.apply(graph, isLabeled, distance, progress);

		var factor = 1.0 + maxStretchFraction;
		var nodes0 = graph.getNumberOfNodes();
		var edges0 = graph.getNumberOfEdges();

		// the edges kept; edges are removed from the graph at the end
		var alive = new HashSet<>(graph.edges());

		// heaviest edge first, then by node ids (deterministic)
		var order = new ArrayList<>(graph.edges());
		order.sort(Comparator
				.comparingInt((IntGraph.Edge e) -> -graph.getWeight(e))
				.thenComparingInt(e -> -Math.min(e.u(), e.v()))
				.thenComparingInt(e -> -Math.max(e.u(), e.v())));

		for (var e : order) {
			if (!alive.contains(e)) continue;
			var aliveWithoutE = new HashSet<>(alive);
			aliveWithoutE.remove(e);
			if (allLabeledPairsWithinStretch(graph, isLabeled, aliveWithoutE, distance, factor))
				alive.remove(e);
			progress.checkForCancel();
		}

		for (var e : new ArrayList<>(graph.edges())) {
			if (!alive.contains(e))
				graph.removeEdge(e);
			progress.checkForCancel();
		}

		var changed = graph.getNumberOfEdges() < edges0;
		if (changed)
			GraphCleanAndSmooth.apply(graph, isLabeled, progress); // dissolve degree-2 / drop dangling auxiliary vertices

		var achieved = maxAchievedStretch(graph, isLabeled, distance);
		System.err.printf("Stretch removal (allow <= %.2f%%): n=%4d e=%4d -> n=%4d e=%4d, max inflation used: %.2f%%%n",
				100 * maxStretchFraction, nodes0, edges0, graph.getNumberOfNodes(), graph.getNumberOfEdges(),
				100 * (achieved - 1.0));
		return changed;
	}

	/**
	 * True iff every labeled pair's shortest path over the {@code alive} edge set is within {@code factor} times
	 * its input distance. An unreachable pair (distance {@code Integer.MAX_VALUE}) fails the test, so removals
	 * that would disconnect taxa are rejected automatically.
	 */
	private static boolean allLabeledPairsWithinStretch(IntGraph graph, IntPredicate isLabeled, Set<IntGraph.Edge> alive,
														BiFunction<Integer, Integer, Integer> distance, double factor) {
		var labeled = graph.nodes().stream().filter(isLabeled::test).toList();
		var shortestPaths = DijkstraAllPairs.allShortestPathLengths(graph, labeled, graph.nodes(), alive);
		for (var u : labeled) {
			for (var v : labeled) {
				if (u < v) {
					var d = distance.apply(u, v);
					var got = shortestPaths.get(u).get(v); // Integer.MAX_VALUE if u,v disconnected over `alive`
					if (got > d * factor + 1e-9)
						return false;
				}
			}
		}
		return true;
	}

	/**
	 * Largest realized ratio {@code d_graph(u,v) / d_input(u,v)} over labeled pairs (1.0 = distances kept exact).
	 * Computed once on the final graph, only to report what fraction of the allowed budget was actually used.
	 */
	private static double maxAchievedStretch(IntGraph graph, IntPredicate isLabeled, BiFunction<Integer, Integer, Integer> distance) {
		var labeled = graph.nodes().stream().filter(isLabeled::test).toList();
		var shortestPaths = DijkstraAllPairs.allShortestPathLengths(graph, labeled, graph.nodes(), null); // null = all current edges
		var max = 1.0;
		for (var u : labeled) {
			for (var v : labeled) {
				if (u < v) {
					var d = distance.apply(u, v);
					var got = shortestPaths.get(u).get(v);
					if (d > 0 && got < Integer.MAX_VALUE)
						max = Math.max(max, (double) got / d);
				}
			}
		}
		return max;
	}
}
