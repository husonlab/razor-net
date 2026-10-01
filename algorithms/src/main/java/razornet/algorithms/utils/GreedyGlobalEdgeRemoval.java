/*
 * GreedyGlobalEdgeRemoval.java Copyright (C) 2026 Daniel H. Huson
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

import java.util.*;
import java.util.function.BiFunction;
import java.util.function.IntPredicate;

/**
 * Greedy edge removal: removes edges, heaviest first, as long as all distances between labeled nodes are
 * preserved; unlabeled nodes of degree at most two are then removed ({@link GraphCleanAndSmooth}).
 */
public class GreedyGlobalEdgeRemoval {
	public static boolean verbose = false;

	/**
	 * Greedily remove redundant edges from an existing graph that realizes D between all labeled nodes
	 *
	 * @param graph     an undirected weighted graph
	 * @param isLabeled is node labeled?
	 * @param distance  target symmetric distances over labeled nodes
	 * @return true, if graph changed
	 */
	public static boolean apply(IntGraph graph, IntPredicate isLabeled, BiFunction<Integer, Integer, Integer> distance, ProgressListener progress) throws CanceledException {
		var nodes = graph.getNumberOfNodes();
		var edges = graph.getNumberOfEdges();

		var alive = new HashSet<>(graph.edges());

		// Sort edges by decreasing weight, then by min node, then by max node
		var order = new ArrayList<>(graph.edges());
		order.sort(Comparator
				.comparingInt((IntGraph.Edge e) -> -graph.getWeight(e))
				.thenComparingInt(e -> -Math.min(e.u(), e.v()))
				.thenComparingInt(e -> -Math.max(e.u(), e.v())));

		for (var e : order) {
			if (!alive.contains(e)) continue;

			var aliveWithoutE = new HashSet<>(alive);
			aliveWithoutE.remove(e);
			if (checkAllDistancesBetweenLabeledNodes(graph, isLabeled, aliveWithoutE, distance)) {
				alive.remove(e);
			}
			progress.checkForCancel();
		}
		for (var e : new ArrayList<>(graph.edges())) {
			if (!alive.contains(e)) {
				graph.removeEdge(e);
			}
			progress.checkForCancel();
		}
		var changed = graph.edges().size() < edges;

		if (changed) {
			if (verbose)
				System.err.printf("E removal: n=%4d e=%4d -> n=%4d e=%4d%n", nodes, edges, graph.getNumberOfNodes(), graph.getNumberOfEdges());

			// clean and removeThroughNodes
			{
				var preCleanNodes = graph.getNumberOfNodes();
				var preCleanEdges = graph.getNumberOfEdges();
				GraphCleanAndSmooth.apply(graph, isLabeled, progress);
				if (verbose && (graph.getNumberOfNodes() < preCleanNodes || graph.getNumberOfEdges() < preCleanEdges))
					System.err.printf("Clean:     n=%4d e=%4d -> n=%4d e=%4d%n", preCleanNodes, preCleanEdges, graph.getNumberOfNodes(), graph.getNumberOfEdges());
			}
		}
		return changed;
	}

	// Validate all taxa pair distances against D.
	private static boolean checkAllDistancesBetweenLabeledNodes(IntGraph graph, IntPredicate isLabeled, Set<IntGraph.Edge> alive, BiFunction<Integer, Integer, Integer> distances) {
		var labeledNodes = graph.nodes().stream().filter(isLabeled::test).toList();

		var shortestPaths = DijkstraAllPairs.allShortestPathLengths(graph, labeledNodes, graph.nodes(), alive);

		for (var u : labeledNodes) {
			for (var v : labeledNodes) {
				if (u < v && !distances.apply(u, v).equals(shortestPaths.get(u).get(v)))
					return false;
			}
		}
		return true;
	}
}