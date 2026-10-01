/*
 * DijkstraAllPairs.java Copyright (C) 2026 Daniel H. Huson
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

/**
 * Shortest path distances using Dijkstra's algorithm
 * Daniel Huson, 9.2025
 */
public final class DijkstraAllPairs {

	public static Map<Integer, Map<Integer, Integer>> allShortestPathLengths(IntGraph graph) {
		return allShortestPathLengths(graph, graph.nodes(), graph.nodes(), null);
	}

	/**
	 * the shortest-path distances between all pairs of the given nodes, using only the allowed nodes and edges; the
	 * graph is undirected
	 *
	 * @param graph        the graph
	 * @param nodes        the nodes between which the distances are wanted
	 * @param allowedNodes the nodes that may be used; must contain {@code nodes}
	 * @param allowedEdges the edges that may be used; null for all edges
	 * @return source -> (target -> distance); Integer.MAX_VALUE if the target cannot be reached
	 */
	public static Map<Integer, Map<Integer, Integer>> allShortestPathLengths(IntGraph graph, Collection<Integer> nodes,
																			 Set<Integer> allowedNodes, Set<IntGraph.Edge> allowedEdges) {
		for (var u : nodes)
			if (!allowedNodes.contains(u))
				throw new IllegalArgumentException("node not an allowed node");

		var result = new HashMap<Integer, Map<Integer, Integer>>();
		for (var s : nodes) {
			result.put(s, shortestPathLengthsFromSource(graph, s, nodes, allowedNodes, allowedEdges));
		}
		return result;
	}

	/**
	 * shortest path lengths from source to all target nodes, only using the allowed nodes and edges
	 *
	 * @param graph        the graph
	 * @param source       the source
	 * @param targetNodes  the targets
	 * @param allowedNodes the allowed nodes, must include source and targets
	 * @param allowedEdges the allowed edges
	 * @return mapping of target nodes to distances from source
	 */
	public static Map<Integer, Integer> shortestPathLengthsFromSource(IntGraph graph, Integer source, Collection<Integer> targetNodes,
																	  Set<Integer> allowedNodes, Set<IntGraph.Edge> allowedEdges) {
		// Distance & predecessor edge maps
		var dist = new HashMap<Integer, Integer>();

		// Initialize distances
		for (var n : allowedNodes) {
			dist.put(n, Integer.MAX_VALUE);
		}
		dist.put(source, 0);

		record PQEntry<Node>(Node node, double dist) {
		}

		// Min-heap by tentative distance
		var pq = new PriorityQueue<PQEntry<Integer>>(Comparator.comparingDouble(e -> e.dist));
		pq.add(new PQEntry<>(source, 0.0));

		// Standard Dijkstra with filtering
		while (!pq.isEmpty()) {
			var cur = pq.poll();
			var u = cur.node;

			if (cur.dist > dist.get(u)) continue; // stale entry

			// Iterate over *all* adjacent edges, then filter to allowed
			for (var e : graph.getAdjacentEdges(u)) {
				if (allowedEdges != null && !allowedEdges.contains(e)) continue;

				var v = e.other(u);  // <-- IMPORTANT: undirected neighbor
				if (!allowedNodes.contains(v)) continue;

				var w = graph.getWeight(e);
				if (w < 0) {
					throw new IllegalArgumentException("Dijkstra requires nonnegative weights: " + w);
				}

				var alt = dist.get(u) + w;
				if (alt < dist.get(v)) {
					dist.put(v, alt);
					pq.add(new PQEntry<>(v, alt));
				}
			}
		}

		// Build per-target PathResult
		var out = new HashMap<Integer, Integer>();
		for (var t : targetNodes) {
			var d = dist.getOrDefault(t, Integer.MAX_VALUE);
			if (d == Integer.MAX_VALUE) {
				out.put(t, Integer.MAX_VALUE);
			} else if (t.equals(source)) {
				out.put(t, 0);
			} else {
				out.put(t, d);
			}
		}
		return out;
	}
}