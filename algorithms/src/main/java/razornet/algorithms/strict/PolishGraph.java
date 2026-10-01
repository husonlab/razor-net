/*
 * PolishGraph.java Copyright (C) 2026 Daniel H. Huson
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
import jloda.util.progress.ProgressListener;
import razornet.algorithms.utils.MutableD;
import razornet.algorithms.utils.IntGraph;
import razornet.algorithms.utils.GraphFromDistances;
import razornet.algorithms.utils.GreedyGlobalEdgeRemoval;

import java.util.*;

/**
 * Polishing for the Strict strategy: for every vertex of degree at least 3, the metric restricted to the vertex
 * and its neighbours is extended at its slack points (with the slack of this package), the neighbourhood is
 * realized anew and put back into the graph, and the global matrix is extended by the new points. See
 * {@code docs/algorithm.md}.
 * <p>
 * Momoko Hayamizu and Daniel Huson, 2026
 */
public final class PolishGraph {
	/**
	 * One polishing round: visits every vertex once, in decreasing order of the largest slack in its neighbourhood.
	 * Returns true if at least one vertex was polished.
	 */
	public static boolean applyRound(IntGraph graph, MutableD distances, Random random, ProgressListener progress) throws CanceledException {
		var order = new ArrayList<>(graph.nodes());

		// sort the vertices by the largest slack in their neighbourhood, descending
		var slackByVertex = new HashMap<Integer, Integer>();
		for (var v : order) {
			if (graph.getDegree(v) < 3) {
				slackByVertex.put(v, 0);
				continue;
			}
			var neighborhood = new HashSet<>(graph.getAdjacentNodes(v));
			neighborhood.add(v);
			var maxSlack = 0;
			for (var x : neighborhood) {
				var s = DistanceMatrixExtension.slackWithArgmin(distances, neighborhood, x);
				if (s.c() > maxSlack) maxSlack = s.c();
			}
			slackByVertex.put(v, maxSlack);
			progress.checkForCancel();
		}
		order.sort((a, b) -> Integer.compare(slackByVertex.get(b), slackByVertex.get(a)));

		var changed = false;
		for (var v : order) {
			if (graph.getDegree(v) < 3) continue;   // can't polish nodes of degree ≤ 2
			if (polishVertex(graph, distances, v, progress)) {
				changed = true;
			}
		}

		return changed;
	}

	/**
	 * Attempt to polish a single vertex v by locally expanding its neighborhood.
	 * Returns true if the graph and distance matrix were modified.
	 */
	private static boolean polishVertex(IntGraph graph, MutableD distances, int v, ProgressListener progress) throws CanceledException {
		// Step 1: collect the neighborhood and find slack vertices within it
		var neighborhood = collectNeighborhood(graph, v);
		var slackMap = computeSlacks(distances, neighborhood, progress);
		if (slackMap.isEmpty()) return false;

		// Step 2: build a local distance matrix and expand it with new aux nodes
		var ids = new IdMap(neighborhood);
		var localMatrix = buildLocalMatrix(distances, ids);
		var newLocalVertices = expandLocalMatrix(localMatrix, slackMap, ids);
		if (newLocalVertices.isEmpty()) return false;

		// Step 3: realize the local subgraph and graft it back into the global graph
		var globalIds = reserveGlobalIds(distances, newLocalVertices);
		replaceSubgraph(graph, neighborhood, localMatrix, ids, newLocalVertices, globalIds);

		// Step 4: extend the global distance matrix with rows for the new vertices
		fillGlobalDistances(distances, neighborhood, localMatrix, ids, newLocalVertices, globalIds);

		return true;
	}

	// -------------------------------------------------------------------------
	// Step 1: neighborhood and slack computation
	// -------------------------------------------------------------------------

	private static Set<Integer> collectNeighborhood(IntGraph graph, int v) {
		var neighborhood = new HashSet<>(graph.getAdjacentNodes(v));
		neighborhood.add(v);
		return neighborhood;
	}

	private static Map<Integer, DistanceMatrixExtension.Slack> computeSlacks(MutableD distances, Set<Integer> neighborhood, ProgressListener progress) throws CanceledException {
		var slackMap = new HashMap<Integer, DistanceMatrixExtension.Slack>();
		for (var x : neighborhood) {
			var slack = DistanceMatrixExtension.slackWithArgmin(distances, neighborhood, x);
			if (slack.c() > 0) slackMap.put(x, slack);
			progress.checkForCancel();
		}
		return slackMap;
	}

	// -------------------------------------------------------------------------
	// Step 2: local matrix, expansion
	// -------------------------------------------------------------------------

	/**
	 * Bidirectional mapping between global node IDs and local indices in M_v.
	 */
	private static final class IdMap {
		final List<Integer> localToGlobal;       // index = local ID, value = global ID
		final Map<Integer, Integer> globalToLocal;

		IdMap(Collection<Integer> globalIds) {
			localToGlobal = new ArrayList<>(globalIds);
			globalToLocal = new HashMap<>();
			for (int i = 0; i < localToGlobal.size(); i++)
				globalToLocal.put(localToGlobal.get(i), i);
		}
	}

	private static MutableD buildLocalMatrix(MutableD distances, IdMap ids) {
		var n = ids.localToGlobal.size();
		var localMatrix = new MutableD(n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				localMatrix.set(i, j, distances.get(ids.localToGlobal.get(i), ids.localToGlobal.get(j)));
			}
		}
		return localMatrix;
	}

	/**
	 * Adds new auxiliary vertices to the local matrix, one per slack vertex, skipping any that would coincide
	 * with an existing vertex. Returns the local ids of the vertices added.
	 */
	private static List<Integer> expandLocalMatrix(MutableD localMatrix, Map<Integer, DistanceMatrixExtension.Slack> slackMap, IdMap ids) {
		var newLocalVertices = new ArrayList<Integer>();
		var orderedSlackKeys = new ArrayList<>(slackMap.keySet());

		for (var x_global : orderedSlackKeys) {
			var slack_global = slackMap.get(x_global);
			var slack_local = new DistanceMatrixExtension.Slack(
					slack_global.c(),
					ids.globalToLocal.get(slack_global.y()),
					ids.globalToLocal.get(slack_global.z())
			);
			var x_local = ids.globalToLocal.get(x_global);

			var columnLocal = DistanceMatrixExtension.buildAuxColumnUsingMaxAll(localMatrix, x_local, slack_local);
			var coincide = DistanceMatrixExtension.findCoincidingVertex(localMatrix, columnLocal);
			if (coincide == null) {
				newLocalVertices.add(localMatrix.appendVertex(columnLocal));
			}
		}
		return newLocalVertices;
	}

	// -------------------------------------------------------------------------
	// Step 3: subgraph replacement
	// -------------------------------------------------------------------------

	/**
	 * Bidirectional map: local ID of a new vertex ↔ its newly reserved global ID.
	 */
	private record GlobalIds(Map<Integer, Integer> localToGlobal, Map<Integer, Integer> globalToLocal,
							 List<Integer> finalLocalToGlobal) {
	}

	private static GlobalIds reserveGlobalIds(MutableD distances, List<Integer> newLocalVertices) {
		// the rows reserved here are filled in fillGlobalDistances and must not be read before
		var localToGlobal = new HashMap<Integer, Integer>();
		var globalToLocal = new HashMap<Integer, Integer>();
		var finalLocalToGlobal = new ArrayList<Integer>();   // appended later

		for (var xp_local : newLocalVertices) {
			var xp_global = distances.newRowCol();
			localToGlobal.put(xp_local, xp_global);
			globalToLocal.put(xp_global, xp_local);
			finalLocalToGlobal.add(xp_global);
		}
		return new GlobalIds(localToGlobal, globalToLocal, finalLocalToGlobal);
	}

	private static void replaceSubgraph(IntGraph graph, Set<Integer> neighborhood, MutableD localMatrix, IdMap ids,
										List<Integer> newLocalVertices, GlobalIds globalIds) throws CanceledException {
		// Realise the new subgraph from the expanded local matrix.
		var subgraphNodes = new HashSet<Integer>();
		for (int i = 0; i < localMatrix.size(); i++) subgraphNodes.add(i);
		var newSubgraph = GraphFromDistances.realizeFromDistances(localMatrix::get, subgraphNodes, null);

		// Build the full local-to-global mapping (original neighborhood IDs + new IDs).
		var fullLocalToGlobal = new ArrayList<>(ids.localToGlobal);
		fullLocalToGlobal.addAll(globalIds.finalLocalToGlobal);

		// Remove existing edges within the original neighborhood.
		// (Edges from neighborhood nodes to outside-neighborhood nodes survive — that's intentional.)
		for (var u : neighborhood) {
			for (var w : neighborhood) {
				if (u < w) graph.removeEdge(u, w);
			}
		}
		// Add edges from the realised subgraph, translating local → global.
		for (var edge : newSubgraph.edges()) {
			var u_global = fullLocalToGlobal.get(edge.u());
			var v_global = fullLocalToGlobal.get(edge.v());
			graph.addEdge(u_global, v_global, newSubgraph.getWeight(edge));
		}
	}

	// -------------------------------------------------------------------------
	// Step 4: global distance matrix update
	// -------------------------------------------------------------------------

	private static void fillGlobalDistances(MutableD distances, Set<Integer> neighborhood,
											MutableD localMatrix, IdMap ids,
											List<Integer> newLocalVertices, GlobalIds globalIds) {
		for (var xp_local : newLocalVertices) {
			var xp_global = globalIds.localToGlobal.get(xp_local);
			var finalColumn = computeFinalColumn(distances, neighborhood, localMatrix, ids,
					globalIds, xp_local, xp_global);

			// Fill the previously reserved row/col for xp_global. The matrix is symmetric,
			// so writing entries [xp_global, i] for i < xp_global covers both halves.
			for (int i = 0; i < xp_global; i++) {
				distances.set(xp_global, i, finalColumn[i]);
			}
		}
	}

	/**
	 * Compute the distance from a newly added vertex xp to every existing global index.
	 * Uses three cases: distance to xp itself (0), distance to another new or neighborhood
	 * node (read from local matrix), distance to outside-neighborhood node (min over
	 * neighborhood path).
	 */
	private static int[] computeFinalColumn(MutableD distances, Set<Integer> neighborhood,
											MutableD localMatrix, IdMap ids, GlobalIds globalIds,
											int xp_local, int xp_global) {
		var finalColumn = new int[distances.size()];

		for (int i = 0; i < distances.size(); i++) {
			if (i == xp_global) {
				finalColumn[i] = 0;
			} else if (globalIds.globalToLocal.containsKey(i)) {
				finalColumn[i] = localMatrix.get(xp_local, globalIds.globalToLocal.get(i));
			} else if (ids.globalToLocal.containsKey(i)) {
				finalColumn[i] = localMatrix.get(xp_local, ids.globalToLocal.get(i));
			} else {
				finalColumn[i] = minDistanceViaNeighborhood(
						distances, localMatrix, ids, neighborhood, xp_local, i);
			}
		}
		return finalColumn;
	}

	/**
	 * For an outside-neighborhood node i, the distance from xp to i is the min over
	 * neighborhood nodes n of (dist_local(xp, n) + dist_global(n, i)).
	 */
	private static int minDistanceViaNeighborhood(MutableD distances, MutableD localMatrix, IdMap ids,
												  Set<Integer> neighborhood, int xp_local, int i) {
		var best = Integer.MAX_VALUE;
		for (var n_global : neighborhood) {
			var n_local = ids.globalToLocal.get(n_global);
			best = Math.min(best, localMatrix.get(xp_local, n_local) + distances.get(n_global, i));
		}
		return best;
	}
}