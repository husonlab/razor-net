/*
 * TightSpanConstruction.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.tightspan;

import jloda.util.CanceledException;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.utils.IntGraph;
import razornet.algorithms.utils.Quantization;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;

/**
 * The Tight Span Walker construction: realizes a metric exactly by walking, for every pair of taxa, a geodesic
 * inside the tight span of the metric; the network is the union of these walks. A point of the tight span is
 * given by its distances to the taxa ({@link TightSpanPoint}). In each step, the Kuratowski graph of the current
 * point f (taxa i,j are adjacent if f[i] + f[j] = d[i][j]) determines the direction, and the step is the largest
 * that keeps the new point in the tight span and reduces the remaining distance by exactly its length. Where no
 * positive step is admissible, or the cap on the number of points is reached, the two points are joined by a
 * direct edge, so that the network always realizes the metric. All arithmetic is in integers.
 * <p>
 * Based on the implementation in PopART (Leigh &amp; Bryant 2015, {@code src/networks/TightSpanWalker.cpp},
 * GPL-3), modified as described above. The tight span is due to Isbell (1964) and Dress (1984).
 *
 * @author Daniel Huson, 8.2026
 */
public class TightSpanConstruction {
	// the three states of a coordinate under a step; the values are the multipliers applied to delta, so that
	// newCoordinate = oldCoordinate + state * delta, and 'green decreases, red increases' is literally the arithmetic
	private static final int GREEN = -1;
	private static final int FROZEN = 0;
	private static final int RED = 1;

	/**
	 * the result of a run
	 *
	 * @param graph          the network, taxa at node ids 0..n-1, weights in quantized units
	 * @param numberOfPoints the number of distinct tight-span points created (taxa included)
	 * @param directEdges    the number of walks that ended in a direct edge
	 * @param hitPointCap    true if the cap on the number of points was reached
	 */
	public record Result(IntGraph graph, int numberOfPoints, int directEdges, boolean hitPointCap) {
	}

	/**
	 * run the construction
	 *
	 * @param distances the quantized distance matrix (every entry an even integer); not modified
	 * @param maxPoints cap on the number of distinct points; once it is reached, the remaining walks end in direct
	 *                  edges
	 */
	public static Result apply(int[][] distances, int maxPoints, ProgressListener progress) throws CanceledException {
		var ntax = distances.length;
		var graph = new IntGraph();

		// the taxa have the node ids 0..ntax-1; the point of a taxon is its row of the matrix
		var points = new ArrayList<TightSpanPoint>();
		var pointToNode = new HashMap<TightSpanPoint, Integer>();
		for (var i = 0; i < ntax; i++) {
			graph.ensureNode(i);
			var point = new TightSpanPoint(distances[i].clone());
			points.add(point);
			// the rows are distinct: quantization has merged taxa at distance 0
			pointToNode.putIfAbsent(point, i);
		}

		var directEdges = 0;
		var hitPointCap = false;

		progress.setMaximum((long) ntax * (ntax - 1) / 2);
		progress.setProgress(0);

		for (var i = 0; i < ntax; i++) {
			for (var j = 0; j < i; j++) {
				if (walk(i, j, distances, graph, points, pointToNode, maxPoints, progress))
					directEdges++;
				if (points.size() >= maxPoints)
					hitPointCap = true;
				progress.incrementProgress();
			}
		}
		return new Result(graph, points.size(), directEdges, hitPointCap);
	}

	/**
	 * walks a geodesic from the point at node {@code fromNode} to the point at node {@code toNode}, adding the
	 * points and edges that it passes through to the graph
	 *
	 * @return true if the walk ended in a direct edge (no positive step available, or the point cap was reached)
	 */
	private static boolean walk(int fromNode, int toNode, int[][] distances, IntGraph graph,
								List<TightSpanPoint> points, HashMap<TightSpanPoint, Integer> pointToNode,
								int maxPoints, ProgressListener progress) throws CanceledException {
		var currentNode = fromNode;
		var target = points.get(toNode);

		while (true) {
			progress.checkForCancel();

			var current = points.get(currentNode);
			var remaining = TightSpanPoint.distance(current, target);
			if (remaining == 0)
				return false;  // already there (only reachable if two walks met at the target)

			// the cap on the number of points is reached: join directly
			if (points.size() >= maxPoints) {
				addEdge(graph, currentNode, toNode, remaining);
				return true;
			}

			var kuratowski = kuratowskiGraph(current, distances);
			var direction = direction(current, target, kuratowski);
			var delta = stepLength(current, target, distances, direction, remaining);

			if (delta <= 0) {
				// no admissible positive step: join directly
				addEdge(graph, currentNode, toNode, remaining);
				return true;
			}
			if (delta >= remaining) {
				addEdge(graph, currentNode, toNode, remaining);
				return false;   // arrived
			}

			var next = step(current, direction, delta);
			var nextNode = pointToNode.get(next);
			if (nextNode == null) {
				nextNode = points.size();
				points.add(next);
				pointToNode.put(next, nextNode);
				graph.ensureNode(nextNode);
			}
			addEdge(graph, currentNode, nextNode, (int) delta);
			currentNode = nextNode;
		}
	}

	/**
	 * the Kuratowski graph at a point f: taxa i,j are adjacent if the pair is tight there, i.e. f[i] + f[j] = d[i][j]
	 */
	private static BitSet[] kuratowskiGraph(TightSpanPoint point, int[][] distances) {
		var ntax = point.size();
		var adjacency = new BitSet[ntax];
		for (var i = 0; i < ntax; i++)
			adjacency[i] = new BitSet(ntax);
		for (var i = 0; i < ntax; i++) {
			for (var j = i + 1; j < ntax; j++) {
				if (point.slack(i, j, distances) == 0) {
					adjacency[i].set(j);
					adjacency[j].set(i);
				}
			}
		}
		return adjacency;
	}

	/**
	 * the direction to move in: green (-1) for the coordinates that decrease, red (+1) for those that increase,
	 * frozen (0) for the rest. A coordinate starts with the sign of f[i] - g[i], where g is the target; a red
	 * coordinate without a green neighbour in the Kuratowski graph becomes frozen, and a frozen coordinate all of
	 * whose neighbours are red becomes green.
	 */
	private static int[] direction(TightSpanPoint point, TightSpanPoint target, BitSet[] kuratowski) {
		var ntax = point.size();
		var direction = new int[ntax];
		for (var i = 0; i < ntax; i++) {
			var difference = point.get(i) - target.get(i);
			direction[i] = (difference > 0 ? GREEN : difference < 0 ? RED : FROZEN);
		}

		// a red coordinate needs a green neighbour; otherwise it is frozen
		for (var i = 0; i < ntax; i++) {
			if (direction[i] == RED && !hasNeighbourOfColour(kuratowski[i], direction, GREEN))
				direction[i] = FROZEN;
		}

		// a frozen coordinate all of whose neighbours are red becomes green; repeat until nothing changes
		for (var changed = true; changed; ) {
			changed = false;
			for (var i = 0; i < ntax; i++) {
				if (direction[i] == FROZEN && !hasNeighbourOfColour(kuratowski[i], direction, FROZEN)
					&& !hasNeighbourOfColour(kuratowski[i], direction, GREEN)) {
					direction[i] = GREEN;
					changed = true;
				}
			}
		}
		return direction;
	}

	private static boolean hasNeighbourOfColour(BitSet neighbours, int[] direction, int colour) {
		for (var j = neighbours.nextSetBit(0); j >= 0; j = neighbours.nextSetBit(j + 1)) {
			if (direction[j] == colour)
				return true;
		}
		return false;
	}

	/**
	 * the largest admissible step: no pair may get negative slack and no coordinate may become negative (C1, C2),
	 * and the remaining distance to the target has to decrease by exactly the step (C3-C6)
	 *
	 * @return the step length, in quantized units; {@code <= 0} if there is no admissible step
	 */
	private static long stepLength(TightSpanPoint point, TightSpanPoint target, int[][] distances,
								   int[] direction, int remaining) {
		var ntax = point.size();
		var delta = (long) remaining;    // C5: never step past the target

		for (var i = 0; i < ntax; i++) {
			var difference = (long) point.get(i) - target.get(i);
			switch (direction[i]) {
				case GREEN -> {
					delta = Math.min(delta, (remaining + difference) / 2);      // C3 (difference >= 0 here)
					for (var j = 0; j < ntax; j++) {
						if (direction[j] == GREEN) {
							if (j <= i)                                         // each pair once; j == i gives f[i]
								delta = Math.min(delta, point.slack(i, j, distances) / 2);   // C1
						} else if (direction[j] == FROZEN)
							delta = Math.min(delta, point.slack(i, j, distances));           // C2
					}
				}
				case RED -> delta = Math.min(delta, (remaining - difference) / 2);           // C4 (difference < 0)
				case FROZEN -> delta = Math.min(delta, remaining - Math.abs(difference));    // C6
				default -> throw new IllegalStateException("Unexpected direction: " + direction[i]);
			}
		}
		return delta;
	}

	/**
	 * the new point reached by moving {@code delta} in the given direction
	 */
	private static TightSpanPoint step(TightSpanPoint point, int[] direction, long delta) {
		var coordinates = new int[point.size()];
		for (var i = 0; i < coordinates.length; i++)
			coordinates[i] = (int) (point.get(i) + direction[i] * delta);
		return new TightSpanPoint(coordinates);
	}

	/**
	 * add an edge, unless the two nodes are already joined (walks share structure, so this happens constantly) or
	 * the weight is not positive
	 */
	private static void addEdge(IntGraph graph, int u, int v, int weight) {
		if (u != v && weight > 0 && !graph.hasEdge(u, v))
			graph.addEdge(u, v, weight);
	}
}
