/*
 * GraphFromDistances.java Copyright (C) 2026 Daniel H. Huson
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
import jloda.util.SetUtils;
import jloda.util.progress.ProgressListener;

import java.util.*;
import java.util.function.ToIntBiFunction;

/**
 * Realization of a metric as a graph on its points: two points are joined by an edge, weighted by their
 * distance, if no third point lies between them.
 */
public class GraphFromDistances {

	public static IntGraph realizeFromDistances(MutableD distances) {
		try {
			return realizeFromDistances(distances::get, SetUtils.fullIndexSet(distances.size()), null);
		} catch (CanceledException ignored) {
			throw new RuntimeException(); // can't happen
		}
	}

	public static IntGraph realizeFromDistances(int[][] distances) {
		return realizeFromDistances(distances, SetUtils.fullIndexSet(distances.length));
	}

	public static IntGraph realizeFromDistances(int[][] distances, Set<Integer> subset) {
		try {
			return realizeFromDistances((a, b) -> distances[a][b], subset, null);
		} catch (CanceledException ignored) {
			throw new RuntimeException(); // can't happen
		}
	}

	public static IntGraph realizeFromDistances(ToIntBiFunction<Integer, Integer> distance, Set<Integer> subset, ProgressListener progress) throws CanceledException {
		var graph = new IntGraph();
		for (int v : subset) graph.ensureNode(v);

		// add an edge {u,v} if no third point w of the subset has d(u,w) + d(w,v) <= d(u,v)
		for (var u : subset) {
			for (var v : subset) {
				if (u < v && !isRedundant(distance, subset, u, v)) {
					graph.addEdge(u, v, distance.applyAsInt(u, v));
				}
			}
			if (progress != null)
				progress.checkForCancel();
		}
		return graph;
	}

	/**
	 * Edge (u,v) is redundant iff some intermediate w yields d(u,w) + d(w,v) ≤ d(u,v),
	 * with both legs strictly positive (excludes trivial paths through coincident vertices).
	 */
	private static boolean isRedundant(ToIntBiFunction<Integer, Integer> distance, Set<Integer> subset, int u, int v) {
		int duv = distance.applyAsInt(u, v);
		for (int w : subset) {
			if (w == u || w == v) continue;
			int duw = distance.applyAsInt(u, w);
			int dwv = distance.applyAsInt(w, v);
			if (duw == 0 || dwv == 0) continue;
			if (duw + dwv <= duv) return true;
		}
		return false;
	}

	/**
	 * the taxa-only realization of a metric: the graph on the taxa with an edge {x,y} of weight d(x,y) for every
	 * pair such that no third taxon lies on a geodesic between x and y. It realizes the metric exactly, without
	 * auxiliary vertices.
	 */
	public static IntGraph realizeTaxaOnly(int[][] distances, ProgressListener progress) throws CanceledException {
		var n = distances.length;
		var graph = new IntGraph();
		for (var i = 0; i < n; i++)
			graph.ensureNode(i);
		for (var x = 0; x < n; x++) {
			if (progress != null)
				progress.checkForCancel();
			for (var y = x + 1; y < n; y++) {
				// coincident taxa (0) are merged by quantization; guard IntGraph.addEdge
				if (distances[x][y] > 0 && isEssentialEdge(distances, x, y))
					graph.addEdge(x, y, distances[x][y]);
			}
		}
		return graph;
	}

	/**
	 * Does no third point lie on a geodesic between x and y? Points at distance 0 from x or from y are skipped. The
	 * test is >=, so that an edge is also dropped if a detour over a third point is shorter, as can happen for an
	 * input that is not quite a metric.
	 */
	private static boolean isEssentialEdge(int[][] distances, int x, int y) {
		var dxy = distances[x][y];
		for (var z = 0; z < distances.length; z++) {
			if (z == x || z == y)
				continue;
			var dxz = distances[x][z];
			var dzy = distances[z][y];
			if (dxz == 0 || dzy == 0)
				continue;
			if (dxy >= dxz + dzy)
				return false;
		}
		return true;
	}
}
