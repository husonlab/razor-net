/*
 * Quantization.java Copyright (C) 2025
 *
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
import java.util.function.Function;
import java.util.function.IntPredicate;
import java.util.function.ToDoubleFunction;

/**
 * Quantization of a distance matrix: the distances are scaled to integers by a power of ten ({@link AutoScaler})
 * and doubled, so that all entries are even; violations of the triangle inequality caused by rounding are
 * repaired, and points at distance 0 are merged.
 *
 * @param size            the number of points after merging
 * @param quantMatrix     the quantized matrix
 * @param mapDistanceBack maps a quantized distance back to the input scale
 * @param mapNodeBack     maps a point to the indices of the input taxa that it represents
 * @param isLabeled       is the node id that of a point of the quantized matrix?
 */
public record Quantization(int size, int[][] quantMatrix, ToDoubleFunction<Integer> mapDistanceBack,
						   Function<Integer, Collection<Integer>> mapNodeBack, IntPredicate isLabeled) {
	/**
	 * quantizes the given distance matrix
	 *
	 * @param m      input distance matrix: square, finite, zero diagonal, non-negative entries
	 * @param digits number of significant digits to keep
	 */
	public static Quantization apply(double[][] m, int digits) {
		if (m == null || m.length == 0) throw new IllegalArgumentException("Matrix must be non-null/non-empty.");
		if (digits < 1) throw new IllegalArgumentException("digits must be >=1.");

		final int n = m.length;
		for (int i = 0; i < n; i++) {
			if (m[i] == null || m[i].length != n) {
				throw new IllegalArgumentException("Matrix must be square; row " + i + " has length " +
												   (m[i] == null ? "null" : m[i].length) + ", expected " + n);
			}
		}

		// Validate values; collect max for overflow-safe scaling.
		final double INT_EPS = 1e-12;
		var maxOffDiag = 0.0;
		var minOffDiag = Double.MAX_VALUE;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				var x = m[i][j];
				if (!Double.isFinite(x))
					throw new IllegalArgumentException("Non-finite value at (" + i + "," + j + "): " + x);
				if (i == j) {
					if (Math.abs(x) > INT_EPS)
						throw new IllegalArgumentException("Diagonal (" + i + "," + j + ") must be 0 but is " + x);
				} else {
					if (x < -INT_EPS)
						throw new IllegalArgumentException("Off-diagonal (" + i + "," + j + ") must be >= 0 but is " + x);
					if (x > maxOffDiag) maxOffDiag = x;
					if (x < minOffDiag) minOffDiag = x;
				}
			}
		}

		// Fewer than two taxa means there is no off-diagonal entry at all, so minOffDiag was never assigned;
		// there is nothing to scale and the (empty) encoding below will not consult the scaler.
		if (n < 2)
			minOffDiag = maxOffDiag = 0.0;

		// min == max is allowed: a matrix on two taxa, or one whose distances are all equal
		var scaler = new AutoScaler(minOffDiag, maxOffDiag, digits);

		// scale to integers, then double, so that all entries are even
		var factor = 2;

		var encoded = new int[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (i == j) {
					encoded[i][j] = 0;
				} else {
					var v = scaler.toInt(m[i][j]);
					var even = v * factor; // multiply ALL values by 2
					checkFitsInt(even, i, j);
					encoded[i][j] = (int) even;
				}
			}
		}

		// Repair triangle inequality violations introduced by rounding
		enforceMetric(encoded);

		// Union-Find: merge any indices whose encoded distance is 0
		var uf = new UnionFind(n);
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if (encoded[i][j] == 0) uf.union(i, j);
			}
		}

		// Components in stable order (by smallest original index)
		Map<Integer, List<Integer>> compToMembers = new LinkedHashMap<>();
		for (int i = 0; i < n; i++) {
			var r = uf.find(i);
			compToMembers.computeIfAbsent(r, k -> new ArrayList<>()).add(i);
		}
		var reps = new ArrayList<>(compToMembers.keySet());
		reps.sort(Comparator.comparingInt(a -> compToMembers.get(a).get(0)));

		final int k = reps.size();
		var reduced = new int[k][k];
		var repIndex = new int[k]; // representative original index per component
		for (int a = 0; a < k; a++) {
			var members = compToMembers.get(reps.get(a));
			repIndex[a] = members.get(0);
		}

		for (int a = 0; a < k; a++) {
			for (int b = 0; b < k; b++) {
				if (a == b) {
					reduced[a][b] = 0;
				} else {
					int i = repIndex[a], j = repIndex[b];
					reduced[a][b] = encoded[i][j];
				}
			}
		}

		// mapNodesBack: reduced id -> unmodifiable sorted list of original ids
		Map<Integer, Collection<Integer>> back = new HashMap<>();
		for (int a = 0; a < k; a++) {
			var list = compToMembers.get(reps.get(a)).stream().sorted().map(t -> t).toList();
			back.put(a, list);
		}
		Function<Integer, Collection<Integer>> mapNodesBack = newId -> {
			var r = back.get(newId);
			if (r == null)
				throw new IllegalArgumentException("Unknown reduced id: " + newId);
			return r;
		};

		return new Quantization(reduced.length, reduced, i -> i == 0 ? 0 : scaler.toDouble(i) / factor, mapNodesBack, id -> id < reduced.length);
	}

	/* -------------------- helpers -------------------- */

	private static void checkFitsInt(long val, int i, int j) {
		if (val > Integer.MAX_VALUE || val < Integer.MIN_VALUE) {
			throw new IllegalArgumentException("Encoded value at (" + i + "," + j + ") overflows int: " + val);
		}
	}

	private static void enforceMetric(int[][] d) {
		int n = d.length;
		long totalAdjustment = 0;
		long totalDistance = 0;
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++)
				totalDistance += d[i][j];

		for (int k = 0; k < n; k++) {
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
					long viaK = (long) d[i][k] + d[k][j];
					if (viaK < d[i][j]) {
						totalAdjustment += d[i][j] - viaK;
						d[i][j] = (int) viaK;
					}
				}
			}
		}
		if (totalDistance > 0) {
			var pct = 100.0 * totalAdjustment / totalDistance;
			if (pct > 0.05) {
				System.err.printf("Metric closure: adjusted %d of %d (%.4f%%)%n", totalAdjustment, totalDistance, pct);
			}
		}
	}

	private static final class UnionFind {
		private final int[] parent;
		private final byte[] rank;

		UnionFind(int n) {
			parent = new int[n];
			rank = new byte[n];
			for (int i = 0; i < n; i++) parent[i] = i;
		}

		int find(int x) {
			if (parent[x] != x) parent[x] = find(parent[x]);
			return parent[x];
		}

		void union(int a, int b) {
			int ra = find(a), rb = find(b);
			if (ra == rb) return;
			if (rank[ra] < rank[rb]) parent[ra] = rb;
			else if (rank[ra] > rank[rb]) parent[rb] = ra;
			else {
				parent[rb] = ra;
				rank[ra]++;
			}
		}
	}


	public double[][] createDoubleMatrix() {
		var doubleMatrix = new double[quantMatrix.length][quantMatrix[0].length];
		for (int i = 0; i < quantMatrix.length; i++) {
			for (int j = 0; j < quantMatrix[0].length; j++) {
				doubleMatrix[i][j] = quantMatrix[i][j];
			}
		}
		return doubleMatrix;
	}

}