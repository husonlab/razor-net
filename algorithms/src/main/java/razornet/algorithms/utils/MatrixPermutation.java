/*
 * MatrixPermutation.java Copyright (C) 2026 Daniel H. Huson
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

import java.util.Arrays;
import java.util.Random;

/**
 * A permutation of the indices of a distance matrix, with the translation of node ids back to the original
 * indices. RazorNet uses the canonical permutation, so that the result does not depend on the order of the input.
 */
public final class MatrixPermutation {

	private final int[] forward;   // forward[i] = which original index is at permuted position i

	/**
	 * Random permutation of indices 0 ... n-1.
	 */
	public static MatrixPermutation random(int n, Random random) {
		var perm = new int[n];
		for (int i = 0; i < n; i++) perm[i] = i;
		// Fisher-Yates shuffle
		for (int i = n - 1; i > 0; i--) {
			int j = random.nextInt(i + 1);
			int t = perm[i];
			perm[i] = perm[j];
			perm[j] = t;
		}
		return new MatrixPermutation(perm);
	}

	/**
	 * Canonical permutation: sort indices by row sum (smallest first), with ties
	 * broken by lexicographic order of the sorted distance vector, then by taxon
	 * id. The result is deterministic and independent of any prior ordering
	 * of the input matrix — applying canonical() to two row/column permutations
	 * of the same labeled metric yields the same canonical matrix.
	 *
	 * @param matrix square symmetric integer distance matrix
	 */
	public static MatrixPermutation canonical(int[][] matrix) {
		var n = matrix.length;

		// Precompute sort keys for each original index.
		var rowSums = new long[n];          // primary: smaller first ("more central" taxa)
		var sortedRows = new int[n][];      // secondary: lexicographic on sorted distance vectors
		for (var i = 0; i < n; i++) {
			var sum = 0L;
			for (var j = 0; j < n; j++) sum += matrix[i][j];
			rowSums[i] = sum;
			sortedRows[i] = matrix[i].clone();
			java.util.Arrays.sort(sortedRows[i]);
		}

		var indices = new Integer[n];
		for (var i = 0; i < n; i++) indices[i] = i;
		java.util.Arrays.sort(indices, (a, b) -> {
			int cmp = Long.compare(rowSums[a], rowSums[b]);
			if (cmp != 0) return cmp;
			for (int k = 0; k < n; k++) {
				cmp = Integer.compare(sortedRows[a][k], sortedRows[b][k]);
				if (cmp != 0) return cmp;
			}
			return Integer.compare(a, b);
		});

		return new MatrixPermutation(Arrays.stream(indices).mapToInt(i -> i).toArray());
	}

	/**
	 * Identity permutation (useful for testing).
	 */
	public static MatrixPermutation identity(int[][] matrix) {
		var n = matrix.length;
		var perm = new int[n];
		for (int i = 0; i < n; i++) perm[i] = i;
		return new MatrixPermutation(perm);
	}

	private MatrixPermutation(int[] forward) {
		this.forward = forward;
	}

	/**
	 * Apply permutation to a square matrix: result[i][j] = matrix[forward[i]][forward[j]].
	 */
	public int[][] apply(int[][] matrix) {
		var n = matrix.length;
		if (n != forward.length)
			throw new IllegalArgumentException(
					"Matrix size " + n + " does not match permutation size " + forward.length);
		var result = new int[n][n];
		for (var i = 0; i < n; i++) {
			for (var j = 0; j < n; j++) {
				result[i][j] = matrix[forward[i]][forward[j]];
			}
		}
		return result;
	}

	/**
	 * Translate a graph computed in this permutation's index space back to the
	 * original index space. Original-taxa nodes (IDs &lt; originalSize) are
	 * un-permuted; aux nodes (IDs &ge; originalSize) are passed through unchanged.
	 *
	 * @param graph        graph in permuted-index space
	 * @param originalSize number of original taxa (= permutation length)
	 * @return new graph with translated node IDs; structure and weights preserved
	 */
	public IntGraph unpermuteGraph(IntGraph graph, int originalSize) {
		if (originalSize != forward.length)
			throw new IllegalArgumentException("originalSize " + originalSize + " does not match permutation size " + forward.length);

		var result = new IntGraph();

		// Translate nodes (preserves isolated nodes)
		for (var v : graph.nodes()) {
			result.ensureNode(unpermuteNode(v, originalSize));
		}

		// Translate edges with weights
		for (var e : graph.edges()) {
			var u = unpermuteNode(e.u(), originalSize);
			var v = unpermuteNode(e.v(), originalSize);
			result.addEdge(u, v, graph.getWeight(e));
		}

		return result;
	}

	/**
	 * Translate a node ID from the permuted index space back to the original index space.
	 * Aux nodes (added by the algorithm beyond the original matrix size) are passed through unchanged.
	 *
	 * @param permutedId         node ID in the permuted matrix
	 * @param originalMatrixSize number of original taxa (= permutation length)
	 */
	public int unpermuteNode(int permutedId, int originalMatrixSize) {
		if (permutedId < originalMatrixSize) {
			return forward[permutedId];
		}
		return permutedId;
	}

	/**
	 * Returns the size of the permutation.
	 */
	public int size() {
		return forward.length;
	}
}