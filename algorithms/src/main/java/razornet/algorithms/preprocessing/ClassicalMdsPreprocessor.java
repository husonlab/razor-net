/*
 * ClassicalMdsPreprocessor.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.preprocessing;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;

/**
 * Classical multidimensional scaling (Torgerson 1952; Gower 1966): embeds the taxa in Euclidean space of the
 * given dimension, using the largest eigenvalues of the doubly centered matrix of squared distances, and returns
 * the Euclidean distances of the embedding. If the matrix does not have enough positive eigenvalues, the input
 * distances are returned unchanged. The reconstruction error is reported on stderr.
 */
public class ClassicalMdsPreprocessor implements DistancesPreprocessor {
	private static final double RELATIVE_EIGEN_THRESHOLD = 1e-12;

	private final int dimension;

	/**
	 * @param dimension the target embedding dimension (typically 2 for a planar network, or 3)
	 */
	public ClassicalMdsPreprocessor(int dimension) {
		if (dimension < 1)
			throw new IllegalArgumentException("dimension must be >= 1, got " + dimension);
		this.dimension = dimension;
	}

	@Override
	public double[][] apply(double[][] distances) {
		var n = distances.length;
		if (dimension > n - 1) {
			System.err.printf("classical MDS: dimension %d exceeds the maximum n-1=%d; keeping the input distances%n", dimension, n - 1);
			return distances;
		}

		// B = double-centered (-1/2 d^2)
		var a = new double[n][n];
		for (var i = 0; i < n; i++) {
			for (var j = 0; j < n; j++) {
				var d = distances[i][j];
				a[i][j] = -0.5 * d * d;
			}
		}
		var b = doubleCenter(a);

		// symmetric eigendecomposition; Jama returns eigenvalues in ascending order, so sort indices descending
		EigenvalueDecomposition evd = new Matrix(b).eig();
		var values = evd.getRealEigenvalues();
		var vectors = evd.getV().getArray(); // vectors[i][k] = i-th component of the k-th eigenvector
		var order = new Integer[n];
		for (var i = 0; i < n; i++)
			order[i] = i;
		java.util.Arrays.sort(order, (p, q) -> Double.compare(values[q], values[p]));

		var lambdaMax = values[order[0]];
		var tau = Math.max(RELATIVE_EIGEN_THRESHOLD, RELATIVE_EIGEN_THRESHOLD * Math.abs(lambdaMax));

		var positiveCount = 0;
		var sumPositive = 0.0;
		for (var lambda : values) {
			if (lambda > tau) {
				positiveCount++;
				sumPositive += lambda;
			}
		}
		if (positiveCount < dimension || values[order[dimension - 1]] <= tau) {
			System.err.printf("classical MDS: only %d positive eigenvalue(s), need %d for a %d-D embedding; keeping the input distances%n",
					positiveCount, dimension, dimension);
			return distances;
		}

		// coordinates: coords[i][k] = eigenvector_k[i] * sqrt(lambda_k) for the top `dimension` eigenvalues
		var coords = new double[n][dimension];
		var capturedMass = 0.0;
		for (var k = 0; k < dimension; k++) {
			var idx = order[k];
			var scale = Math.sqrt(values[idx]);
			capturedMass += values[idx];
			for (var i = 0; i < n; i++)
				coords[i][k] = vectors[i][idx] * scale;
		}

		var out = new double[n][n];
		for (var i = 0; i < n; i++) {
			for (var j = i + 1; j < n; j++) {
				var sumSq = 0.0;
				for (var k = 0; k < dimension; k++) {
					var diff = coords[i][k] - coords[j][k];
					sumSq += diff * diff;
				}
				var dist = Math.sqrt(sumSq);
				out[i][j] = dist;
				out[j][i] = dist;
			}
		}

		// the embedded distances are not rounded here; the quantization (significant digits) does that
		reportError(distances, out, tau, positiveCount, n, sumPositive > 0 ? capturedMass / sumPositive : 0);
		return out;
	}

	/**
	 * writes the embedding quality (variance captured, and how far the embedded distances are from the input) to stderr
	 */
	private void reportError(double[][] din, double[][] dout, double tau, int positiveCount, int n, double capturedFraction) {
		var m = n * (n - 1) / 2;
		var sumSqErr = 0.0;
		var sumSqIn = 0.0;
		var maxRel = 0.0;
		var anyRel = false;
		for (var i = 0; i < n; i++) {
			for (var j = i + 1; j < n; j++) {
				var diff = dout[i][j] - din[i][j];
				sumSqErr += diff * diff;
				sumSqIn += din[i][j] * din[i][j];
				if (din[i][j] > tau) {
					var rel = Math.abs(diff) / din[i][j];
					if (!anyRel || rel > maxRel)
						maxRel = rel;
					anyRel = true;
				}
			}
		}
		var relativeRmse = (sumSqIn > 0) ? Math.sqrt(sumSqErr / sumSqIn) : 0;
		System.err.printf("classical MDS to %d-D: %d positive eigenvalue(s), top %d capture %.1f%% of the positive eigen-mass%n",
				dimension, positiveCount, dimension, 100 * capturedFraction);
		System.err.printf("classical MDS reconstruction: relativeRMSE=%.4g%s (0 = the input was already %d-D Euclidean)%n",
				relativeRmse, anyRel ? String.format(", maxRelError=%.4g", maxRel) : "", dimension);
	}

	/** B = H A H with the centering matrix H = I - (1/n) 11^T. */
	static double[][] doubleCenter(double[][] a) {
		var n = a.length;
		var invN = 1.0 / n;
		var rowMean = new double[n];
		var colMean = new double[n];
		var grand = 0.0;
		for (var i = 0; i < n; i++) {
			for (var j = 0; j < n; j++) {
				rowMean[i] += a[i][j];
				colMean[j] += a[i][j];
				grand += a[i][j];
			}
		}
		for (var i = 0; i < n; i++) {
			rowMean[i] *= invN;
			colMean[i] *= invN;
		}
		grand *= invN * invN;

		var b = new double[n][n];
		for (var i = 0; i < n; i++) {
			for (var j = 0; j < n; j++) {
				b[i][j] = a[i][j] - rowMean[i] - colMean[j] + grand;
			}
		}
		return b;
	}
}
