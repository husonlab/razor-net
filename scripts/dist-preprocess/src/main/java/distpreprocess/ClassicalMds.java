package distpreprocess;

import org.ojalgo.matrix.decomposition.Eigenvalue;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.R064Store;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Classical MDS (Torgerson / Gower principal coordinates).
 * If fewer than {@code dimension} eigenvalues exceed the positivity threshold,
 * this class fails instead of silently reducing the dimension.
 */
public final class ClassicalMds {
	private ClassicalMds() {
	}

	public record EigenReport(
			double[] eigenvaluesDescending,
			int positiveCount,
			int nearZeroCount,
			int negativeCount,
			double tau,
			double sumPositive,
			double sumNegative,
			double negativeMassRatio,
			double[] proportionPositive
	) {
		public void print(PrintStream err) {
			err.printf("eigenvaluesDescending=%s%n", Arrays.toString(eigenvaluesDescending));
			err.printf("positiveCount=%d nearZeroCount=%d negativeCount=%d tau=%.3e%n",
					positiveCount, nearZeroCount, negativeCount, tau);
			err.printf("sumPositive=%.8g sumNegative=%.8g negativeMassRatio=%.8g%n",
					sumPositive, sumNegative, negativeMassRatio);
			for (int i = 0; i < proportionPositive.length; i++) {
				err.printf("proportionPositive[%d]=%.8g%n", i + 1, proportionPositive[i]);
			}
		}
	}

	public record ReconstructionError(
			double RMSE,
			double relativeRMSE,
			double maxAbsError,
			Double maxRelError
	) {
		public void print(PrintStream err) {
			err.printf("RMSE=%.8g relativeRMSE=%.8g maxAbsError=%.8g", RMSE, relativeRMSE, maxAbsError);
			if (maxRelError == null) {
				err.printf(" maxRelError=undefined (no pair with D_in > tau)%n");
			} else {
				err.printf(" maxRelError=%.8g%n", maxRelError);
			}
		}
	}

	public record Result(double[][] distances, EigenReport eigen, ReconstructionError error) {
	}

	public static Result apply(double[][] din, int dimension, PrintStream err) {
		int n = din.length;
		if (dimension < 1) {
			throw new DistPreprocessException("dimension must be >= 1, found " + dimension);
		}
		if (dimension > n - 1) {
			throw new DistPreprocessException(
					"dimension d=" + dimension + " exceeds the classical-MDS maximum n-1=" + (n - 1));
		}

		double[][] a = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				double dij = din[i][j];
				a[i][j] = -0.5 * dij * dij;
			}
		}

		double[][] b = doubleCenter(a);
		EigenDecomp evd = decomposeSymmetric(b);

		double lambdaMax = evd.values[0];
		double tau = Math.max(1e-12, 1e-12 * Math.abs(lambdaMax));

		int positiveCount = 0;
		int nearZeroCount = 0;
		int negativeCount = 0;
		double sumPositive = 0;
		double sumNegative = 0;
		for (double lambda : evd.values) {
			if (lambda > tau) {
				positiveCount++;
				sumPositive += lambda;
			} else if (lambda < -tau) {
				negativeCount++;
				sumNegative += lambda;
			} else {
				nearZeroCount++;
			}
		}

		if (sumPositive == 0) {
			throw new DistPreprocessException(
					"sumPositive is 0; all eigenvalues are <= tau=" + tau + ". eigenvaluesDescending="
							+ Arrays.toString(evd.values));
		}
		if (positiveCount < dimension || evd.values[dimension - 1] <= tau) {
			throw new DistPreprocessException(String.format(
					"Need %d positive eigenvalues (lambda > tau=%.3e) for classical MDS, but positiveCount=%d. eigenvaluesDescending=%s",
					dimension, tau, positiveCount, Arrays.toString(evd.values)));
		}

		double negativeMassRatio = Math.abs(sumNegative) / sumPositive;
		double[] proportionPositive = new double[dimension];
		for (int i = 0; i < dimension; i++) {
			proportionPositive[i] = evd.values[i] / sumPositive;
		}

		EigenReport eigen = new EigenReport(
				evd.values, positiveCount, nearZeroCount, negativeCount, tau,
				sumPositive, sumNegative, negativeMassRatio, proportionPositive);
		eigen.print(err);

		double[][] coords = new double[n][dimension];
		for (int i = 0; i < n; i++) {
			for (int aIdx = 0; aIdx < dimension; aIdx++) {
				coords[i][aIdx] = evd.vectors[i][aIdx] * Math.sqrt(evd.values[aIdx]);
			}
		}

		double[][] dout = new double[n][n];
		for (int i = 0; i < n; i++) {
			dout[i][i] = 0;
			for (int j = i + 1; j < n; j++) {
				double sumSq = 0;
				for (int aIdx = 0; aIdx < dimension; aIdx++) {
					double diff = coords[i][aIdx] - coords[j][aIdx];
					sumSq += diff * diff;
				}
				double dist = Math.sqrt(sumSq);
				dout[i][j] = dist;
				dout[j][i] = dist;
			}
		}

		ReconstructionError error = reconstructionError(din, dout, tau);
		error.print(err);
		return new Result(dout, eigen, error);
	}

	static ReconstructionError reconstructionError(double[][] din, double[][] dmds, double tau) {
		int n = din.length;
		int m = n * (n - 1) / 2;
		if (m == 0) {
			throw new DistPreprocessException("Need at least 2 taxa to compute reconstruction error");
		}
		double sumSqErr = 0;
		double sumSqIn = 0;
		double maxAbs = 0;
		double maxRel = 0;
		boolean anyRel = false;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double diff = dmds[i][j] - din[i][j];
				sumSqErr += diff * diff;
				sumSqIn += din[i][j] * din[i][j];
				double abs = Math.abs(diff);
				if (abs > maxAbs) {
					maxAbs = abs;
				}
				if (din[i][j] > tau) {
					double rel = abs / din[i][j];
					if (!anyRel || rel > maxRel) {
						maxRel = rel;
					}
					anyRel = true;
				}
			}
		}
		double rmse = Math.sqrt(sumSqErr / m);
		if (sumSqIn == 0) {
			throw new DistPreprocessException("relativeRMSE is undefined: all input pairwise distances are 0");
		}
		double relativeRmse = Math.sqrt(sumSqErr / sumSqIn);
		return new ReconstructionError(rmse, relativeRmse, maxAbs, anyRel ? maxRel : null);
	}

	/** B = H A H with H = I - (1/n) 11^T. */
	static double[][] doubleCenter(double[][] a) {
		int n = a.length;
		double invN = 1.0 / n;
		double[] rowMean = new double[n];
		double[] colMean = new double[n];
		double grand = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				rowMean[i] += a[i][j];
				colMean[j] += a[i][j];
				grand += a[i][j];
			}
		}
		for (int i = 0; i < n; i++) {
			rowMean[i] *= invN;
			colMean[i] *= invN;
		}
		grand *= invN * invN;

		double[][] b = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				b[i][j] = a[i][j] - rowMean[i] - colMean[j] + grand;
			}
		}
		return b;
	}

	private record EigenDecomp(double[] values, double[][] vectors) {
	}

	private static EigenDecomp decomposeSymmetric(double[][] b) {
		int n = b.length;
		R064Store store = R064Store.FACTORY.make(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				store.set(i, j, b[i][j]);
			}
		}

		Eigenvalue<Double> evd = Eigenvalue.R064.make(store, true);
		if (!evd.decompose(store)) {
			throw new DistPreprocessException("Symmetric eigenvalue decomposition failed");
		}

		double[] rawValues = new double[n];
		for (int i = 0; i < n; i++) {
			rawValues[i] = evd.getEigenvalues().doubleValue(i);
		}
		MatrixStore<Double> v = evd.getV();
		double[][] rawVectors = new double[n][n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				rawVectors[i][j] = v.doubleValue(i, j);
			}
		}

		Integer[] order = new Integer[n];
		for (int i = 0; i < n; i++) {
			order[i] = i;
		}
		Arrays.sort(order, (p, q) -> Double.compare(rawValues[q], rawValues[p]));

		double[] values = new double[n];
		double[][] vectors = new double[n][n];
		for (int dest = 0; dest < n; dest++) {
			int src = order[dest];
			values[dest] = rawValues[src];
			for (int i = 0; i < n; i++) {
				vectors[i][dest] = rawVectors[i][src];
			}
		}
		return new EigenDecomp(values, vectors);
	}
}
