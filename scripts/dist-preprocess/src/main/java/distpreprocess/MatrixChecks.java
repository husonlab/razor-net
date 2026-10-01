package distpreprocess;

import java.io.PrintStream;
import java.util.List;

/**
 * Diagonal / symmetry / tiny-negative checks, and triangle-inequality inspection.
 * Values outside the numerical tolerance are errors. Triangle violations are
 * reported and never repaired.
 */
public final class MatrixChecks {
	public static final double EPSILON = 1e-12;

	private MatrixChecks() {
	}

	public record SanitizeReport(
			int diagonalFixed,
			double maxDiagonalAbs,
			int pairsSymmetrized,
			double maxAsymmetry,
			int tinyNegativesFixed,
			double maxTinyNegativeAbs
	) {
		public void print(PrintStream err) {
			err.printf("sanitize: diagonalFixed=%d maxDiagonalAbs=%.3e pairsSymmetrized=%d maxAsymmetry=%.3e tinyNegativesFixed=%d maxTinyNegativeAbs=%.3e%n",
					diagonalFixed, maxDiagonalAbs, pairsSymmetrized, maxAsymmetry, tinyNegativesFixed, maxTinyNegativeAbs);
		}
	}

	public record TriangleReport(
			int triangleInequalityViolations,
			double maxTriangleSlack,
			int worstI,
			int worstJ,
			int worstK,
			double worstDij,
			double worstDik,
			double worstDkj
	) {
		public void print(PrintStream err, String stage, List<String> labels) {
			err.printf("%s: triangleInequalityViolations=%d maxTriangleSlack=%.8g%n",
					stage, triangleInequalityViolations, maxTriangleSlack);
			if (triangleInequalityViolations > 0) {
				String li = labels.get(worstI);
				String lj = labels.get(worstJ);
				String lk = labels.get(worstK);
				err.printf("%s: worstTriple labels=(%s, %s, %s) D(%s,%s)=%.8g D(%s,%s)=%.8g D(%s,%s)=%.8g%n",
						stage, li, lj, lk,
						li, lj, worstDij,
						li, lk, worstDik,
						lk, lj, worstDkj);
			}
		}
	}

	/**
	 * Copy {@code input}, fix only entries within {@link #EPSILON}, and fail on larger defects.
	 */
	public static double[][] sanitize(double[][] input, List<String> labels, PrintStream err) {
		int n = input.length;
		double[][] d = new double[n][n];
		for (int i = 0; i < n; i++) {
			System.arraycopy(input[i], 0, d[i], 0, n);
		}

		int diagonalFixed = 0;
		double maxDiagonalAbs = 0;
		for (int i = 0; i < n; i++) {
			double abs = Math.abs(d[i][i]);
			if (abs > maxDiagonalAbs) {
				maxDiagonalAbs = abs;
			}
			if (abs > EPSILON) {
				throw new DistPreprocessException(String.format(
						"Diagonal entry D[%s,%s]=%.8g exceeds tolerance %g and will not be treated as zero",
						labels.get(i), labels.get(i), d[i][i], EPSILON));
			}
			if (abs > 0) {
				d[i][i] = 0;
				diagonalFixed++;
			}
		}

		int pairsSymmetrized = 0;
		double maxAsymmetry = 0;
		int tinyNegativesFixed = 0;
		double maxTinyNegativeAbs = 0;

		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double a = d[i][j];
				double b = d[j][i];
				double asym = Math.abs(a - b);
				if (asym > maxAsymmetry) {
					maxAsymmetry = asym;
				}
				if (asym > EPSILON) {
					throw new DistPreprocessException(String.format(
							"Asymmetry |D[%s,%s]-D[%s,%s]|=%.8g exceeds tolerance %g",
							labels.get(i), labels.get(j), labels.get(j), labels.get(i), asym, EPSILON));
				}
				if (asym > 0) {
					double mean = 0.5 * (a + b);
					d[i][j] = mean;
					d[j][i] = mean;
					pairsSymmetrized++;
					a = mean;
				} else {
					d[j][i] = a;
				}

				if (a < -EPSILON) {
					throw new DistPreprocessException(String.format(
							"Negative distance D[%s,%s]=%.8g is below -%g",
							labels.get(i), labels.get(j), a, EPSILON));
				}
				if (a < 0) {
					double abs = Math.abs(a);
					if (abs > maxTinyNegativeAbs) {
						maxTinyNegativeAbs = abs;
					}
					d[i][j] = 0;
					d[j][i] = 0;
					tinyNegativesFixed++;
				}
			}
		}

		SanitizeReport report = new SanitizeReport(
				diagonalFixed, maxDiagonalAbs, pairsSymmetrized, maxAsymmetry, tinyNegativesFixed, maxTinyNegativeAbs);
		report.print(err);
		return d;
	}

	/**
	 * Count triangle-inequality violations. Does not modify {@code d}.
	 */
	public static TriangleReport checkTriangleInequality(double[][] d) {
		int n = d.length;
		int violations = 0;
		double maxSlack = 0;
		int wi = -1, wj = -1, wk = -1;
		double wdij = 0, wdik = 0, wdkj = 0;

		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				double dij = d[i][j];
				for (int k = 0; k < n; k++) {
					if (k == i || k == j) {
						continue;
					}
					double dik = d[i][k];
					double dkj = d[k][j];
					double slack = dij - (dik + dkj);
					if (dij > dik + dkj + EPSILON) {
						violations++;
						if (slack > maxSlack) {
							maxSlack = slack;
							wi = i;
							wj = j;
							wk = k;
							wdij = dij;
							wdik = dik;
							wdkj = dkj;
						}
					}
				}
			}
		}
		return new TriangleReport(violations, maxSlack, wi, wj, wk, wdij, wdik, wdkj);
	}
}
