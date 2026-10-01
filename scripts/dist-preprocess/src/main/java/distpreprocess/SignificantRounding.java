package distpreprocess;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * Rounds distances to a given number of significant digits (HALF_UP).
 * Only the upper triangle is rounded; the lower triangle is copied and the
 * diagonal is left at 0 so symmetry is exact.
 */
public final class SignificantRounding {
	private SignificantRounding() {
	}

	public static BigDecimal roundValue(double value, int significantDigits) {
		if (significantDigits < 1) {
			throw new DistPreprocessException("significantDigits must be >= 1, found " + significantDigits);
		}
		if (value == 0.0) {
			return BigDecimal.ZERO;
		}
		BigDecimal bd = new BigDecimal(Double.toString(value));
		MathContext mc = new MathContext(significantDigits, RoundingMode.HALF_UP);
		return bd.round(mc);
	}

	public static String formatRounded(double value, int significantDigits) {
		return roundValue(value, significantDigits).stripTrailingZeros().toPlainString();
	}

	/**
	 * Returns a new symmetric matrix with rounded off-diagonal entries.
	 */
	public static double[][] apply(double[][] input, int significantDigits) {
		int n = input.length;
		double[][] out = new double[n][n];
		for (int i = 0; i < n; i++) {
			out[i][i] = 0;
			for (int j = i + 1; j < n; j++) {
				double rounded = roundValue(input[i][j], significantDigits).doubleValue();
				out[i][j] = rounded;
				out[j][i] = rounded;
			}
		}
		return out;
	}
}
