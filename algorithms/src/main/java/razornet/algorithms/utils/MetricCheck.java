/*
 * MetricCheck.java Copyright (C) 2025
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

import jloda.util.StringUtils;

/**
 * Checks whether a distance matrix is a metric, i.e. satisfies the triangle inequality. A matrix that is not a
 * metric has no realization. The check warns and never refuses.
 * <p>
 * Daniel Huson and Momoko Hayamizu, 8.2026
 */
public final class MetricCheck {
	/** Rough ceiling on triple tests, so the check stays negligible beside the search itself. */
	private static final long WORK_BUDGET = 200_000_000L;

	private MetricCheck() {
	}

	/**
	 * @param exhaustive  true if every triple was tested; when false only {@code sampledK} of the n possible
	 *                    middle points were used, so {@code violations} is a LOWER bound and "none found" is
	 *                    not proof of a metric
	 */
	public record Result(boolean exhaustive, int sampledK, long violations, double worst,
						 int worstI, int worstJ, int worstK) {
		public boolean isMetric() {
			return violations == 0;
		}
	}

	/**
	 * Counts triples with {@code d(i,j) > d(i,k) + d(k,j)}. Cost is |K|*n^2; for large n the middle point k
	 * is sampled on a fixed stride (deterministic, so a run is reproducible) to stay inside the work budget.
	 */
	public static Result check(double[][] d) {
		var n = d.length;
		if (n < 3)
			return new Result(true, n, 0, 0, -1, -1, -1);
		var maxK = (int) Math.max(1, Math.min(n, WORK_BUDGET / ((long) n * n)));
		var stride = Math.max(1, n / maxK);
		var exhaustive = (stride == 1);
		long violations = 0;
		var worst = 0.0;
		int wi = -1, wj = -1, wk = -1, used = 0;
		var tol = 1e-9;
		for (var k = 0; k < n; k += stride) {
			used++;
			for (var i = 0; i < n; i++) {
				var ik = d[i][k];
				for (var j = 0; j < n; j++) {
					var excess = d[i][j] - (ik + d[k][j]);
					if (excess > tol) {
						violations++;
						if (excess > worst) {
							worst = excess;
							wi = i;
							wj = j;
							wk = k;
						}
					}
				}
			}
		}
		return new Result(exhaustive, used, violations, worst, wi, wj, wk);
	}

	/**
	 * Runs the check and prints a warning if the input is not a metric.
	 *
	 * @param context what is being checked, e.g. the algorithm name, so the user can tell which step complained
	 * @return the result, in case the caller wants to record it
	 */
	public static Result checkAndWarn(double[][] d, String context) {
		var r = check(d);
		if (!r.isMetric()) {
			System.err.printf("WARNING (%s): the input is NOT a metric -- %d triangle-inequality violation%s%s, worst by %s.%n",
					context, r.violations(), r.violations() == 1 ? "" : "s",
					r.exhaustive() ? "" : " (sampled: %d of %d middle points, so this is a lower bound)".formatted(r.sampledK(), d.length),
					StringUtils.trim(r.worst()));
			if (r.worstI() >= 0)
				System.err.printf("         e.g. d(%d,%d) = %s > d(%d,%d) + d(%d,%d) = %s (0-based indices).%n",
						r.worstI(), r.worstJ(), StringUtils.trim(d[r.worstI()][r.worstJ()]),
						r.worstI(), r.worstK(), r.worstK(), r.worstJ(),
						StringUtils.trim(d[r.worstI()][r.worstK()] + d[r.worstK()][r.worstJ()]));
			System.err.printf("         No realization exists, so the network cannot reproduce all distances. Consider%n");
			System.err.printf("         preprocessing (e.g. ClassicalMDS2D). Continuing anyway.%n");
		}
		return r;
	}
}
