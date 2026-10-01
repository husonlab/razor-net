/*
 * SplitDecompositionPreprocessor.java Copyright (C) 2026 Daniel H. Huson
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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Split pruning: computes the split decomposition of the input metric (Bandelt &amp; Dress 1992),
 * {@code d = sum_S alpha_S * delta_S + residue}, and subtracts the contributions of the splits of smallest weight,
 * up to a given fraction of the total split weight. The result is again a metric, with {@code 0 <= d' <= d}. The
 * split decomposition is a port of the one in SplitsTree6.
 */
public class SplitDecompositionPreprocessor implements DistancesPreprocessor {
	private static final int MAX_TAXA = 500;      // split decomposition is superpolynomial; above this, skip
	private static final double EPSILON = 1e-7;

	private final double removeWeightFraction;    // remove the smallest splits that together make up this fraction of the total split weight

	/**
	 * @param removeWeightFraction remove the smallest-weight splits whose weights sum to at most this fraction of
	 *                             the total split weight (e.g. 0.01 = the bottom 1%)
	 */
	public SplitDecompositionPreprocessor(double removeWeightFraction) {
		this.removeWeightFraction = removeWeightFraction;
	}

	@Override
	public double[][] apply(double[][] distances) {
		var n = distances.length;
		if (n < 4 || removeWeightFraction <= 0) // <4 taxa: only trivial splits, nothing to prune
			return distances;
		if (n > MAX_TAXA) {
			System.err.printf("split-decomposition preprocessing skipped: %d taxa exceeds the limit of %d%n", n, MAX_TAXA);
			return distances;
		}

		var splits = decompose(distances);
		if (splits.isEmpty())
			return distances;
		var total = 0.0;
		for (var s : splits)
			total += s.weight;
		if (total <= 0)
			return distances;

		// smallest first, remove while the cumulative removed weight stays within the requested fraction of the total
		splits.sort((a, b) -> Double.compare(a.weight, b.weight));
		var out = new double[n][n];
		for (var i = 0; i < n; i++)
			out[i] = distances[i].clone();

		var budget = removeWeightFraction * total;
		var removedWeight = 0.0;
		var removedCount = 0;
		for (var s : splits) {
			if (removedWeight + s.weight > budget)
				break;
			removedWeight += s.weight;
			removedCount++;
			// subtract this split's weight from every pair it separates (side A vs side B); taxa are 1-based in the BitSet
			for (var i = s.a.nextSetBit(1); i >= 0; i = s.a.nextSetBit(i + 1))
				for (var j = 1; j <= n; j++)
					if (!s.a.get(j)) {
						out[i - 1][j - 1] -= s.weight;
						out[j - 1][i - 1] -= s.weight;
					}
		}

		System.err.printf("split-decomposition preprocessing: %d splits (total weight %.4g); removed the %d smallest "
						+ "(%.4g, %.2f%% of total) to keep the strong splits and the residue%n",
				splits.size(), total, removedCount, removedWeight, total > 0 ? 100 * removedWeight / total : 0);
		return out;
	}

	private record Split(BitSet a, double weight) {
	}

	/**
	 * canonical split decomposition (Bandelt &amp; Dress 1992), ported from splitstree6; returns the positive-weight
	 * weakly-compatible splits (side A is a 1-based taxon set). Weight-0 trivial splits are omitted (they do not
	 * affect the reconstruction).
	 */
	private static List<Split> decompose(double[][] d) {
		var n = d.length;
		List<Split> previousSplits = new ArrayList<>();
		var previousTaxa = new BitSet();
		previousTaxa.set(1);

		for (var t = 2; t <= n; t++) {
			List<Split> nextSplits = new ArrayList<>();

			// does {t} vs the previous taxa form a split?
			{
				var at = new BitSet();
				at.set(t);
				var wgt = isolationIndex(t, at, previousTaxa, d);
				if (wgt > 0)
					nextSplits.add(new Split(at, wgt));
			}

			// extend each previously found split by placing t on either side
			for (var prev : previousSplits) {
				var a = prev.a;                       // side A within taxa 1..t-1 (does not contain t)
				var b = complement(a, t - 1);         // side B within taxa 1..t-1

				var aWithT = (BitSet) a.clone();
				aWithT.set(t);
				var wgt1 = Math.min(prev.weight, isolationIndex(t, aWithT, b, d));
				if (wgt1 > 0)
					nextSplits.add(new Split(aWithT, wgt1));

				var bWithT = (BitSet) b.clone();
				bWithT.set(t);
				var wgt2 = Math.min(prev.weight, isolationIndex(t, bWithT, a, d));
				if (wgt2 > 0)
					nextSplits.add(new Split(bWithT, wgt2));
			}

			previousSplits = nextSplits;
			previousTaxa.set(t);
		}
		return previousSplits;
	}

	/** isolation index of the split A vs B (both within taxa 1..t, with t in A): the min over quartets. */
	private static double isolationIndex(int t, BitSet a, BitSet b, double[][] d) {
		var min = Double.MAX_VALUE;
		for (var i = 1; i <= t; i++) {
			if (a.get(i))
				for (var j = 1; j <= t; j++) {
					if (b.get(j))
						for (var k = j; k <= t; k++) {
							if (b.get(k)) {
								var val = quartetIndex(t, i, j, k, d);
								if (val < min) {
									if (val <= EPSILON)
										return 0;
									min = val;
								}
							}
						}
				}
		}
		return min;
	}

	/** the quartet isolation index of i,j vs k,m (1-based taxa). */
	private static double quartetIndex(int i, int j, int k, int m, double[][] d) {
		return 0.5 * (Math.max(di(d, i, k) + di(d, j, m), di(d, i, m) + di(d, j, k)) - di(d, i, j) - di(d, k, m));
	}

	private static double di(double[][] d, int a, int b) {
		return d[a - 1][b - 1];
	}

	private static BitSet complement(BitSet a, int ntax) {
		var result = new BitSet();
		for (var t = a.nextClearBit(1); t != -1 && t <= ntax; t = a.nextClearBit(t + 1))
			result.set(t);
		return result;
	}
}
