/*
 * QuartetDisplay.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.testsuite;

import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Random;

/**
 * How strongly a network displays each of the three splits of a quartet of taxa, and how that compares with the
 * truth.
 * <p>
 * In a graph, the support of the quartet split ab|cd is the distance between the shortest paths from a to b and
 * the shortest paths from c to d: the length of the part of the graph that separates the two pairs, and zero if a
 * shortest path between a and b meets one between c and d. In a tree this is the length of the path that joins
 * the path from a to b to the path from c to d. In a set of weighted splits, the support is the total weight of
 * the splits that have a and b on one side and c and d on the other.
 *
 * @author Daniel Huson, 2026
 */
public final class QuartetDisplay {
	/** reports, for four taxa given by their indices, the supports of the splits ab|cd, ac|bd and ad|bc */
	public interface Support {
		void get(int a, int b, int c, int d, double[] supports);
	}

	/** the supports displayed by a graph, for the taxa in the given order; a taxon that is not in the graph has none */
	public static Support ofGraph(TaxonGraph graph, List<String> taxa) {
		var n = taxa.size();
		var nodes = graph.getNumberOfNodes();
		var node = new int[n];
		var distance = new double[n][];
		for (var i = 0; i < n; i++) {
			node[i] = graph.nodeOf(taxa.get(i));
			if (node[i] >= 0)
				distance[i] = graph.distancesFrom(node[i]);
		}
		// for each pair of taxa i<j: the nodes on the shortest paths between them, and the distance of every node
		// from these
		var interval = new int[n][n][];
		var fromInterval = new float[n][n][];
		for (var i = 0; i < n; i++) {
			for (var j = i + 1; j < n; j++) {
				if (node[i] < 0 || node[j] < 0)
					continue;
				var total = distance[i][node[j]];
				var tolerance = 1e-9 * Math.max(1.0, total);
				var members = new int[nodes];
				var count = 0;
				for (var v = 0; v < nodes; v++)
					if (distance[i][v] + distance[j][v] <= total + tolerance)
						members[count++] = v;
				interval[i][j] = Arrays.copyOf(members, count);
				fromInterval[i][j] = graph.distancesFrom(interval[i][j]);
			}
		}
		return (a, b, c, d, supports) -> {
			if (node[a] < 0 || node[b] < 0 || node[c] < 0 || node[d] < 0) {
				supports[0] = supports[1] = supports[2] = 0;
				return;
			}
			supports[0] = between(interval, fromInterval, a, b, c, d);
			supports[1] = between(interval, fromInterval, a, c, b, d);
			supports[2] = between(interval, fromInterval, a, d, b, c);
		};
	}

	/** the distance between the shortest paths from a to b and the shortest paths from c to d */
	private static double between(int[][][] interval, float[][][] fromInterval, int a, int b, int c, int d) {
		var from = a < b ? fromInterval[a][b] : fromInterval[b][a];
		var best = Float.POSITIVE_INFINITY;
		for (var v : (c < d ? interval[c][d] : interval[d][c])) {
			if (from[v] < best) {
				best = from[v];
				if (best <= 0)
					break;
			}
		}
		return best;
	}

	/**
	 * the supports displayed by a set of weighted splits
	 *
	 * @param sides   for each split, the indices of the taxa on one of its sides
	 * @param weights the weights of the splits
	 */
	public static Support ofSplits(List<BitSet> sides, double[] weights, int numberOfTaxa) {
		// for each taxon, the splits that have it on the given side
		var splitsOf = new BitSet[numberOfTaxa];
		for (var i = 0; i < numberOfTaxa; i++)
			splitsOf[i] = new BitSet(sides.size());
		for (var s = 0; s < sides.size(); s++)
			for (var i = sides.get(s).nextSetBit(0); i >= 0; i = sides.get(s).nextSetBit(i + 1))
				if (i < numberOfTaxa)
					splitsOf[i].set(s);
		return (a, b, c, d, supports) -> {
			supports[0] = separating(splitsOf, weights, a, b, c, d);
			supports[1] = separating(splitsOf, weights, a, c, b, d);
			supports[2] = separating(splitsOf, weights, a, d, b, c);
		};
	}

	/** the total weight of the splits with a and b on one side and c and d on the other */
	private static double separating(BitSet[] splitsOf, double[] weights, int a, int b, int c, int d) {
		var total = 0.0;
		var one = (BitSet) splitsOf[a].clone();
		one.and(splitsOf[b]);
		one.andNot(splitsOf[c]);
		one.andNot(splitsOf[d]);
		for (var s = one.nextSetBit(0); s >= 0; s = one.nextSetBit(s + 1))
			total += weights[s];
		var other = (BitSet) splitsOf[c].clone();
		other.and(splitsOf[d]);
		other.andNot(splitsOf[a]);
		other.andNot(splitsOf[b]);
		for (var s = other.nextSetBit(0); s >= 0; s = other.nextSetBit(s + 1))
			total += weights[s];
		return total;
	}

	/**
	 * the supports that the distances themselves give to the splits of a quartet: for ab|cd, half the amount by which
	 * the larger of d(a,c)+d(b,d) and d(a,d)+d(b,c) exceeds d(a,b)+d(c,d), if it does. No graph that displays these
	 * distances can display less.
	 */
	public static Support ofDistances(double[][] distances) {
		return (a, b, c, d, supports) -> {
			var abcd = distances[a][b] + distances[c][d];
			var acbd = distances[a][c] + distances[b][d];
			var adbc = distances[a][d] + distances[b][c];
			supports[0] = Math.max(0, 0.5 * (Math.max(acbd, adbc) - abcd));
			supports[1] = Math.max(0, 0.5 * (Math.max(abcd, adbc) - acbd));
			supports[2] = Math.max(0, 0.5 * (Math.max(abcd, acbd) - adbc));
		};
	}

	/** the sum of the given supports */
	public static Support sum(List<Support> parts) {
		return (a, b, c, d, supports) -> {
			var total = new double[3];
			for (var part : parts) {
				part.get(a, b, c, d, supports);
				for (var j = 0; j < 3; j++)
					total[j] += supports[j];
			}
			System.arraycopy(total, 0, supports, 0, 3);
		};
	}

	/**
	 * The comparison of the supports displayed by a network with the true supports, summed over quartets. For a
	 * quartet, the split with the largest true support is its main split; support for another split is conflict.
	 *
	 * @param quartets         quartets compared
	 * @param trueSupport      true support, all splits
	 * @param shownSupport     displayed support, all splits
	 * @param sharedSupport    for each split the smaller of its true and its displayed support, summed
	 * @param trueConflict     true support of the splits that are not the main split of their quartet
	 * @param shownConflict    displayed support of these splits
	 * @param sharedConflict   for each of these splits the smaller of its true and its displayed support, summed
	 * @param falseSupport     displayed support of the splits that have no true support at all
	 * @param trueSplits       splits with true support that are not the main split of their quartet
	 * @param foundSplits      of these, the network displays support
	 * @param falseSplits      splits without true support for which the network displays support
	 * @param conflictQuartets quartets with true conflict
	 * @param conflictDetected of these, the network displays support for a true conflicting split
	 * @param simpleQuartets   quartets without true conflict
	 * @param simpleClean      of these, the network displays no support for a conflicting split
	 */
	public record Scores(long quartets, double trueSupport, double shownSupport, double sharedSupport, double trueConflict,
						 double shownConflict, double sharedConflict, double falseSupport, long trueSplits, long foundSplits,
						 long falseSplits, long conflictQuartets, long conflictDetected, long simpleQuartets, long simpleClean) {
		public static final String HEADER = "quartets\ttrueSupport\tshownSupport\tsharedSupport\ttrueConflict\tshownConflict"
											+ "\tsharedConflict\tfalseSupport\ttrueSplits\tfoundSplits\tfalseSplits"
											+ "\tconflictQuartets\tconflictDetected\tsimpleQuartets\tsimpleClean";

		public String toTsv() {
			return "%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%d\t%d\t%d\t%d\t%d\t%d\t%d".formatted(quartets, trueSupport,
					shownSupport, sharedSupport, trueConflict, shownConflict, sharedConflict, falseSupport, trueSplits,
					foundSplits, falseSplits, conflictQuartets, conflictDetected, simpleQuartets, simpleClean);
		}
	}

	/**
	 * compares the supports displayed by a network with the true supports on all quartets of the taxa in use, or on
	 * a random sample of them
	 *
	 * @param use         the taxa to draw quartets from
	 * @param maxQuartets if there are more quartets than this, this many are drawn at random
	 */
	public static Scores compare(Support truth, Support network, boolean[] use, long maxQuartets, long seed) {
		var taxa = new int[use.length];
		var n = 0;
		for (var i = 0; i < use.length; i++)
			if (use[i])
				taxa[n++] = i;
		var counts = new long[8];
		var sums = new double[7];
		if (n >= 4) {
			var total = (long) n * (n - 1) * (n - 2) * (n - 3) / 24;
			if (total <= maxQuartets) {
				for (var a = 0; a < n; a++)
					for (var b = a + 1; b < n; b++)
						for (var c = b + 1; c < n; c++)
							for (var d = c + 1; d < n; d++)
								add(truth, network, taxa[a], taxa[b], taxa[c], taxa[d], counts, sums);
			} else {
				var random = new Random(seed);
				var q = new int[4];
				for (var k = 0L; k < maxQuartets; k++) {
					for (var i = 0; i < 4; i++) {
						boolean repeated;
						do {
							q[i] = taxa[random.nextInt(n)];
							repeated = false;
							for (var j = 0; j < i; j++)
								repeated |= q[j] == q[i];
						} while (repeated);
					}
					add(truth, network, q[0], q[1], q[2], q[3], counts, sums);
				}
			}
		}
		return new Scores(counts[0], sums[0], sums[1], sums[2], sums[3], sums[4], sums[5], sums[6], counts[5], counts[6],
				counts[7], counts[1], counts[2], counts[3], counts[4]);
	}

	private static final double EPSILON = 1e-9;

	private static void add(Support truth, Support network, int a, int b, int c, int d, long[] counts, double[] sums) {
		var t = new double[3];
		var s = new double[3];
		truth.get(a, b, c, d, t);
		network.get(a, b, c, d, s);
		var main = t[0] >= t[1] && t[0] >= t[2] ? 0 : (t[1] >= t[2] ? 1 : 2);
		if (t[main] <= EPSILON)
			main = -1; // the truth does not resolve this quartet: all displayed support is conflict
		var hasConflict = false;
		var detected = false;
		var clean = true;
		counts[0]++;
		for (var j = 0; j < 3; j++) {
			sums[0] += t[j];
			sums[1] += s[j];
			sums[2] += Math.min(t[j], s[j]);
			if (j != main) {
				sums[3] += t[j];
				sums[4] += s[j];
				sums[5] += Math.min(t[j], s[j]);
				if (t[j] > EPSILON) {
					hasConflict = true;
					counts[5]++;
					if (s[j] > EPSILON) {
						detected = true;
						counts[6]++;
					}
				} else if (s[j] > EPSILON) {
					sums[6] += s[j];
					counts[7]++;
				}
				if (s[j] > EPSILON)
					clean = false;
			}
		}
		if (hasConflict) {
			counts[1]++;
			if (detected)
				counts[2]++;
		} else {
			counts[3]++;
			if (clean)
				counts[4]++;
		}
	}
}
