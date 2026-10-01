/*
 * AnalysisLog.java Copyright (C) 2025
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

import jloda.phylo.PhyloGraph;
import jloda.util.SetUtils;
import jloda.util.StringUtils;
import splitstree6.data.TaxaBlock;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Writes the analysis log of a RazorNet run: the input metric, the initial slack of every point under both
 * definitions (the slack of the Strict strategy and the positive slack of the Relaxed strategy), the recursion
 * trace of one restart and the distribution of the lengths over the restarts (both from {@link ExtensionTrace}),
 * the result, and what the strategies did.
 * <p>
 * Daniel Huson and Momoko Hayamizu, 8.2026
 */
public final class AnalysisLog {
	private AnalysisLog() {
	}

	/**
	 * @param inputDistances the ORIGINAL input distances, before quantization
	 * @param maxError       max realization error, or NaN to have this method compute it from the graph
	 * @param tracing        whether {@link ExtensionTrace} recorded this run
	 * @param runSummary     the caller's account of the run (strategies, restarts, winner, counters); may be empty
	 */
	public static void write(String logFile, String algorithmName, String optionsSummary, TaxaBlock taxaBlock,
							 double[][] inputDistances, int significantDigits, PhyloGraph graph,
							 double seconds, double maxError, boolean tracing, String runSummary) {
		try (var w = new BufferedWriter(new FileWriter(logFile))) {
			var n = inputDistances.length;
			w.write("RazorNet analysis log%n".formatted());
			w.write("=====================%n%n".formatted());
			w.write("algorithm  : %s%n".formatted(algorithmName));
			w.write("taxa       : %d%n".formatted(n));
			w.write("options    : %s%n%n".formatted(optionsSummary));

			var minD = Double.MAX_VALUE;
			var maxD = 0.0;
			for (var i = 0; i < n; i++)
				for (var j = 0; j < n; j++)
					if (i != j) {
						minD = Math.min(minD, inputDistances[i][j]);
						maxD = Math.max(maxD, inputDistances[i][j]);
					}
			var q = Quantization.apply(inputDistances, significantDigits);
			w.write("1. INPUT METRIC%n".formatted());
			w.write("   off-diagonal range      : %s .. %s%n".formatted(StringUtils.trim(minD), StringUtils.trim(maxD)));
			w.write("   quantized size          : %d (points at distance 0 are merged by Quantization)%n".formatted(q.size()));
			w.write("   1 quantized unit        : %s input units%n%n".formatted(StringUtils.trim(q.mapDistanceBack().applyAsDouble(1))));

			var M = new MutableD(q.quantMatrix());
			var size = M.size();
			w.write("2. INITIAL SLACK over S = all points, under both definitions%n".formatted());
			w.write("   c(x)  = min{ (D(x,y)+D(x,z)-D(y,z))/2 : y,z in S\\{x}, y!=z }   the Strict strategy's slack; 0 iff x lies on a geodesic of some pair%n".formatted());
			w.write("   c+(x) = the smallest POSITIVE value of that expression, 0 if none   the Relaxed strategy's slack%n".formatted());
			w.write("   (#zero: pairs {y,z} with x between them; #min: pairs attaining c(x); quantized units, c(x) also in input units)%n%n".formatted());
			w.write("   %-12s %10s %10s %8s %8s %10s%n".formatted("point", "c(x)", "c(x) input", "#zero", "#min", "c+(x)"));
			w.write("   %s%n".formatted("-".repeat(64)));
			var trueSlackPoints = 0;
			var positiveSlackPoints = 0;
			var trueValues = new ArrayList<Integer>();
			for (var x = 0; x < size; x++) {
				var min = Integer.MAX_VALUE;
				var minPositive = Integer.MAX_VALUE;
				var zero = 0;
				var attaining = 0;
				for (var y = 0; y < size; y++) {
					if (y == x) continue;
					for (var z = y + 1; z < size; z++) {
						if (z == x) continue;
						var value = (M.get(x, y) + M.get(x, z) - M.get(y, z)) / 2;
						if (value == 0) zero++;
						if (value > 0 && value < minPositive) minPositive = value;
						if (value < min) {
							min = value;
							attaining = 1;
						} else if (value == min)
							attaining++;
					}
				}
				var c = (min == Integer.MAX_VALUE ? 0 : min);
				var cPlus = (minPositive == Integer.MAX_VALUE ? 0 : minPositive);
				if (c > 0) {
					trueSlackPoints++;
					trueValues.add(c);
				}
				if (cPlus > 0)
					positiveSlackPoints++;
				w.write("   %-12s %10d %10s %8d %8d %10d%n".formatted(label(taxaBlock, x), c,
						StringUtils.trim(q.mapDistanceBack().applyAsDouble(c)), zero, attaining, cPlus));
			}
			trueValues.sort(Comparator.reverseOrder());
			w.write("   %s%n".formatted("-".repeat(64)));
			w.write("   slack points: %d of %d under c(x) (Strict), %d of %d under c+(x) (Relaxed)%n".formatted(
					trueSlackPoints, size, positiveSlackPoints, size));
			w.write("   c(x) values in the order the Strict Step 1 processes them: %s%n%n".formatted(trueValues));

			w.write("3. RECURSION TRACE (one restart; extendSubset calls in order)%n".formatted());
			var lines = ExtensionTrace.callLines();
			if (!tracing) {
				w.write("   not recorded (tracing was not enabled for this run)%n%n".formatted());
			} else if (lines.isEmpty()) {
				w.write("   no calls recorded%n%n".formatted());
			} else {
				w.write("   %s%n".formatted(ExtensionTrace.callHeader()));
				w.write("   %s%n".formatted("-".repeat(140)));
				for (var line : lines)
					w.write("   %s%n".formatted(line));
				w.write("   largest |D| over the whole run: %d%n%n".formatted(ExtensionTrace.maxMatrixSize()));
			}

			w.write("4. RESTART LENGTH DISTRIBUTION%n".formatted());
			var lengths = ExtensionTrace.restartLengths();
			if (!tracing || lengths.isEmpty()) {
				w.write("   not recorded%n%n".formatted());
			} else {
				var sorted = new ArrayList<>(lengths);
				Collections.sort(sorted);
				var sum = 0.0;
				for (var v : sorted) sum += v;
				w.write("   restarts recorded : %d%n".formatted(sorted.size()));
				w.write("   min / mean / max  : %s / %s / %s%n".formatted(StringUtils.trim(sorted.get(0)),
						StringUtils.trim(sum / sorted.size()), StringUtils.trim(sorted.get(sorted.size() - 1))));
				w.write("   distinct values   : %d%n".formatted(new TreeSet<>(sorted).size()));
				w.write("   histogram (value x count):%n".formatted());
				var counts = new TreeMap<Double, Integer>();
				for (var v : sorted) counts.merge(v, 1, Integer::sum);
				var shown = 0;
				for (var e : counts.entrySet()) {
					w.write("     %12s  %s (%d)%n".formatted(StringUtils.trim(e.getKey()),
							"*".repeat(Math.min(60, e.getValue())), e.getValue()));
					if (++shown >= 40) {
						w.write("     ... %d further values%n".formatted(counts.size() - shown));
						break;
					}
				}
				w.write("%n".formatted());
			}

			w.write("5. RESULT%n".formatted());
			if (graph != null) {
				var length = graph.edgeStream().mapToDouble(graph::getWeight).sum();
				var aux = (int) graph.nodeStream().filter(v -> !graph.hasTaxa(v)).count();
				w.write("   length              : %s%n".formatted(StringUtils.trim(length)));
				w.write("   nodes / edges       : %d / %d%n".formatted(graph.getNumberOfNodes(), graph.getNumberOfEdges()));
				w.write("   auxiliary vertices  : %d%n".formatted(aux));
			}
			var err = (Double.isNaN(maxError) && graph != null) ? maxRealizationError(graph, inputDistances) : maxError;
			w.write("   max realization err : %s%n".formatted(Double.isNaN(err) ? "not checked" : "%.3g".formatted(err)));
			if (!Double.isNaN(err) && err > 1e-6)
				w.write("   *** WARNING: this network does NOT realize the input metric ***%n".formatted());
			w.write("   seconds             : %.2f%n".formatted(seconds));
			if (runSummary != null && !runSummary.isBlank()) {
				w.write("%n6. STRATEGIES%n".formatted());
				for (var line : runSummary.strip().split("\\R"))
					w.write("   %s%n".formatted(line));
			}
			System.err.printf("Wrote analysis log: %s%n", logFile);
		} catch (IOException ex) {
			System.err.printf("Warning: could not write %s: %s%n", logFile, ex.getMessage());
		}
	}

	/**
	 * the largest absolute difference between an input distance and the corresponding shortest-path distance
	 */
	private static double maxRealizationError(PhyloGraph graph, double[][] inputDistances) {
		var nodes = new ArrayList<jloda.graph.Node>();
		var taxonOf = new ArrayList<Integer>();
		graph.nodeStream().filter(graph::hasTaxa).forEach(v -> {
			nodes.add(v);
			taxonOf.add(graph.getTaxon(v));
		});
		if (nodes.size() != inputDistances.length)
			return Double.POSITIVE_INFINITY; // wrong number of labeled nodes: not a valid realization
		var maxError = 0.0;
		for (var a = 0; a < nodes.size(); a++) {
			// one shortest-path computation per taxon
			var distanceFromA = splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances
					.graphDistances(graph, nodes.get(a));
			for (var b = a + 1; b < nodes.size(); b++) {
				var in = inputDistances[taxonOf.get(a) - 1][taxonOf.get(b) - 1];
				double out = distanceFromA.get(nodes.get(b));
				maxError = Math.max(maxError, Math.abs(in - out));
			}
		}
		return maxError;
	}

	private static String label(TaxaBlock taxaBlock, int index) {
		return (taxaBlock != null && index < taxaBlock.getNtax()) ? taxaBlock.getLabel(index + 1) : ("p" + index);
	}
}
