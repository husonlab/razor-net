/*
 * CompareStrategies.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.graph.Node;
import jloda.phylo.PhyloGraph;
import jloda.util.ProgramExecutorService;
import jloda.util.progress.ProgressSilent;
import razornet.algorithms.RazorNet;
import razornet.algorithms.utils.MatrixIO;
import splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;

import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Stream;

/**
 * Runs the three strategies of RazorNet (Strict, Relaxed and Both) on the same {@code .dist} instances and
 * reports, per instance, the lengths, the times and whether each network realizes its input. The exit status is 0
 * if every network realizes its input and {@code Both} is never longer than the shorter of the two single
 * strategies. By default the runs are restart-bound ({@code --maxSeconds 0}) and single-threaded, which makes the
 * comparison reproducible.
 * <p>
 * Usage: {@code CompareStrategies [--input DIR] [--maxSeconds N] [--restarts N] [--maxSize N] [--seed N] [--threads N]}
 *
 * @author Daniel Huson, 2026
 */
public class CompareStrategies {

	/** One run on one instance; {@code error} non-null means the run threw. */
	private record Result(double length, int nodes, int edges, int auxiliary, double maxDiff, double seconds,
						  String error) {
		boolean realized(double tolerance) {
			return error == null && maxDiff <= tolerance;
		}
	}

	public static void main(String[] args) throws Exception {
		razornet.algorithms.utils.AdvancedOptions.enable(); // every option stays reachable by name
		var input = "tests";
		var maxSeconds = 0;      // restart-bound by default: a budget makes the comparison timing-dependent
		var restarts = 100;
		var maxSize = Integer.MAX_VALUE;
		var seed = 0;
		var threads = 1;         // one thread, for reproducibility
		var realizeTol = 1e-3;

		for (var i = 0; i < args.length; i++) {
			try {
				switch (args[i]) {
					case "--input", "-i" -> input = args[++i];
					case "--maxSeconds", "-t" -> maxSeconds = Integer.parseInt(args[++i]);
					case "--restarts", "-r" -> restarts = Integer.parseInt(args[++i]);
					case "--maxSize", "-m" -> maxSize = Integer.parseInt(args[++i]);
					case "--seed", "-s" -> seed = Integer.parseInt(args[++i]);
					case "--threads" -> threads = Integer.parseInt(args[++i]);
					case "--realizeTol" -> realizeTol = Double.parseDouble(args[++i]);
					case "--help", "-h" -> {
						printUsage();
						return;
					}
					default -> {
						System.err.println("Unknown argument: " + args[i]);
						printUsage();
						System.exit(2);
					}
				}
			} catch (ArrayIndexOutOfBoundsException | NumberFormatException ex) {
				System.err.println("Bad argument near '" + args[Math.max(0, i)] + "': " + ex.getMessage());
				System.exit(2);
			}
		}

		if (threads > 0)
			ProgramExecutorService.setNumberOfCoresToUse(threads);
		final var tolerance = realizeTol;

		final List<Path> files;
		try (Stream<Path> walk = Files.walk(Path.of(input))) {
			files = walk.filter(p -> p.getFileName().toString().endsWith(".dist"))
					.sorted(Comparator.comparing(Path::toString)).toList();
		}
		if (files.isEmpty()) {
			System.err.println("No .dist files found under: " + input);
			System.exit(2);
		}

		final var strategies = List.of(RazorNet.Strategy.Strict, RazorNet.Strategy.Relaxed, RazorNet.Strategy.Both);

		final var report = System.out;
		report.printf("RazorNet strategies: %d instances under '%s' (maxSeconds=%d, restarts=%d, seed=%d, threads=%d)%n",
				files.size(), input, maxSeconds, restarts, seed, threads);
		report.println("Verdict OK: every network realizes its input and Both is no longer than the shorter single strategy.");
		report.println();
		report.printf("%-11s %5s  %12s %12s %12s  %9s %9s %9s  %-8s %-9s %s%n",
				"TYPE", "N", "STRICT", "RELAXED", "BOTH", "T-STRICT", "T-RELAXED", "T-BOTH", "REALIZED", "VERDICT", "SHORTEST");
		report.println("-".repeat(118));

		var ok = 0;
		var bad = 0;
		var skipped = 0;
		var notes = new ArrayList<String>();

		for (var file : files) {
			var label = file.getFileName().toString().replace(".dist", "");
			var type = label.replaceAll("-n\\d+$", "");
			final double[][] matrix;
			try {
				matrix = MatrixIO.read(file.toString());
			} catch (Exception ex) {
				// a file that cannot be read is reported and skipped
				skipped++;
				notes.add("%s: not read (%s)".formatted(label, ex.getMessage()));
				continue;
			}
			var n = matrix.length;
			if (n > maxSize) {
				skipped++;
				continue;
			}

			var results = new ArrayList<Result>();
			for (var strategy : strategies) {
				var algorithm = new RazorNet();
				algorithm.setOptionStrategy(strategy);
				algorithm.optionMaxSecondsProperty().set(maxSeconds);
				algorithm.optionRandomRestartsProperty().set(restarts);
				algorithm.optionRandomSeedProperty().set(seed);
				results.add(run(algorithm, matrix));
			}
			var strict = results.get(0);
			var relaxed = results.get(1);
			var both = results.get(2);

			var realizedFlags = new StringBuilder();
			for (var i = 0; i < strategies.size(); i++)
				realizedFlags.append(results.get(i).realized(tolerance) ? strategies.get(i).name().charAt(0) : '-');
			var allRealized = results.stream().allMatch(r -> r.realized(tolerance));

			var shortestSingle = Math.min(strict.length(), relaxed.length());
			// restart-bound, Both is exactly the minimum; under a time budget the phases share the clock
			var bothNotLonger = both.error() == null && both.length() <= shortestSingle + 1e-9;

			String verdict;
			if (!allRealized) {
				verdict = "UNREALIZED";
				bad++;
				notes.add("%s: a network does not realize its input (%s)".formatted(label, realizedFlags));
			} else if (!bothNotLonger) {
				verdict = maxSeconds > 0 ? "BUDGET" : "WORSE";
				bad++;
				notes.add("%s: Both %s is longer than the shorter single strategy %s%s".formatted(label,
						lengthString(both), lengthString(shortestSingle),
						maxSeconds > 0 ? " (budget-bound run: the phases shared the clock)" : ""));
			} else {
				verdict = "OK";
				ok++;
			}

			var shortest = new ArrayList<String>();
			var best = results.stream().filter(r -> r.error() == null).mapToDouble(Result::length).min().orElse(Double.NaN);
			for (var i = 0; i < strategies.size(); i++)
				if (results.get(i).error() == null && Math.abs(results.get(i).length() - best) <= 1e-9)
					shortest.add(strategies.get(i).name());

			report.printf("%-11s %5d  %12s %12s %12s  %8.1fs %8.1fs %8.1fs  %-8s %-9s %s%n",
					type, n, lengthString(strict), lengthString(relaxed), lengthString(both),
					strict.seconds(), relaxed.seconds(), both.seconds(), realizedFlags, verdict,
					shortest.size() == strategies.size() ? "tie" : String.join("=", shortest));
		}

		report.println("-".repeat(118));
		report.printf("%d ok, %d not ok%s.%n", ok, bad, skipped > 0 ? ", " + skipped + " skipped (--maxSize or unreadable)" : "");
		if (!notes.isEmpty()) {
			report.println("Details:");
			for (var s : notes)
				report.println("  - " + s);
		}
		System.exit(bad == 0 ? 0 : 1);
	}

	/** Runs one configured algorithm on one matrix, muting the (chatty) algorithm output. */
	private static Result run(RazorNet algorithm, double[][] matrix) {
		var n = matrix.length;

		var taxa = new TaxaBlock();
		taxa.setNtax(n);
		var names = new ArrayList<String>();
		for (var t = 1; t <= n; t++)
			names.add("t" + t);
		taxa.addTaxaByNames(names);

		var distances = new DistancesBlock();
		distances.setNtax(n);
		distances.set(matrix);

		var network = new NetworkBlock();

		var start = System.currentTimeMillis();
		String error = null;
		var savedOut = System.out;
		var savedErr = System.err;
		try (var mute = new PrintStream(OutputStream.nullOutputStream())) {
			System.setOut(mute);
			System.setErr(mute);
			algorithm.compute(new ProgressSilent(), taxa, distances, network);
		} catch (Throwable ex) {
			error = ex.getClass().getSimpleName() + (ex.getMessage() != null ? ": " + ex.getMessage() : "");
		} finally {
			System.setOut(savedOut);
			System.setErr(savedErr);
		}
		var seconds = (System.currentTimeMillis() - start) / 1000.0;

		if (error != null)
			return new Result(Double.NaN, 0, 0, 0, Double.NaN, seconds, error);

		var graph = network.getGraph();
		var length = graph.edgeStream().mapToDouble(graph::getWeight).sum();
		var nodes = graph.getNumberOfNodes();
		var edges = graph.getNumberOfEdges();
		var auxiliary = (int) graph.nodeStream().filter(v -> !graph.hasTaxa(v)).count();
		return new Result(length, nodes, edges, auxiliary, realizationMaxDiff(graph, matrix), seconds, null);
	}

	/** Largest absolute difference between an input distance and the corresponding shortest-path distance. */
	private static double realizationMaxDiff(PhyloGraph graph, double[][] d) {
		var nodes = new ArrayList<Node>();
		var taxonOf = new ArrayList<Integer>();
		graph.nodeStream().filter(graph::hasTaxa).forEach(v -> {
			nodes.add(v);
			taxonOf.add(graph.getTaxon(v));
		});
		if (nodes.size() != d.length)
			return Double.POSITIVE_INFINITY; // wrong number of labeled nodes => not a valid realization
		var maxDiff = 0.0;
		for (var a = 0; a < nodes.size(); a++) {
			// one shortest-path computation per taxon
			var distanceFromA = CheckPairwiseDistances.graphDistances(graph, nodes.get(a));
			for (var b = a + 1; b < nodes.size(); b++) {
				var in = d[taxonOf.get(a) - 1][taxonOf.get(b) - 1];
				double out = distanceFromA.get(nodes.get(b));
				maxDiff = Math.max(maxDiff, Math.abs(in - out));
			}
		}
		return maxDiff;
	}

	private static String lengthString(Result r) {
		return r.error() != null ? "ERROR" : lengthString(r.length());
	}

	private static String lengthString(double length) {
		return length == Math.rint(length) ? "%d".formatted((long) length) : "%.4f".formatted(length);
	}

	private static void printUsage() {
		System.err.println("""
				CompareStrategies - run RazorNet's Strict, Relaxed and Both strategies over the same .dist
				instances and compare them

				Options:
				  -i, --input DIR       directory to scan recursively for .dist files (default: tests)
				  -t, --maxSeconds N    per-run time budget (default: 0 = none, so the comparison is restart-bound)
				  -r, --restarts N      random restarts per strategy (default: 100)
				  -m, --maxSize N       skip instances with more than N taxa (for quick runs)
				  -s, --seed N          random seed (default: 0)
				      --threads N       size of the parallel pool (default: 1, for reproducibility)
				      --realizeTol F    tolerance for the realization check (default: 1e-3)
				  -h, --help            show this help

				Exit status is 0 iff every network realizes its input and Both is never longer than the shorter of
				the two single strategies. Which strategy is shortest where is reported, not judged.""");
	}
}
