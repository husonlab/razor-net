/*
 * RunTestSuite.java Copyright (C) 2026 Daniel H. Huson
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
 * Runs RazorNet over a directory of generated {@code .dist} test instances (see
 * {@link GenerateTestInstances}) and turns each into a pass/fail check:
 * <ul>
 *   <li><b>realization</b> (all tiers): every pairwise shortest-path distance in RazorNet's network must
 *       equal the input distance (within tolerance). This is the fundamental correctness requirement.</li>
 *   <li><b>length</b> (tier-specific, read from each instance's {@code # ... referenceLength=...} header):
 *       <ul>
 *         <li>tree &mdash; the reference is the optimum, so realized length must be &le; reference*(1+tol);</li>
 *         <li>general &mdash; the reference is a realizing graph (upper bound), so length must be &le; reference*(1+tol);</li>
 *         <li>circle &mdash; the reference is a split weight, not a network length, so length is reported for
 *             information only and does not affect pass/fail.</li>
 *       </ul></li>
 * </ul>
 * Exits with status 0 iff every instance passes, so it doubles as a CI regression check.
 * <p>
 * Usage: {@code RunTestSuite [--input DIR] [--maxSeconds N] [--restarts N] [--lengthTolerance F] [--maxSize N]}
 *
 * @author Daniel Huson, 2026
 */
public class RunTestSuite {

	private record Meta(String type, int n, double referenceLength, boolean referenceIsOptimal) {
	}

	public static void main(String[] args) throws Exception {
		razornet.algorithms.utils.AdvancedOptions.enable(); // every option stays reachable by name
		var input = "tests";
		var maxSeconds = 60;
		var restarts = 100;
		var lengthTol = 0.02;    // 2% slack on length checks
		var realizeTol = 1e-3;   // absolute tolerance on realized distances
		var maxSize = Integer.MAX_VALUE;

		for (var i = 0; i < args.length; i++) {
			try {
				switch (args[i]) {
					case "--input", "-i" -> input = args[++i];
					case "--maxSeconds", "-t" -> maxSeconds = Integer.parseInt(args[++i]);
					case "--restarts", "-r" -> restarts = Integer.parseInt(args[++i]);
					case "--lengthTolerance", "-l" -> lengthTol = Double.parseDouble(args[++i]);
					case "--maxSize", "-m" -> maxSize = Integer.parseInt(args[++i]);
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

		final List<Path> files;
		try (Stream<Path> walk = Files.walk(Path.of(input))) {
			// only the generated test instances, which carry a '# type=... referenceLength=...' header
			files = walk.filter(p -> p.getFileName().toString().endsWith(".dist"))
					.filter(RunTestSuite::isTestInstance)
					.sorted(Comparator.comparing(Path::toString)).toList();
		}
		if (files.isEmpty()) {
			System.err.println("No .dist files found under: " + input);
			System.exit(2);
		}

		// keep a handle to the real stdout: RazorNet is chatty, so we mute stdout/stderr during each run
		final var report = System.out;
		report.printf("RazorNet test suite: %d instances under '%s' (maxSeconds=%d, restarts=%d, lengthTol=%.0f%%)%n%n",
				files.size(), input, maxSeconds, restarts, lengthTol * 100);
		report.printf("%-8s %5s  %-8s %-9s %11s %11s %7s  %-22s %6s%n",
				"TYPE", "N", "REALIZED", "maxDiff", "LENGTH", "REFERENCE", "RATIO", "LENGTH-CHECK", "TIME");
		report.println("-".repeat(100));

		var passed = 0;
		var failed = 0;
		var skipped = 0;
		var failures = new ArrayList<String>();

		for (var file : files) {
			var meta = parseMeta(file);
			if (meta.n() > maxSize) {
				skipped++;
				continue;
			}

			var matrix = MatrixIO.read(file.toString());
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
			var razor = new RazorNet();
			razor.optionMaxSecondsProperty().set(maxSeconds);
			razor.optionRandomRestartsProperty().set(restarts);

			var start = System.currentTimeMillis();
			String error = null;
			var savedOut = System.out;
			var savedErr = System.err;
			try (var mute = new PrintStream(OutputStream.nullOutputStream())) {
				System.setOut(mute);
				System.setErr(mute);
				razor.compute(new ProgressSilent(), taxa, distances, network);
			} catch (Throwable ex) {
				error = ex.getClass().getSimpleName() + (ex.getMessage() != null ? ": " + ex.getMessage() : "");
			} finally {
				System.setOut(savedOut);
				System.setErr(savedErr);
			}
			var elapsed = (System.currentTimeMillis() - start) / 1000.0;

			if (error != null) {
				report.printf("%-8s %5d  %-8s %-9s %11s %11s %7s  %-22s %5.1fs%n",
						meta.type(), meta.n(), "ERROR", "-", "-", refString(meta), "-", trunc(error, 22), elapsed);
				failed++;
				failures.add("%s-n%04d: %s".formatted(meta.type(), meta.n(), error));
				continue;
			}

			var graph = network.getGraph();
			var maxDiff = realizationMaxDiff(graph, matrix);
			var realized = maxDiff <= realizeTol;
			var length = graph.edgeStream().mapToDouble(graph::getWeight).sum();

			// length verdict
			boolean lengthOk;
			String lengthCheck;
			String ratio;
			if (Double.isNaN(meta.referenceLength()) || meta.referenceLength() <= 0) {
				lengthOk = true;
				lengthCheck = "n/a (no reference)";
				ratio = "-";
			} else if (meta.type().equals("circle")) {
				lengthOk = true; // split weight is not a network length
				lengthCheck = "n/a (split weight)";
				ratio = "-";
			} else {
				var r = length / meta.referenceLength();
				ratio = "%.3f".formatted(r);
				var excess = 100.0 * (r - 1.0);
				lengthOk = length <= meta.referenceLength() * (1.0 + lengthTol);
				if (meta.referenceIsOptimal())
					lengthCheck = lengthOk ? (excess <= 1e-6 ? "OK (optimal)" : "OK (+%.2f%%)".formatted(excess))
							: "SUBOPTIMAL (+%.2f%%)".formatted(excess);
				else
					lengthCheck = lengthOk ? "OK (<= reference)" : "OVER (+%.2f%%)".formatted(excess);
			}

			var ok = realized && lengthOk;
			report.printf("%-8s %5d  %-8s %-9.2e %11s %11s %7s  %-22s %5.1fs%n",
					meta.type(), meta.n(), realized ? "yes" : "NO", maxDiff,
					"%.1f".formatted(length), refString(meta), ratio, lengthCheck, elapsed);

			if (ok) {
				passed++;
			} else {
				failed++;
				var why = new ArrayList<String>();
				if (!realized) why.add("not realized (maxDiff=%.3g)".formatted(maxDiff));
				if (!lengthOk) why.add(lengthCheck);
				failures.add("%s-n%04d: %s".formatted(meta.type(), meta.n(), String.join(", ", why)));
			}
		}

		report.println("-".repeat(100));
		report.printf("%d passed, %d failed%s.%n", passed, failed, skipped > 0 ? ", " + skipped + " skipped (--maxSize)" : "");
		if (!failures.isEmpty()) {
			report.println("Failures:");
			for (var f : failures)
				report.println("  - " + f);
		}
		System.exit(failed == 0 ? 0 : 1);
	}

	/** Largest absolute difference between an input distance and the corresponding shortest-path distance in RazorNet's network. */
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

	/** True iff the file carries the '# type=...' header that marks a generated test instance. */
	private static boolean isTestInstance(Path file) {
		try {
			for (var line : Files.readAllLines(file)) {
				line = line.trim();
				if (line.startsWith("# type="))
					return true;
				if (!line.isEmpty() && !line.startsWith("#"))
					return false;  // reached the size line without a header
			}
		} catch (java.io.IOException ignored) {
		}
		return false;
	}

	/** Reads the '# type=... n=... referenceLength=... referenceIsOptimal=...' header, falling back to the parent directory name. */
	private static Meta parseMeta(Path file) throws java.io.IOException {
		String type = null;
		var n = -1;
		var referenceLength = Double.NaN;
		var optimal = false;
		for (var line : Files.readAllLines(file)) {
			line = line.trim();
			if (line.startsWith("# type=")) {
				for (var tok : line.substring(1).trim().split("\\s+")) {
					var kv = tok.split("=", 2);
					if (kv.length != 2) continue;
					switch (kv[0]) {
						case "type" -> type = kv[1];
						case "n" -> n = Integer.parseInt(kv[1]);
						case "referenceLength" -> referenceLength = Double.parseDouble(kv[1]);
						case "referenceIsOptimal" -> optimal = Boolean.parseBoolean(kv[1]);
					}
				}
				break;
			}
			if (!line.isEmpty() && !line.startsWith("#"))
				break; // hit the size line without a header
		}
		if (type == null) {
			var parent = file.getParent() == null ? "" : file.getParent().getFileName().toString();
			type = parent.replace("-distances", "");
		}
		if (n < 0) {
			try {
				n = MatrixIO.read(file.toString()).length;
			} catch (Exception ignored) {
			}
		}
		return new Meta(type, n, referenceLength, optimal);
	}

	private static String refString(Meta meta) {
		if (Double.isNaN(meta.referenceLength()))
			return "-";
		return "%.0f%s".formatted(meta.referenceLength(),
				meta.type().equals("circle") ? " (sw)" : meta.referenceIsOptimal() ? " (opt)" : "");
	}

	private static String trunc(String s, int max) {
		return s.length() <= max ? s : s.substring(0, max - 1) + "…";
	}

	private static void printUsage() {
		System.err.println("""
				RunTestSuite - run RazorNet over a directory of .dist test instances and check the results

				Options:
				  -i, --input DIR           directory to scan recursively for .dist files (default: tests)
				  -t, --maxSeconds N        per-instance time budget passed to RazorNet (default: 60)
				  -r, --restarts N          random restarts per instance (default: 100)
				  -l, --lengthTolerance F   fractional slack on length checks, e.g. 0.02 = 2% (default: 0.02)
				  -m, --maxSize N           skip instances with more than N taxa (for quick runs)
				  -h, --help                show this help

				Checks per instance: distances realized exactly (all tiers); realized length <= reference*(1+tol)
				for tree (optimum) and general (upper bound); circle length is reported only. Exit status is 0
				iff all (non-skipped) instances pass.""");
	}
}
