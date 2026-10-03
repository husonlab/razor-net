/*
 * CompareDistanceMethods.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.phylo.PhyloSplitsGraph;
import jloda.util.Basic;
import jloda.util.FileUtils;
import jloda.util.ProgramExecutorService;
import jloda.util.progress.ProgressSilent;
import splitstree6.algorithms.distances.distances2splits.NeighborNet;
import splitstree6.algorithms.network.network2network.SimplifyFilter;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.SplitsBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.ImportManager;
import splitstree6.io.writers.ExportManager;
import splitstree6.layout.splits.algorithms.ConvexHull;
import splitstree6.layout.splits.algorithms.EqualAngle;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * Compares NeighborNet with RazorNet, and with RazorNet followed by the Simplify filter, on distance matrices,
 * without a GUI. For each input and method it reports one row: the size and length of the network, how well the
 * distances it displays agree with the input, the time and, if the trees that generated the distances are given,
 * how the quartet supports it displays compare with those of the trees ({@link QuartetDisplay}). The method
 * {@code Distances} reports the quartet supports that the input distances give by themselves.
 *
 * <pre>
 * CompareDistanceMethods [--truth name] [--methods Distances,NeighborNet,RazorNet,RazorNet+Simplify] [--strategy S]
 *                        [--restarts N] [--maxSeconds N] [--seed N] [--threads N] [--maxQuartets N]
 *                        [--set option=value]* [--tag text] [--report file] [--networks directory] [--silent] file...
 * </pre>
 * {@code --set} sets a further option of RazorNet by name, and {@code --tag} is appended to its method name in the
 * report. {@code --timeLimit} ends the Simplify filter after that many seconds, which is then reported as a timeout.
 * {@code --truth} names a file next to each input that holds the generating trees in Newick format, one per line,
 * with the number of mutations on each edge as its length.
 *
 * @author Daniel Huson, 2026
 */
public class CompareDistanceMethods {
	private static final String HEADER = "file\tntax\tmethod\tstatus\tseconds\tnodes\tedges\tcycles\tlength\tsplits\tfit\tminRatio\tmaxRatio";

	/**
	 * a network as a graph, with the number of splits behind it (0 if it is not a split network) and its time; the
	 * graph is null for the row that reports the input distances themselves
	 */
	private record Outcome(String method, TaxonGraph graph, QuartetDisplay.Support display, int splits, double seconds) {
	}

	public static void main(String[] args) throws Exception {
		razornet.algorithms.utils.AdvancedOptions.enable();
		var truthName = "";
		var methods = List.of("Distances", "NeighborNet", "RazorNet", "RazorNet+Simplify");
		var strategy = "";
		var restarts = 100;
		var maxSeconds = 0;
		var seed = 0;
		var threads = 1;
		var maxQuartets = 100000L;
		var timeLimit = 0;
		var reportFile = "";
		var networksDirectory = "";
		var silent = false;
		var tag = "";
		var settings = new java.util.LinkedHashMap<String, String>(); // further options of RazorNet, by name
		var inputs = new ArrayList<String>();
		for (var i = 0; i < args.length; i++) {
			switch (args[i]) {
				case "--truth" -> truthName = args[++i];
				case "--tag" -> tag = args[++i];
				case "--timeLimit" -> timeLimit = Integer.parseInt(args[++i]);
				case "--set" -> {
					var assignment = args[++i].split("=", 2);
					settings.put(assignment[0], assignment[1]);
				}
				case "--methods" -> methods = List.of(args[++i].split(","));
				case "--strategy" -> strategy = args[++i];
				case "--restarts", "-r" -> restarts = Integer.parseInt(args[++i]);
				case "--maxSeconds", "-t" -> maxSeconds = Integer.parseInt(args[++i]);
				case "--seed" -> seed = Integer.parseInt(args[++i]);
				case "--threads" -> threads = Integer.parseInt(args[++i]);
				case "--maxQuartets" -> maxQuartets = Long.parseLong(args[++i]);
				case "--report" -> reportFile = args[++i];
				case "--networks" -> networksDirectory = args[++i];
				case "--silent", "-s" -> silent = true;
				default -> {
					if (args[i].startsWith("-")) {
						System.err.println("Unknown argument: " + args[i]);
						System.exit(2);
					}
					inputs.add(args[i]);
				}
			}
		}
		if (inputs.isEmpty()) {
			System.err.println("Usage: CompareDistanceMethods [--truth name] [--methods Distances,NeighborNet,RazorNet,RazorNet+Simplify] [--strategy S] [--restarts N] [--maxSeconds N] [--seed N] [--threads N] [--maxQuartets N] [--report file] [--networks directory] [--silent] file...");
			System.exit(2);
		}
		ProgramExecutorService.setNumberOfCoresToUse(threads);
		if (!networksDirectory.isEmpty())
			Files.createDirectories(new File(networksDirectory).toPath());
		var out = reportFile.isEmpty() ? System.out : new PrintStream(reportFile);
		if (silent)
			Basic.hideSystemErr();
		out.println(HEADER + (truthName.isEmpty() ? "" : "\t" + QuartetDisplay.Scores.HEADER));

		for (var input : inputs) {
			var taxa = new TaxaBlock();
			var distances = new DistancesBlock();
			var reader = ImportManager.getInstance().getReader(input);
			if (reader == null || !DistancesBlock.class.isAssignableFrom(reader.getToClass()))
				throw new IOException("Not distances: " + input);
			((splitstree6.io.utils.DataReaderBase<DistancesBlock>) reader).read(new ProgressSilent(), input, taxa, distances);
			var labels = taxa.getLabels();
			var n = labels.size();

			// the truth: the quartet supports of the generating trees
			QuartetDisplay.Support truth = null;
			var use = new boolean[n];
			java.util.Arrays.fill(use, true);
			if (!truthName.isEmpty()) {
				var displays = new ArrayList<QuartetDisplay.Support>();
				for (var line : Files.readAllLines(new File(new File(input).getParentFile(), truthName).toPath())) {
					if (!line.isBlank() && !line.startsWith("#"))
						displays.add(QuartetDisplay.ofGraph(TaxonGraph.ofNewick(line), labels));
				}
				truth = QuartetDisplay.sum(displays);
			}

			var outcomes = new ArrayList<Outcome>();
			NetworkBlock razorNetwork = null;
			var razorSeconds = 0.0;
			for (var method : methods) {
				var start = System.currentTimeMillis();
				switch (method) {
					case "Distances" -> {
						var matrix = new double[n][n];
						for (var i = 0; i < n; i++)
							for (var j = 0; j < n; j++)
								matrix[i][j] = distances.get(i + 1, j + 1);
						outcomes.add(new Outcome(method, null, QuartetDisplay.ofDistances(matrix), 0, 0));
					}
					case "NeighborNet" -> {
						var splits = new SplitsBlock();
						new NeighborNet().compute(new ProgressSilent(), taxa, distances, splits);
						var seconds = (System.currentTimeMillis() - start) / 1000.0;
						var graph = new PhyloSplitsGraph();
						var used = new BitSet();
						EqualAngle.apply(new ProgressSilent(), true, taxa, splits, graph, new BitSet(), used);
						if (used.cardinality() < splits.getNsplits())
							ConvexHull.apply(new ProgressSilent(), taxa, splits, graph, used);
						var sides = new ArrayList<BitSet>();
						var weights = new double[splits.getNsplits()];
						for (var split : splits.getSplits()) {
							var side = new BitSet();
							for (var t = split.getA().nextSetBit(0); t >= 0; t = split.getA().nextSetBit(t + 1))
								side.set(t - 1);
							weights[sides.size()] = split.getWeight();
							sides.add(side);
						}
						outcomes.add(new Outcome(method, TaxonGraph.of(graph, taxa::getLabel, e -> splits.get(graph.getSplit(e)).getWeight()),
								QuartetDisplay.ofSplits(sides, weights, n), splits.getNsplits(), seconds));
					}
					case "RazorNet", "RazorNet+Simplify" -> {
						if (razorNetwork == null) {
							var algorithm = new razornet.algorithms.RazorNet();
							algorithm.optionRandomRestartsProperty().set(restarts);
							algorithm.optionRandomSeedProperty().set(seed);
							algorithm.optionMaxSecondsProperty().set(maxSeconds);
							if (!strategy.isEmpty())
								algorithm.optionStrategyProperty().set(razornet.algorithms.RazorNet.Strategy.valueOf(strategy));
							for (var entry : settings.entrySet()) {
								var option = splitstree6.options.Option.getName2Options(algorithm).get(entry.getKey());
								if (option == null)
									throw new IOException("Unknown option of RazorNet: " + entry.getKey());
								((javafx.beans.property.Property<Object>) option.getProperty()).setValue(
										option.getOptionValueType() == splitstree6.options.OptionValueType.Enum
												? option.getEnumValueForName(entry.getValue())
												: splitstree6.options.OptionValueType.parseType(option.getOptionValueType(), entry.getValue()));
							}
							razorNetwork = new NetworkBlock();
							algorithm.compute(new ProgressSilent(), taxa, distances, razorNetwork);
							razorSeconds = (System.currentTimeMillis() - start) / 1000.0;
							start = System.currentTimeMillis();
						}
						var network = razorNetwork;
						var seconds = razorSeconds;
						if (method.endsWith("+Simplify")) {
							var input0 = razorNetwork;
							network = TimeLimit.call(() -> {
								var filter = new SimplifyFilter();
								filter.optionReportTableProperty().set(false);
								var result = new NetworkBlock();
								filter.compute(new ProgressSilent(), taxa, input0, result);
								return result;
							}, timeLimit);
							seconds += (System.currentTimeMillis() - start) / 1000.0;
							if (network == null) { // the filter exceeded the time limit
								outcomes.add(new Outcome(method + tag, null, null, 0, seconds));
								continue;
							}
						}
						var graph = network.getGraph();
						var taxonGraph = TaxonGraph.of(graph, taxa::getLabel, graph::getWeight);
						outcomes.add(new Outcome(method + tag, taxonGraph, QuartetDisplay.ofGraph(taxonGraph, labels), 0, seconds));
						if (!networksDirectory.isEmpty()) {
							var file = new File(networksDirectory, FileUtils.replaceFileSuffix(new File(input).getName(), "-" + method.replace('+', '-') + ".nex"));
							try (var w = FileUtils.getOutputWriterPossiblyZIPorGZIP(file.getPath())) {
								ExportManager.getInstance().write(taxa, network, "NexusWithTaxa", w);
							}
						}
					}
					default -> throw new IOException("Unknown method: " + method);
				}
			}

			for (var outcome : outcomes) {
				if (outcome.display() == null) {
					out.printf("%s\t%d\t%s\ttimeout\t%.2f%n", input, n, outcome.method(), outcome.seconds());
					continue;
				}
				if (outcome.graph() == null) {
					out.printf("%s\t%d\t%s\tok\t\t\t\t\t\t\t\t\t", input, n, outcome.method());
					if (truth != null)
						out.print("\t" + QuartetDisplay.compare(truth, outcome.display(), use, maxQuartets, seed).toTsv());
					out.println();
					continue;
				}
				var shown = outcome.graph().taxonDistances(labels);
				double squares = 0, total = 0, minRatio = Double.POSITIVE_INFINITY, maxRatio = 0;
				for (var i = 0; i < n; i++) {
					for (var j = i + 1; j < n; j++) {
						var d = distances.get(i + 1, j + 1);
						squares += (shown[i][j] - d) * (shown[i][j] - d);
						total += d * d;
						if (d > 0) {
							minRatio = Math.min(minRatio, shown[i][j] / d);
							maxRatio = Math.max(maxRatio, shown[i][j] / d);
						}
					}
				}
				out.printf("%s\t%d\t%s\tok\t%.2f\t%d\t%d\t%d\t%s\t%d\t%.3f\t%.4f\t%.4f", input, n, outcome.method(), outcome.seconds(),
						outcome.graph().getNumberOfNodes(), outcome.graph().getNumberOfEdges(), outcome.graph().cycles(),
						jloda.util.StringUtils.trim("%.6f", outcome.graph().length()), outcome.splits(),
						100.0 * (1.0 - squares / total), minRatio, maxRatio);
				if (truth != null)
					out.print("\t" + QuartetDisplay.compare(truth, outcome.display(), use, maxQuartets, seed).toTsv());
				out.println();
				out.flush();
			}
		}
		if (out != System.out)
			out.close();
		System.exit(0); // the thread pool would keep the process alive
	}
}
