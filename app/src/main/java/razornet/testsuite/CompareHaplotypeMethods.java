/*
 * CompareHaplotypeMethods.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.Basic;
import jloda.util.FileUtils;
import jloda.util.ProgramExecutorService;
import jloda.util.progress.ProgressSilent;
import splitstree6.algorithms.characters.characters2network.Characters2Network;
import splitstree6.algorithms.characters.characters2network.MedianJoining;
import splitstree6.algorithms.network.network2network.SimplifyFilter;
import splitstree6.data.CharactersBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.ImportManager;
import splitstree6.io.writers.ExportManager;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * Compares haplotype network methods on alignments, without a GUI: median joining (SplitsTree), TCS, RazorHaploNet
 * and RazorHaploNet followed by the Simplify filter.
 * <p>
 * All methods get the same input: the alignment without the columns that contain a gap, a missing or an ambiguous
 * character, and with identical sequences merged into one haplotype. For each input and method the tool reports
 * one row: the size of the network, the number of mutations it draws, the excess of the mutations along its
 * shortest paths over the sequence differences, and the time. If the genealogy that generated the sequences is
 * given, it also reports how many of the edges of the true haplotype tree the network contains
 * ({@link HaplotypeTruth}) and how the quartet supports it displays compare with the true ones
 * ({@link QuartetDisplay}). The method {@code Distances} reports the quartet supports that the sequence
 * differences give by themselves, and the method {@code Truth} measures the true haplotype tree in the same way
 * as a network.
 *
 * <pre>
 * CompareHaplotypeMethods [--methods Distances,Truth,MedianJoining,TCS,RazorHaploNet,RazorHaploNet+Simplify]
 *                         [--truth treeName sequencesName] [--epsilon N] [--strategy S] [--restarts N]
 *                         [--maxSeconds N] [--timeLimit N] [--seed N] [--threads N] [--maxQuartets N]
 *                         [--report file] [--networks directory] [--silent] alignment...
 * </pre>
 * {@code --maxSeconds} is the time budget of RazorHaploNet and TCS; {@code --timeLimit} ends median joining or the
 * Simplify filter after that many seconds, which is then reported as a timeout.
 * {@code --truth} names two files next to each input: the genealogy in Newick format with all nodes named, and the
 * sequences of all its nodes in FASTA format.
 *
 * @author Daniel Huson, 2026
 */
public class CompareHaplotypeMethods {
	private static final String HEADER = "file\tsequences\thaplotypes\tsites\tmethod\tstatus\tseconds\tnodes\tauxiliary\tedges\tcycles"
										 + "\tmutations\texcess\texactPairs\tmaxRatio";
	private static final String TRUTH_HEADER = "trueMutations\ttrueEdges\tshownEdges\tsharedEdges\tunplaced";

	public static void main(String[] args) throws Exception {
		razornet.algorithms.utils.AdvancedOptions.enable();
		var methods = List.of("MedianJoining", "TCS", "RazorHaploNet", "RazorHaploNet+Simplify");
		var truthTree = "";
		var truthSequences = "";
		var epsilon = 0;
		var strategy = "";
		var restarts = 1000;
		var maxSeconds = 60;
		var seed = 0;
		var threads = 0;
		var maxQuartets = 100000L;
		var timeLimit = 0;
		var reportFile = "";
		var networksDirectory = "";
		var silent = false;
		var inputs = new ArrayList<String>();
		for (var i = 0; i < args.length; i++) {
			switch (args[i]) {
				case "--methods" -> methods = List.of(args[++i].split(","));
				case "--timeLimit" -> timeLimit = Integer.parseInt(args[++i]);
				case "--truth" -> {
					truthTree = args[++i];
					truthSequences = args[++i];
				}
				case "--epsilon" -> epsilon = Integer.parseInt(args[++i]);
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
			System.err.println("Usage: CompareHaplotypeMethods [--methods Distances,Truth,MedianJoining,TCS,RazorHaploNet,RazorHaploNet+Simplify] [--truth treeName sequencesName] [--epsilon N] [--strategy S] [--restarts N] [--maxSeconds N] [--timeLimit N] [--seed N] [--threads N] [--maxQuartets N] [--report file] [--networks directory] [--silent] alignment...");
			System.exit(2);
		}
		ProgramExecutorService.setNumberOfCoresToUse(threads);
		if (!networksDirectory.isEmpty())
			Files.createDirectories(new File(networksDirectory).toPath());
		var out = reportFile.isEmpty() ? System.out : new PrintStream(reportFile);
		if (silent)
			Basic.hideSystemErr();
		var hasTruth = !truthTree.isEmpty();
		out.println(HEADER + (hasTruth ? "\t" + TRUTH_HEADER + "\t" + QuartetDisplay.Scores.HEADER : ""));

		for (var input : inputs) {
			var base = FileUtils.replaceFileSuffix(new File(input).getName(), "");
			// the common input: complete deletion, identical sequences merged
			var allTaxa = new TaxaBlock();
			var allCharacters = new CharactersBlock();
			read(input, allTaxa, allCharacters);
			// a gap, a missing character, N or an IUPAC ambiguity code; any other symbol is a state (so are 0 and 1)
			var incomplete = "-?NRYSWKMBDHV" + Character.toUpperCase(allCharacters.getGapCharacter())
							 + Character.toUpperCase(allCharacters.getMissingCharacter());
			var kept = new ArrayList<Integer>();
			for (var j = 1; j <= allCharacters.getNchar(); j++) {
				var complete = true;
				for (var t = 1; t <= allTaxa.getNtax() && complete; t++)
					complete = incomplete.indexOf(Character.toUpperCase(allCharacters.get(t, j))) < 0;
				if (complete)
					kept.add(j);
			}
			var haplotypes = new LinkedHashMap<String, String>(); // label -> sequence
			var sequence2label = new LinkedHashMap<String, String>();
			var members = new LinkedHashMap<String, List<String>>();
			for (var t = 1; t <= allTaxa.getNtax(); t++) {
				var buf = new StringBuilder();
				for (var j : kept)
					buf.append(Character.toUpperCase(allCharacters.get(t, j)));
				var label = sequence2label.get(buf.toString());
				if (label == null) {
					label = allTaxa.getLabel(t).replaceAll("[^A-Za-z0-9_.-]", "_");
					while (haplotypes.containsKey(label))
						label += "_";
					sequence2label.put(buf.toString(), label);
					haplotypes.put(label, buf.toString());
					members.put(label, new ArrayList<>());
				}
				members.get(label).add(allTaxa.getLabel(t));
			}
			var workDirectory = networksDirectory.isEmpty() ? Files.createTempDirectory("haplotypes").toFile() : new File(networksDirectory);
			var haplotypeFile = new File(workDirectory, base + ".haplotypes.fasta");
			try (var w = new PrintWriter(haplotypeFile)) {
				for (var entry : haplotypes.entrySet())
					w.printf(">%s%n%s%n", entry.getKey(), entry.getValue());
			}
			try (var w = new PrintWriter(new File(workDirectory, base + ".haplotypes.tsv"))) {
				for (var entry : members.entrySet())
					w.printf("%s\t%d\t%s%n", entry.getKey(), entry.getValue().size(), String.join(",", entry.getValue()));
			}
			var taxa = new TaxaBlock();
			var characters = new CharactersBlock();
			read(haplotypeFile.getPath(), taxa, characters);
			var labels = new ArrayList<>(haplotypes.keySet());
			var n = labels.size();
			if (taxa.getNtax() != n)
				throw new IOException("Haplotype file read back with %d instead of %d sequences".formatted(taxa.getNtax(), n));
			var hamming = new double[n][n];
			for (var i = 0; i < n; i++)
				for (var j = 0; j < n; j++)
					hamming[i][j] = TaxonGraph.differences(haplotypes.get(labels.get(i)), haplotypes.get(labels.get(j)));

			HaplotypeTruth truth = null;
			QuartetDisplay.Support trueSupport = null;
			var use = new boolean[n];
			java.util.Arrays.fill(use, true);
			if (hasTruth) {
				var directory = new File(input).getParentFile();
				var keptSites = kept;
				truth = new HaplotypeTruth(Files.readString(new File(directory, truthTree).toPath()),
						SimulateTreeMixtures.readFasta(new File(directory, truthSequences)), full -> {
					var buf = new StringBuilder();
					for (var j : keptSites)
						buf.append(Character.toUpperCase(full.charAt(j - 1)));
					return buf.toString();
				}, haplotypes);
				trueSupport = QuartetDisplay.ofGraph(truth.getGraph(), labels);
				for (var i = 0; i < n; i++)
					use[i] = !truth.getUnplaced().contains(labels.get(i));
			}
			var prefix = "%s\t%d\t%d\t%d".formatted(input, allTaxa.getNtax(), n, kept.size());

			NetworkBlock razorNetwork = null;
			var razorSeconds = 0.0;
			for (var method : methods) {
				var start = System.currentTimeMillis();
				NetworkBlock network = null;
				TaxonGraph graph = null;
				double seconds;
				switch (method) {
					case "Distances" -> {
						out.printf("%s\t%s\tok\t\t\t\t\t\t\t\t\t", prefix, method);
						if (hasTruth)
							out.printf("\t%s\t%d\t\t\t%d\t%s", jloda.util.StringUtils.trim("%.1f", truth.getGraph().length()),
									truth.getEdges().size(), truth.getUnplaced().size(),
									QuartetDisplay.compare(trueSupport, QuartetDisplay.ofDistances(hamming), use, maxQuartets, seed).toTsv());
						out.println();
						continue;
					}
					case "Truth" -> {
						if (!hasTruth)
							continue;
						graph = truth.getGraph();
						seconds = 0;
					}
					case "MedianJoining" -> {
						var algorithm = new MedianJoining();
						algorithm.setOptionEpsilon(epsilon);
						network = TimeLimit.call(() -> compute(algorithm, taxa, characters), timeLimit);
						seconds = (System.currentTimeMillis() - start) / 1000.0;
					}
					case "TCS" -> {
						var algorithm = new razornet.algorithms.TCSHaplotypeNetwork();
						algorithm.optionMaxSecondsProperty().set(maxSeconds);
						network = TimeLimit.call(() -> compute(algorithm, taxa, characters), timeLimit);
						seconds = (System.currentTimeMillis() - start) / 1000.0;
					}
					case "RazorHaploNet", "RazorHaploNet+Simplify" -> {
						if (razorNetwork == null) {
							var algorithm = new razornet.algorithms.RazorHaploNet();
							algorithm.optionRandomRestartsProperty().set(restarts);
							algorithm.optionRandomSeedProperty().set(seed);
							algorithm.optionMaxSecondsProperty().set(maxSeconds);
							if (!strategy.isEmpty())
								algorithm.optionStrategyProperty().set(razornet.algorithms.RazorNet.Strategy.valueOf(strategy));
							razorNetwork = compute(algorithm, taxa, characters);
							razorSeconds = (System.currentTimeMillis() - start) / 1000.0;
							start = System.currentTimeMillis();
						}
						network = razorNetwork;
						seconds = razorSeconds;
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
						}
					}
					default -> throw new IOException("Unknown method: " + method);
				}
				if (network == null && graph == null) { // the method exceeded the time limit
					out.printf("%s\t%s\ttimeout\t%.2f%n", prefix, method, seconds);
					out.flush();
					continue;
				}
				if (network != null && !networksDirectory.isEmpty()) {
					try (var w = FileUtils.getOutputWriterPossiblyZIPorGZIP(new File(networksDirectory, base + "-" + method.replace('+', '-') + ".nex").getPath())) {
						ExportManager.getInstance().write(taxa, network, "NexusWithTaxa", w);
					}
				}

				if (graph == null) {
					var block = network;
					graph = TaxonGraph.ofHaplotypeNetwork(network.getGraph(), taxa::getLabel,
							v -> block.getNodeData(v).get(NetworkBlock.NODE_STATES_KEY));
				}
				// the true tree is measured on the haplotypes that it places
				var shown = graph.taxonDistances(labels);
				double excess = 0, maxRatio = 0;
				var exact = 0;
				var pairs = 0;
				for (var i = 0; i < n; i++) {
					for (var j = i + 1; j < n; j++) {
						if (Double.isNaN(shown[i][j]))
							continue;
						pairs++;
						excess += shown[i][j] - hamming[i][j];
						if (Math.abs(shown[i][j] - hamming[i][j]) < 1e-9)
							exact++;
						maxRatio = Math.max(maxRatio, shown[i][j] / hamming[i][j]);
					}
				}
				var auxiliary = 0;
				for (var v = 0; v < graph.getNumberOfNodes(); v++)
					if (!graph.hasTaxon(v))
						auxiliary++;
				out.printf("%s\t%s\tok\t%.2f\t%d\t%d\t%d\t%d\t%s\t%s\t%.4f\t%.4f", prefix, method, seconds, graph.getNumberOfNodes(), auxiliary,
						graph.getNumberOfEdges(), graph.cycles(), jloda.util.StringUtils.trim("%.1f", graph.length()),
						jloda.util.StringUtils.trim("%.1f", excess), pairs > 0 ? (double) exact / pairs : 1.0, maxRatio);
				if (hasTruth) {
					var edges = method.equals("Truth") ? truth.getEdges() : HaplotypeTruth.edgesBySequence(graph);
					var shared = edges.stream().filter(truth.getEdges()::contains).count();
					out.printf("\t%s\t%d\t%d\t%d\t%d\t%s", jloda.util.StringUtils.trim("%.1f", truth.getGraph().length()),
							truth.getEdges().size(), edges.size(), shared, truth.getUnplaced().size(),
							QuartetDisplay.compare(trueSupport, QuartetDisplay.ofGraph(graph, labels), use, maxQuartets, seed).toTsv());
				}
				out.println();
				out.flush();
			}
		}
		if (out != System.out)
			out.close();
		System.exit(0); // the thread pool would keep the process alive
	}

	private static void read(String file, TaxaBlock taxa, CharactersBlock characters) throws IOException {
		var reader = ImportManager.getInstance().getReader(file);
		if (reader == null || !CharactersBlock.class.isAssignableFrom(reader.getToClass()))
			throw new IOException("Not an alignment: " + file);
		((splitstree6.io.utils.DataReaderBase<CharactersBlock>) reader).read(new ProgressSilent(), file, taxa, characters);
	}

	private static NetworkBlock compute(Characters2Network algorithm, TaxaBlock taxa, CharactersBlock characters) throws IOException {
		var network = new NetworkBlock();
		algorithm.compute(new ProgressSilent(), taxa, characters, network);
		return network;
	}
}
