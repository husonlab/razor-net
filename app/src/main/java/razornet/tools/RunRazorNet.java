/*
 * RunRazorNet.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.tools;

import javafx.beans.property.Property;
import jloda.fx.util.ArgsOptions;
import jloda.fx.util.ProgramProperties;
import jloda.graph.Node;
import jloda.phylo.PhyloGraph;
import jloda.util.*;
import jloda.util.progress.ProgressListener;
import jloda.util.progress.ProgressPercentage;
import jloda.util.progress.ProgressSilent;
import razornet.algorithms.TCS;
import razornet.algorithms.utils.AnalysisLog;
import razornet.algorithms.utils.ExtensionTrace;
import splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances;
import splitstree6.algorithms.network.network2network.StretchFilter;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.ImportManager;
import splitstree6.io.writers.ExportManager;
import splitstree6.options.Option;
import splitstree6.options.OptionValueType;
import splitstree6.workflow.Algorithm;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Runs RazorNet or TCS on one or more distance files, without a GUI, and writes the
 * resulting networks and an optional report. Naming several algorithms runs each of them on every input. The
 * options of the algorithms are obtained from their {@code listOptions()} and offered as command-line options.
 * {@code --maxStretch} applies SplitsTree's stretch filter to each network, which then no longer realizes its
 * input. The exit status is 0 if every network was computed and realizes its input, 1 if not, and 2 on a usage
 * error.
 *
 * @author Daniel Huson, 2026
 */
public class RunRazorNet {
	/** short names for the options that are used most; the other options are named after the option */
	private static final Map<String, String[]> OPTION_ALIASES = Map.of(
			"RandomRestarts", new String[]{"-rs", "restarts"},
			"RandomSeed", new String[]{"-sd", "seed"},
			"MaxSeconds", new String[]{"-t", "maxSeconds"});

	/** Absolute tolerance when checking that the network realizes the input distances; as in {@code RunTestSuite}. */
	private static final double REALIZE_TOLERANCE = 1e-3;

	/** Exporter name for a bare NETWORK block, i.e. Nexus output without the header and the TAXA block. */
	private static final String NEXUS_BLOCK_ONLY = "NexusBlockOnly";

	/** What one engine did on one input; a row of the report. */
	private record Result(String file, String algorithm, int ntax, double length, int nodes, int edges,
						  int steinerNodes, double maxError, double seconds, int seed, int restarts,
						  boolean budgetHit, String error) {
		boolean failed() {
			return error != null;
		}
	}

	public static void main(String[] args) {
		try {
			razornet.algorithms.utils.AdvancedOptions.enable(); // the command line offers every option, basic and advanced
			ProgramProperties.setProgramName("RunRazorNet");
			ProgramProperties.setProgramVersion(razornet.main.Version.SHORT_DESCRIPTION);
			PeakMemoryUsageMonitor.start();
			var failures = (new RunRazorNet()).run(args);
			PeakMemoryUsageMonitor.report();
			System.exit(failures == 0 ? 0 : 1);
		} catch (UsageException ex) {
			if (!ex.getMessage().startsWith("Help"))
				System.err.println("Usage error: " + ex.getMessage());
			System.exit(2);
		} catch (Throwable ex) {
			// Throwable, so that an Error also ends the process
			System.err.println("Fatal error: " + ex);
			System.exit(2);
		}
	}

	/**
	 * Parses the arguments, runs every (input, engine) combination, and returns the number of failures.
	 */
	private int run(String[] args) throws Exception {
		final var options = new ArgsOptions(args, RunRazorNet.class,
				"Computes RazorNet networks from distances, without a GUI");
		options.setVersion(ProgramProperties.getProgramVersion());
		options.setLicense("This is free software, licensed under the terms of the GNU General Public License, Version 3.");
		options.setAuthors("Daniel H. Huson and Momoko Hayamizu");

		options.comment("Input Output:");
		final var inputFiles = options.getOptionMandatory("-i", "input", "File(s) containing a distance matrix (or directory)", new String[0]);
		final var outputFiles = options.getOption("-o", "output", "Output file(s), a directory, stdout, or .gz", new String[]{"stdout"});
		final var inputFormat = options.getOption("-f", "format", "Input format", ImportManager.getInstance().getAllFileFormats(), ImportManager.UNKNOWN_FORMAT);
		final var exportFormat = options.getOption("-e", "exporter", "Output format (Nexus writes a complete file, with the #nexus header and a TAXA block)",
				CollectionUtils.concatenate(ExportManager.getInstance().getExporterNames(new NetworkBlock()),
						List.of(NEXUS_BLOCK_ONLY)), "Nexus");
		final var reportFile = options.getOption("-rp", "report", "Write a per-instance report (TSV) to this file", "");

		options.comment("Algorithm:");
		// legal values are checked in createAlgorithms rather than by ArgsOptions, so that both
		// '-a RazorNet TCS' and '-a RazorNet,TCS' are accepted
		final var algorithmNames = expandCommas(options.getOption("-a", "algorithm",
				"Algorithm(s) to run: RazorNet, TCS; naming both runs and compares them (RazorNet's strategy: -strategy)",
				new String[]{"RazorNet"}));
		final var algorithms = createAlgorithms(algorithmNames);
		setupAlgorithmOptions(options, algorithms);

		options.comment(ArgsOptions.OTHER);
		final var inputExtension = options.getOption("-x", "inputExt", "File extension of input files (when input is a directory)", "");
		final var recursive = options.getOption("-r", "recursive", "Recursively visit all sub-directories (when input is a directory)", false);
		final var check = options.getOption("-c", "check", "Check that each network realizes its input distances", true);
		final var maxStretch = options.getOption("-ms", "maxStretch", "Post-processing with the stretch filter: remove edges while every pairwise distance stays within this percent of its value (0 = off)", 0.0);
		ProgramExecutorService.setNumberOfCoresToUse(options.getOption("-th", "threads", "Maximum number of threads to use (0=all available)", 0));
		final var propertiesFile = options.getOption("-p", "propertiesFile", "Properties file, to run with the GUI's saved option defaults", "");
		final var silent = options.getOption("-s", "silent", "Silent mode (hide all stderr output)", false);
		final var writeLog = options.getOption("-log", "log", "Write a detailed analysis log to <inputfile>.log for each input", false);
		options.done();

		// not loaded by default, so that the defaults of a run do not depend on the options last set in the GUI
		if (!propertiesFile.isBlank()) {
			ProgramProperties.load(propertiesFile);
			algorithms.clear();
			algorithms.addAll(createAlgorithms(algorithmNames)); // re-create so the loaded defaults take effect
			setupAlgorithmOptions(options, algorithms);
		}
		if (silent)
			Basic.hideSystemErr();

		final var inputs = collectInputFiles(inputFiles, inputExtension, recursive);
		final var outputs = resolveOutputFiles(inputs, outputFiles, exportFormat, algorithms.size() > 1);

		System.err.printf("Input files: %,d%n", inputs.size());
		System.err.printf("Algorithms:  %s%n", StringUtils.toString(algorithmNames, ", "));

		final var results = new ArrayList<Result>();
		for (var i = 0; i < inputs.size(); i++) {
			final var inputFile = inputs.get(i);
			for (var algorithm : algorithms) {
				// one output per input, or a single output for all inputs
				var output = outputs.get(outputs.size() > 1 ? i : 0);
				if (algorithms.size() > 1 && !output.equals("stdout"))
					output = FileUtils.replaceFileSuffix(output, "-" + algorithm.getClass().getSimpleName() + FileUtils.getFileSuffix(output));
				results.add(runOne(inputFile, output, inputFormat, exportFormat, algorithm, maxStretch, check, silent, writeLog));
			}
		}

		if (!reportFile.isBlank())
			writeReport(reportFile, results);

		var failures = (int) results.stream().filter(Result::failed).count();
		System.err.printf("%nComputed: %,d, failed: %,d%n", results.size() - failures, failures);
		for (var result : results) {
			if (result.failed())
				System.err.printf("  - %s (%s): %s%n", result.file(), result.algorithm(), result.error());
		}
		return failures;
	}

	/**
	 * Runs one algorithm on one input file and writes its network. Any exception is caught and recorded in the
	 * result, so that a batch continues.
	 */
	private Result runOne(String inputFile, String outputFile, String inputFormat, String exportFormat,
						  Algorithm<DistancesBlock, NetworkBlock> algorithm, double maxStretch, boolean check, boolean silent, boolean writeLog) {
		// the class name, as getName() is the display name
		final var name = algorithm.getClass().getSimpleName();
		final var seed = intOption(algorithm, "RandomSeed");
		final var restarts = intOption(algorithm, "RandomRestarts");
		final var maxSeconds = intOption(algorithm, "MaxSeconds");

		var taxaBlock = new TaxaBlock();
		var distancesBlock = new DistancesBlock();
		var networkBlock = new NetworkBlock();
		var start = System.currentTimeMillis();
		try {
			read(inputFile, inputFormat, taxaBlock, distancesBlock);
			System.err.printf("%s on %s (%d taxa)...%n", name, FileUtils.getFileNameWithoutPath(inputFile), taxaBlock.getNtax());

			// RazorNet writes its analysis log itself (optionLogFile); for the other algorithms the log is written here
			var ownLogging = writeLog && algorithm instanceof razornet.algorithms.RazorNet r && setLogFile(r, inputFile + ".log");
			if (writeLog && !ownLogging)
				ExtensionTrace.start();
			try (ProgressListener progress = (silent ? new ProgressSilent() : new ProgressPercentage())) {
				algorithm.compute(progress, taxaBlock, distancesBlock, networkBlock);
				if (maxStretch > 0) {
					var filter = new StretchFilter();
					filter.setOptionMaxStretchPercent(maxStretch);
					var filtered = new NetworkBlock();
					filter.compute(progress, taxaBlock, networkBlock, filtered);
					networkBlock = filtered;
					System.err.printf("  stretch filter: %s%n", filter.getShortDescription());
				}
			}
			var seconds = (System.currentTimeMillis() - start) / 1000.0;

			var graph = networkBlock.getGraph();
			// the summary line that the program shows under the drawing
			if (!silent)
				System.err.printf("  %s%n", splitstree6.view.network.AdditionalConsoleOutput.infoFor(networkBlock, distancesBlock, null));
			var maxError = check ? maxRealizationError(graph, distancesBlock) : Double.NaN;
			var length = graph.edgeStream().mapToDouble(graph::getWeight).sum();
			var steiner = (int) graph.nodeStream().filter(v -> !graph.hasTaxa(v)).count();

			write(outputFile, exportFormat, taxaBlock, networkBlock);

			if (writeLog && !ownLogging) {
				var n = distancesBlock.getNtax();
				var raw = new double[n][n];
				for (var a = 0; a < n; a++)
					for (var b = 0; b < n; b++) raw[a][b] = distancesBlock.get(a + 1, b + 1);
				ExtensionTrace.stop();
				AnalysisLog.write(inputFile + ".log", name, optionsSummary(algorithm), taxaBlock, raw,
						intOption(algorithm, "SignificantDigits"), graph, seconds, maxError, true, "");
			}

			// after the stretch filter the network is not exact, so the realization check does not apply
			var error = (check && maxStretch == 0 && maxError > REALIZE_TOLERANCE)
					? "does not realize its distances (max error %.3g)".formatted(maxError) : null;

			return new Result(inputFile, name, taxaBlock.getNtax(), length, graph.getNumberOfNodes(),
					graph.getNumberOfEdges(), steiner, maxError, seconds, seed, restarts,
					maxSeconds > 0 && seconds >= maxSeconds, error);
		} catch (Exception ex) {
			var seconds = (System.currentTimeMillis() - start) / 1000.0;
			var message = ex.getClass().getSimpleName() + (ex.getMessage() != null ? ": " + ex.getMessage() : "");
			return new Result(inputFile, name, taxaBlock.getNtax(), 0, 0, 0, 0, Double.NaN, seconds,
					seed, restarts, false, message);
		}
	}

	/** Points RazorNet at its own log file; returns true so the caller knows not to write one itself. */
	private static boolean setLogFile(razornet.algorithms.RazorNet algorithm, String logFile) {
		algorithm.setOptionLogFile(logFile);
		return true;
	}

	/** All {@code listOptions()} values of the algorithm, so a log says exactly what produced it. */
	private static String optionsSummary(Algorithm<DistancesBlock, NetworkBlock> algorithm) {
		var name2option = Option.getName2Options(algorithm);
		var buf = new StringBuilder();
		for (var name : algorithm.listOptions()) {
			var option = name2option.get(name.replaceAll("^option", ""));
			if (option != null)
				buf.append(buf.isEmpty() ? "" : ", ").append(name.replaceAll("^option", ""))
						.append("=").append(option.getProperty().getValue());
		}
		return buf.toString();
	}

	/**
	 * Reads a distance matrix, using the reader that {@link ImportManager} determines from the file (or the
	 * one named by {@code --format}).
	 */
	private void read(String inputFile, String inputFormat, TaxaBlock taxaBlock, DistancesBlock distancesBlock) throws IOException {
		var reader = (inputFormat.equals(ImportManager.UNKNOWN_FORMAT) ? ImportManager.getInstance().getReader(inputFile)
				: ImportManager.getInstance().getImporterByDataTypeAndFileFormat(DistancesBlock.class, inputFormat));
		if (reader == null)
			throw new IOException("Can't determine input format: " + inputFile);
		if (!DistancesBlock.class.isAssignableFrom(reader.getToClass()))
			throw new IOException("Not distances (%s): %s".formatted(reader.getToClass().getSimpleName(), inputFile));
		((splitstree6.io.utils.DataReaderBase<DistancesBlock>) reader).read(new ProgressSilent(), inputFile, taxaBlock, distancesBlock);
	}

	/**
	 * Writes the network. Nexus output is a complete file, with the "#nexus" header and a TAXA block; a bare NETWORK
	 * block is written with {@code -e NexusBlockOnly}.
	 */
	private void write(String outputFile, String exportFormat, TaxaBlock taxaBlock, NetworkBlock networkBlock) throws IOException {
		var exporter = switch (exportFormat) {
			case "Nexus" -> "NexusWithTaxa";
			case NEXUS_BLOCK_ONLY -> "Nexus";
			default -> exportFormat;
		};
		try (Writer w = FileUtils.getOutputWriterPossiblyZIPorGZIP(outputFile)) {
			ExportManager.getInstance().write(taxaBlock, networkBlock, exporter, w);
		}
	}

	/**
	 * Largest absolute difference between an input distance and the corresponding shortest-path distance in
	 * the network. This is the fundamental correctness property: the network must realize the metric.
	 */
	private static double maxRealizationError(PhyloGraph graph, DistancesBlock distancesBlock) {
		var nodes = new ArrayList<Node>();
		var taxonOf = new ArrayList<Integer>();
		graph.nodeStream().filter(graph::hasTaxa).forEach(v -> {
			nodes.add(v);
			taxonOf.add(graph.getTaxon(v));
		});
		if (nodes.size() != distancesBlock.getNtax())
			return Double.POSITIVE_INFINITY; // wrong number of labeled nodes => not a valid realization
		var maxError = 0.0;
		for (var a = 0; a < nodes.size(); a++) {
			// one shortest-path computation per taxon
			var distanceFromA = CheckPairwiseDistances.graphDistances(graph, nodes.get(a));
			for (var b = a + 1; b < nodes.size(); b++) {
				var in = distancesBlock.get(taxonOf.get(a), taxonOf.get(b));
				double out = distanceFromA.get(nodes.get(b));
				maxError = Math.max(maxError, Math.abs(in - out));
			}
		}
		return maxError;
	}

	/** Splits any comma-separated entries, so that a list can be written either way. */
	private static String[] expandCommas(String[] values) {
		var result = new ArrayList<String>();
		for (var value : values) {
			for (var part : value.split(",")) {
				if (!part.isBlank())
					result.add(part.trim());
			}
		}
		return result.toArray(new String[0]);
	}

	/** Instantiates the named engines. All are {@code Distances2Network} algorithms and are driven identically. */
	private static ArrayList<Algorithm<DistancesBlock, NetworkBlock>> createAlgorithms(String[] names) throws UsageException {
		var algorithms = new ArrayList<Algorithm<DistancesBlock, NetworkBlock>>();
		for (var name : names) {
			switch (name) {
				case "RazorNet" -> algorithms.add(new razornet.algorithms.RazorNet());
				case "TCS" -> algorithms.add(new TCS());
				default -> throw new UsageException("Unknown algorithm: " + name);
			}
		}
		return algorithms;
	}

	/**
	 * Declares one command-line option for each option of the algorithms and applies the values given. With several
	 * algorithms, the union of their options is offered, and a value is applied to every algorithm that has the
	 * option.
	 */
	private static void setupAlgorithmOptions(ArgsOptions options, List<Algorithm<DistancesBlock, NetworkBlock>> algorithms) throws UsageException {
		final var name2options = new LinkedHashMap<String, List<Option>>();
		for (var algorithm : algorithms) {
			for (var option : Option.getAllOptions(algorithm)) {
				name2options.computeIfAbsent(option.getName(), k -> new ArrayList<>()).add(option);
			}
		}

		for (var entry : name2options.entrySet()) {
			final var name = entry.getKey();
			final var alias = OPTION_ALIASES.get(name);
			final var longKey = (alias != null ? alias[1] : name.substring(0, 1).toLowerCase() + name.substring(1));
			// a generated option has no short form
			final var shortKey = (alias != null ? alias[0] : "-" + longKey);
			final var first = entry.getValue().get(0);
			// an option without a tool tip: use its name as the description
			var description = first.getToolTipText();
			if (description == null || description.equalsIgnoreCase("option" + name))
				description = StringUtils.fromCamelCase(name);

			final Object value;
			switch (first.getOptionValueType()) {
				case Integer -> value = options.getOption(shortKey, longKey, description, (Integer) first.getProperty().getValue());
				case Double -> value = options.getOption(shortKey, longKey, description, (Double) first.getProperty().getValue());
				case Boolean -> value = options.getOption(shortKey, longKey, description, (Boolean) first.getProperty().getValue());
				case Enum -> {
					var chosen = options.getOption(shortKey, longKey, description, first.getLegalValues(),
							first.getProperty().getValue().toString());
					value = first.getEnumValueForName(chosen);
				}
				// a String is passed as it is: the round trip below would wrap it in quotes
				case String -> value = options.getOption(shortKey, longKey, description, (String) first.getProperty().getValue());
				default -> value = OptionValueType.parseType(first.getOptionValueType(),
						options.getOption(shortKey, longKey, description,
								OptionValueType.toStringType(first.getOptionValueType(), first.getProperty().getValue())));
			}
			if (value != null) {
				for (var option : entry.getValue()) {
					((Property<Object>) option.getProperty()).setValue(value);
				}
			}
		}
	}

	/**
	 * Expands the input arguments: files as given, directories by their contents (see {@code -x} and {@code -r}). A
	 * file found in a directory is skipped if it does not contain distances.
	 */
	private static ArrayList<String> collectInputFiles(String[] inputFiles, String extension, boolean recursive) throws IOException {
		var files = new ArrayList<String>();
		var skipped = 0;
		for (var input : inputFiles) {
			if (FileUtils.isDirectory(input)) {
				var found = FileUtils.getAllFilesInDirectory(input, recursive,
						extension.isBlank() ? new String[0] : new String[]{extension});
				for (var file : found.stream().sorted().toList()) {
					if (DistancesBlock.class.equals(ImportManager.getInstance().determineInputType(file)))
						files.add(file);
					else
						skipped++;
				}
			} else {
				FileUtils.checkFileReadableNonEmpty(input);
				files.add(input);
			}
		}
		if (skipped > 0)
			System.err.printf("Skipped %,d file(s) that do not contain distances%n", skipped);
		if (files.isEmpty())
			throw new IOException("No input files found");
		return files;
	}

	/**
	 * Determines the output for each input: one output per input, a directory, or stdout. With several algorithms,
	 * the output has to be a directory or stdout.
	 */
	private static ArrayList<String> resolveOutputFiles(List<String> inputs, String[] outputFiles, String exportFormat,
													    boolean severalAlgorithms) throws IOException {
		var outputs = new ArrayList<String>();
		if (outputFiles.length == 1 && FileUtils.isDirectory(outputFiles[0])) {
			var exporter = ExportManager.getInstance().getExporterByName(NetworkBlock.class, exportFormat);
			var extension = "." + (exporter != null && !exporter.getFileExtensions().isEmpty() ? exporter.getFileExtensions().get(0) : "txt");
			for (var input : inputs) {
				// note: not FileUtils.getFilePath, which resolves against a reference *file* and so would
				// write into the parent of the requested directory
				outputs.add(new File(outputFiles[0],
						FileUtils.replaceFileSuffix(FileUtils.getFileNameWithoutPath(input), "-out" + extension)).getPath());
			}
		} else if (outputFiles.length == 1) {
			if (severalAlgorithms && !outputFiles[0].equals("stdout"))
				throw new IOException("Several algorithms: --output must be a directory or stdout, got: " + outputFiles[0]);
			outputs.add(outputFiles[0]);
		} else if (outputFiles.length == inputs.size()) {
			outputs.addAll(List.of(outputFiles));
		} else {
			throw new IOException("Number of output files (%d) must be 1 or match the number of input files (%d)"
					.formatted(outputFiles.length, inputs.size()));
		}
		return outputs;
	}

	/** writes the report: one row per input and algorithm */
	private static void writeReport(String reportFile, List<Result> results) throws IOException {
		try (var w = FileUtils.getOutputWriterPossiblyZIPorGZIP(reportFile)) {
			w.write("file\talgorithm\tntax\tlength\tnodes\tedges\tsteinerNodes\tmaxRealizationError\tseconds\tseed\trestarts\tbudgetHit\terror\n");
			for (var r : results) {
				w.write("%s\t%s\t%d\t%s\t%d\t%d\t%d\t%s\t%.2f\t%d\t%d\t%s\t%s%n".formatted(
						r.file(), r.algorithm(), r.ntax(), StringUtils.removeTrailingZerosAfterDot("%.4f", r.length()),
						r.nodes(), r.edges(), r.steinerNodes(),
						Double.isNaN(r.maxError()) ? "NA" : "%.3g".formatted(r.maxError()),
						r.seconds(), r.seed(), r.restarts(), r.budgetHit(), r.error() == null ? "" : r.error()));
			}
		}
		System.err.println("Report: " + reportFile);
	}

	private static int intOption(Algorithm<?, ?> algorithm, String name) {
		var option = Option.getName2Options(algorithm).get(name);
		return option == null ? 0 : ((Number) option.getProperty().getValue()).intValue();
	}
}
