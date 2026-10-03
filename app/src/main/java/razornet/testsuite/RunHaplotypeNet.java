/*
 * RunHaplotypeNet.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.progress.ProgressPercentage;
import jloda.util.progress.ProgressSilent;
import splitstree6.data.CharactersBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.ImportManager;

import java.util.ArrayList;
import java.util.TreeMap;

/**
 * Runs {@code RazorHaploNet} on an alignment, without a GUI, and reports the network: its length and size, the
 * agreement of the mutations drawn on each edge with its weight, and the excess of the displayed distances over
 * the sequence differences. {@code --engine} selects RazorNet (the default) or TCS ({@code tcs});
 * {@code --strategy} selects the strategy of RazorNet.
 *
 * <pre>
 * RunHaplotypeNet alignment [--engine razornet|tcs] [--strategy Strict|Relaxed|Both] [--seed N]
 *                 [--maxSeconds N] [--restarts N] [--labeling MinExcess|Parsimony] [--raw|--noRemove|--noContract]
 *                 [--missingAsState] [--silent] [--set optionName=value]*
 * </pre>
 * {@code --raw} turns off the removal of superfluous edges and the contraction of unlabeled edges.
 *
 * @author Daniel Huson, 2026
 */
public class RunHaplotypeNet {
	/**
	 * Reports the excess (the mutations shown along the shortest paths minus the sequence differences, summed over
	 * all pairs of taxa) and the largest shortfall, which must be 0. The sequence differences are counted in the same
	 * way as by the algorithm: over the columns without missing or ambiguous data, or over all columns with
	 * {@code missingAsState}.
	 */
	private static void reportRealization(TaxaBlock taxa, CharactersBlock characters, NetworkBlock networkBlock,
										  boolean missingAsState) {
		var graph = networkBlock.getGraph();
		var gap = Character.toUpperCase(characters.getGapCharacter());
		var missing = Character.toUpperCase(characters.getMissingCharacter());
		var kept = new ArrayList<Integer>();
		for (var j = 1; j <= characters.getNchar(); j++) {
			var ok = true;
			if (!missingAsState)
				for (var i = 1; i <= taxa.getNtax() && ok; i++) {
					var c = Character.toUpperCase(characters.get(i, j));
					if (c != gap && (c == missing || c == 'N' || "RYSWKMBDHV".indexOf(c) >= 0)) ok = false;
				}
			if (ok) kept.add(j);
		}
		var nodes = new ArrayList<jloda.graph.Node>();
		var taxonOf = new ArrayList<Integer>();
		graph.nodeStream().filter(graph::hasTaxa).forEach(v -> {
			nodes.add(v);
			taxonOf.add(graph.getTaxon(v));
		});
		var shortfall = 0.0;
		var worst = "";
		var totalExcess = 0.0;
		var pairs = 0;
		for (var a = 0; a < nodes.size(); a++) {
			// one shortest-path computation per taxon
			var distanceFromA = splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances
					.graphDistances(graph, nodes.get(a));
			for (var b = a + 1; b < nodes.size(); b++) {
				var s = taxonOf.get(a);
				var t = taxonOf.get(b);
				var d = 0;
				for (var j : kept)
					if (Character.toUpperCase(characters.get(s, j)) != Character.toUpperCase(characters.get(t, j))) d++;
				double shown = distanceFromA.get(nodes.get(b));
				totalExcess += shown - d;
				pairs++;
				if (d - shown > shortfall) {
					shortfall = d - shown;
					worst = "%s vs %s: differ at %d sites, network shows %s".formatted(taxa.getLabel(s), taxa.getLabel(t), d, shown);
				}
			}
		}
		System.out.printf("  excess (mutations shown - sequence differences) over %,d pairs: %.1f (measured over %,d of %,d sites)%n",
				pairs, totalExcess, kept.size(), characters.getNchar());
		System.out.printf("  worst SHORTFALL (must be 0): %s%s%n",
				shortfall == Math.rint(shortfall) ? String.valueOf((long) shortfall) : String.valueOf(shortfall),
				shortfall > 1e-9 ? "   <-- BUG: " + worst : "");
	}

	/**
	 * sets an option by name; an algorithm that does not have the option ignores it
	 */
	private static void setOptionByName(splitstree6.options.IOptionsCarrier algorithm, String name, Object value) {
		var option = splitstree6.options.Option.getName2Options(algorithm).get(name);
		if (option == null)
			return;   // this engine has no such option, e.g. TCS has no restarts or seed
		var toSet = (option.getOptionValueType() == splitstree6.options.OptionValueType.Enum
				? option.getEnumValueForName(String.valueOf(value)) : value);
		((javafx.beans.property.Property<Object>) option.getProperty()).setValue(toSet);
	}

	public static void main(String[] args) throws Exception {
		razornet.algorithms.utils.AdvancedOptions.enable(); // --set reaches every option, basic and advanced
		if (args.length == 0 || args[0].startsWith("-")) {
			System.err.println("Usage: RunHaplotypeNet alignment [--engine razornet|tcs] [--strategy Strict|Relaxed|Both] [--seed N] [--maxSeconds N] [--restarts N] [--labeling M] [--raw|--noRemove|--noContract] [--missingAsState] [--silent] [--set optionName=value]*");
			System.exit(2);
		}
		var input = args[0];
		var maxSeconds = 60;
		var restarts = 100;
		var labeling = "MinExcess";
		var raw = false;
		var noRemove = false;
		var noContract = false;
		var missingAsState = false;
		var silent = false;
		var engine = "razornet";   // RazorNet; "tcs" TCS
		var strategy = "";         // RazorNet's strategy; empty = the algorithm's default (Both)
		var seed = 0;
		// any further option, as name=value, so it needs no flag of its own here
		var settings = new java.util.LinkedHashMap<String, String>();
		for (var i = 1; i < args.length; i++) {
			switch (args[i]) {
				case "--engine", "-e" -> engine = args[++i];
				case "--strategy" -> strategy = args[++i];
				case "--seed" -> seed = Integer.parseInt(args[++i]);
				case "--set" -> {
					var assignment = args[++i];
					var eq = assignment.indexOf('=');
					if (eq < 1) {
						System.err.println("Expected optionName=value, got: " + assignment);
						System.exit(2);
					}
					settings.put(assignment.substring(0, eq), assignment.substring(eq + 1));
				}
				case "--maxSeconds", "-t" -> maxSeconds = Integer.parseInt(args[++i]);
				case "--restarts", "-r" -> restarts = Integer.parseInt(args[++i]);
				case "--labeling", "-l" -> labeling = args[++i];
				case "--raw" -> raw = true;
				case "--noRemove" -> noRemove = true;
				case "--noContract" -> noContract = true;
				case "--missingAsState", "-m" -> missingAsState = true;
				case "--silent", "-s" -> silent = true;
				default -> {
					System.err.println("Unknown argument: " + args[i]);
					System.exit(2);
				}
			}
		}

		var taxa = new TaxaBlock();
		var characters = new CharactersBlock();
		var reader = ImportManager.getInstance().getReader(input);
		if (reader == null || !CharactersBlock.class.isAssignableFrom(reader.getToClass()))
			throw new java.io.IOException("Not an alignment: " + input);
		((splitstree6.io.utils.DataReaderBase<CharactersBlock>) reader).read(new ProgressSilent(), input, taxa, characters);
		System.out.printf("Input: %s%n  taxa=%d sites=%d%n", input, taxa.getNtax(), characters.getNchar());

		// the options are set by name, as the algorithms have no common supertype that declares them
		splitstree6.algorithms.characters.characters2network.Characters2Network algorithm = switch (engine.toLowerCase()) {
			case "tcs" -> new razornet.algorithms.TCSHaplotypeNetwork();
			default -> new razornet.algorithms.RazorHaploNet();
		};
		setOptionByName(algorithm, "MaxSeconds", maxSeconds);
		setOptionByName(algorithm, "RandomRestarts", restarts);
		setOptionByName(algorithm, "RandomSeed", seed);
		if (!strategy.isEmpty())
			setOptionByName(algorithm, "Strategy", strategy);
		setOptionByName(algorithm, "Labeling", labeling);
		setOptionByName(algorithm, "MissingAsState", missingAsState);
		if (raw || noRemove)
			setOptionByName(algorithm, "RemoveSuperfluous", false);
		if (raw || noContract)
			setOptionByName(algorithm, "ContractUnlabeled", false);
		// further options are set by name
		for (var entry : settings.entrySet()) {
			var option = splitstree6.options.Option.getName2Options(algorithm).get(entry.getKey().replaceAll("^option", ""));
			if (option == null) {
				System.err.println("Unknown option: " + entry.getKey() + "; this algorithm offers " + algorithm.listOptions());
				System.exit(2);
			}
			var value = (option.getOptionValueType() == splitstree6.options.OptionValueType.Enum
					? option.getEnumValueForName(entry.getValue())
					: splitstree6.options.OptionValueType.parseType(option.getOptionValueType(), entry.getValue()));
			if (value == null) {
				System.err.println("Bad value for " + entry.getKey() + ": " + entry.getValue()
								   + (option.getLegalValues() != null ? "; legal values are " + option.getLegalValues() : ""));
				System.exit(2);
			}
			((javafx.beans.property.Property<Object>) option.getProperty()).setValue(value);
		}
		System.out.printf("Engine: %s%n", algorithm.getClass().getSimpleName());

		var networkBlock = new NetworkBlock();
		// no workflow here, so the sites are numbered as given
		try (var progress = silent ? new ProgressSilent() : new ProgressPercentage()) {
			algorithm.compute(progress, taxa, characters, networkBlock);
		}

		var graph = networkBlock.getGraph();
		var length = graph.edgeStream().mapToDouble(graph::getWeight).sum();
		var auxiliary = (int) graph.nodeStream().filter(v -> !graph.hasTaxa(v)).count();
		// the summary line that the program shows under the drawing
		System.out.printf("  viewer info: %s%n",
				splitstree6.view.network.AdditionalConsoleOutput.infoFor(networkBlock, null, characters));
		System.out.printf("Result (labeling=%s%s):%n  length=%s nodes=%d edges=%d auxiliary=%d%n",
				labeling, raw ? ", raw: no removal/contraction" : (noRemove ? ", no removal" : (noContract ? ", no contraction" : "")),
				length == Math.rint(length) ? String.valueOf((long) length) : String.valueOf(length),
				graph.getNumberOfNodes(), graph.getNumberOfEdges(), auxiliary);

		// Per-edge agreement between the mutations drawn and the length displayed.
		var mismatch = new TreeMap<String, Integer>();
		int labelled = 0, sites = 0, agree = 0;
		double weightOfUnlabelled = 0;
		for (var e : graph.edges()) {
			var label = networkBlock.getEdgeData(e).get(NetworkBlock.EDGE_SITES_KEY);
			var drawn = (label == null || label.isBlank()) ? 0 : label.split(",").length;
			var w = graph.getWeight(e);
			if (drawn > 0) {
				labelled++;
				sites += drawn;
			} else
				weightOfUnlabelled += w;
			if (Math.abs(w - drawn) < 1e-9)
				agree++;
			else
				mismatch.merge(drawn < w ? "fewer mutations than length" : "more mutations than length", 1, Integer::sum);
		}
		System.out.printf("  edges with site labels=%d of %d, mutations drawn=%d%n", labelled, graph.getNumberOfEdges(), sites);
		System.out.printf("  edges where mutations == weight: %d of %d%n", agree, graph.getNumberOfEdges());
		for (var entry : mismatch.entrySet())
			System.out.printf("  edges with %s: %d%n", entry.getKey(), entry.getValue());
		// the excess and the shortfall of the final network, after edge removal and contraction
		reportRealization(taxa, characters, networkBlock, missingAsState);

		if (weightOfUnlabelled > 0)
			System.out.printf("  total weight carried by UNLABELED edges: %s%n",
					weightOfUnlabelled == Math.rint(weightOfUnlabelled) ? String.valueOf((long) weightOfUnlabelled) : String.valueOf(weightOfUnlabelled));
	}
}
