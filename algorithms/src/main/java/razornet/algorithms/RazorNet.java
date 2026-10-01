/*
 * RazorNet.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms;

import javafx.beans.property.*;
import jloda.fx.util.ProgramProperties;
import jloda.graph.Node;
import jloda.phylo.PhyloGraph;
import jloda.util.*;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.preprocessing.PreprocessingMethod;
import razornet.algorithms.utils.*;
import splitstree6.options.Option;
import splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances;
import splitstree6.algorithms.distances.distances2network.Distances2Network;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;

import java.util.*;

/**
 * RazorNet: computes a network that realizes the distances exactly, with as little total edge length as it can
 * find, using auxiliary vertices where they shorten the network.
 * <p>
 * The input is quantized to even integers and the taxa-only realization is recorded. A best-of-N random-restart
 * search then extends the distance matrix by auxiliary points, realizes it as a graph, polishes it and removes
 * redundant edges; each strategy's best graph gets a finishing polish; an optional stretch filter removes
 * further edges. The extension has two strategies ({@link Strategy}): Strict compactifies exactly the points
 * that lie on no geodesic between two others, Relaxed also points that lie between others; Both runs both, the
 * likely winner first, and keeps the shorter network. The result is exact unless the stretch filter is on.
 * <p>
 * Source: Hayamizu et al., RazorNet, in preparation; see {@code docs/algorithm.md}.
 *
 * @author Daniel Huson and Momoko Hayamizu, 2026
 */
public class RazorNet extends Distances2Network {
	/**
	 * The extension strategy: Strict compactifies exactly the points that lie on no geodesic between two others,
	 * all in one pass; Relaxed also compactifies points that lie between others, one per call; Both runs both and
	 * keeps the shorter network.
	 */
	public enum Strategy {Strict, Relaxed, Both}

	private boolean verbose = true;

	private static final int MAX_FINISH_ROUNDS = 100;   // cap on the cycles of a finishing polish
	private static final int FINISH_FALLBACK_SECONDS = 30; // deadline of a finishing polish when there is no time budget
	private static final double SEARCH_BUDGET_FRACTION = 0.95; // share of maxSeconds for the search; the rest is for the stretch filter and the final check
	// in Both, the first phase (restarts and finishing polish) may use this share of the search budget
	private static final double FIRST_PHASE_BUDGET_FRACTION = 0.8;
	// within a phase, the restarts may use this share of the phase's window; the finishing polish gets the rest
	private static final double RESTART_SHARE_OF_PHASE = 0.8;
	// in Both, the relaxed strategy runs first if at least this fraction of the taxa lie between two others
	private static final double RELAXED_FIRST_INTERIOR_FRACTION = 0.1;
	private static final int INTERIOR_COUNT_MAX_TAXA = 1000; // above this, the O(n^3) interior count is skipped and the strict strategy runs first
	private static final int STRETCH_UNBOUNDED_MAX_NODES = 2000; // up to this size, the stretch filter runs without a time limit
	private static final int EXPLOSION_NODES_PER_TAXON = 6;      // above this node-to-taxon ratio, the user is advised

	private final ObjectProperty<Strategy> optionStrategy = new SimpleObjectProperty<>(this, "optionStrategy", Strategy.Both);

	// optional transformation of the input distances before the search
	private final ObjectProperty<PreprocessingMethod> optionPreprocessing = new SimpleObjectProperty<>(this, "optionPreprocessing", PreprocessingMethod.None);

	// for PruneSmallSplits: remove the smallest splits making up this percent of the total split weight
	private final DoubleProperty optionSplitFilterPercent = new SimpleDoubleProperty(this, "optionSplitFilterPercent", 1.0);

	private final IntegerProperty optionSignificantDigits = new SimpleIntegerProperty(this, "optionSignificantDigits", 6);
	// restarts per strategy
	private final IntegerProperty optionRandomRestarts = new SimpleIntegerProperty(this, "optionRandomRestarts", 1000);

	// Strict strategy: the restarts run in batches of this size, and the phase ends when a batch does not shorten the
	// best network (0 = run all restarts)
	private final IntegerProperty optionMaxNonImprovingRestarts = new SimpleIntegerProperty(this, "optionMaxNonImprovingRestarts", 3);

	private final BooleanProperty optionRemoveRedundant = new SimpleBooleanProperty(this, "optionRemoveRedundant", true);

	// lossy post-processing: remove edges while no pairwise distance grows by more than this percent (0 = off)
	private final DoubleProperty optionMaxStretchPercent = new SimpleDoubleProperty(this, "optionMaxStretchPercent", 0.0);

	// extend+polish cycles per restart
	private final IntegerProperty optionMaxPolishRounds = new SimpleIntegerProperty(this, "optionMaxPolishRounds", 1);

	// after the search, polish each strategy's best graph until its length stops decreasing
	private final BooleanProperty optionFinishingPolish = new SimpleBooleanProperty(this, "optionFinishingPolish", true);

	// if set, an analysis log is written to this file; not saved in the program properties
	private final StringProperty optionLogFile = new SimpleStringProperty(this, "optionLogFile", "");

	private final IntegerProperty optionRandomSeed = new SimpleIntegerProperty(this, "optionRandomSeed", 0);

	private final IntegerProperty optionMaxSeconds = new SimpleIntegerProperty(this, "optionMaxSeconds", 60);

	/** prefix of the keys under which the options are saved in the program properties */
	private static final String KEY = "razornet.algorithm.RazorNet";

	{
		ProgramProperties.track(KEY + optionSignificantDigits.getName(), optionSignificantDigits, 6);
		ProgramProperties.track(KEY + optionRandomRestarts.getName(), optionRandomRestarts, 1000);
		ProgramProperties.track(KEY + optionMaxNonImprovingRestarts.getName(), optionMaxNonImprovingRestarts, 3);
		ProgramProperties.track(KEY + optionMaxStretchPercent.getName(), optionMaxStretchPercent, 0.0);
		ProgramProperties.track(KEY + optionSplitFilterPercent.getName(), optionSplitFilterPercent, 1.0);
		ProgramProperties.track(KEY + optionRandomSeed.getName(), optionRandomSeed, 0);
		ProgramProperties.track(KEY + optionMaxSeconds.getName(), optionMaxSeconds, 60);
		ProgramProperties.track(KEY + optionFinishingPolish.getName(), optionFinishingPolish, true);
		// not tracked, so they start at their defaults in every session: optionStrategy and optionPreprocessing
		// (enums), optionMaxPolishRounds, optionRemoveRedundant and optionLogFile
	}

	@Override
	public String getShortDescription() {
		return "Computes a network that realizes the distances exactly, with as little total edge length as it can find.";
	}

	/**
	 * the basic options and, in advanced mode ({@link AdvancedOptions}), a separator followed by the advanced options
	 */
	@Override
	public List<String> listOptions() {
		var list = new ArrayList<>(List.of(optionPreprocessing.getName(), optionStrategy.getName(), optionMaxStretchPercent.getName(),
				optionSignificantDigits.getName(), optionMaxSeconds.getName()));
		if (AdvancedOptions.enabled()) {
			list.add(Option.SEPARATOR);
			list.addAll(List.of(optionRandomRestarts.getName(), optionMaxNonImprovingRestarts.getName(), optionMaxPolishRounds.getName(),
					optionFinishingPolish.getName(), optionSplitFilterPercent.getName(), optionRandomSeed.getName(), optionLogFile.getName()));
		}
		return list;
	}

	/** all option properties, basic and advanced */
	private List<javafx.beans.property.Property<?>> allOptionProperties() {
		return List.of(optionPreprocessing, optionStrategy, optionMaxStretchPercent, optionSignificantDigits, optionMaxSeconds,
				optionRandomRestarts, optionMaxNonImprovingRestarts, optionMaxPolishRounds, optionFinishingPolish,
				optionSplitFilterPercent, optionRandomSeed, optionLogFile);
	}

	@Override
	public String getToolTip(String optionName) {
		if (optionStrategy.getName().equals(optionName))
			return "which extension strategy drives the search: Strict compactifies every taxon that lies on no geodesic between "
				   + "two others (nearly deterministic, provably safe at the top level); Relaxed also compactifies taxa that lie between "
				   + "others, one per call, and explores the many ties with random restarts (what grids and Hamming distances with "
				   + "sampled intermediates need); Both (default) runs both, the likely winner first (Relaxed when at least 10% of "
				   + "the taxa lie between two others), and keeps the shorter network";
		else if (optionMaxNonImprovingRestarts.getName().equals(optionName))
			return "Strict strategy: its restarts run in batches of this size, and the phase stops when a batch fails to shorten "
				   + "the best network (0 = run all restarts). The Relaxed strategy always runs all restarts";
		else if (optionPreprocessing.getName().equals(optionName))
			return "optional distance preprocessing before the search. Classical MDS to 2-D/3-D turns messy, non-Euclidean distances such as driving distances into an embeddable metric (pair with low significant digits and a max stretch for a sparse, map-like network). PruneSmallSplits denoises via split decomposition, dropping the weakest splits (see split filter percent) so a generic metric no longer over-expands the search";
		else if (optionSplitFilterPercent.getName().equals(optionName))
			return "for the PruneSmallSplits preprocessing: remove the smallest splits that together make up this percent of the total split weight (default 1)";
		else if (optionMaxPolishRounds.getName().equals(optionName))
			return "number of extend+polish cycles per random restart (default 1: under a time budget, more restarts beat deeper polishing)";
		else if (optionLogFile.getName().equals(optionName))
			return "diagnostic: if set to a file path, write an analysis log there -- the initial slack values under both definitions, a recursion trace, the distribution of lengths over the restarts, what each strategy achieved, and the result. Leave empty for no log";
		else if (optionFinishingPolish.getName().equals(optionName))
			return "after the restart search, keep polishing the single best network until its length stops decreasing; recovers metrics that need several extend+polish cycles (default on)";
		else if (optionMaxStretchPercent.getName().equals(optionName))
			return "lossy post-processing: remove edges while keeping every pairwise distance within this percent of the input (0 = keep distances exact)";
		else if (optionRandomRestarts.getName().equals(optionName))
			return "number of random restarts per strategy (the Strict strategy may stop earlier, see max non-improving restarts)";
		else if (optionSignificantDigits.getName().equals(optionName))
			return "number of significant digits kept when the distances are scaled to integers; fewer digits give a coarser, sparser network";
		else if (optionRandomSeed.getName().equals(optionName))
			return "random number seed";
		else if (optionMaxSeconds.getName().equals(optionName))
			return "maximum time in seconds to search before returning the best result found so far (0 = no limit)";
		else return super.getToolTip(optionName);
	}

	/**
	 * One strategy's core, as the pipeline uses it: the cap on the working matrix, the extension, one round of
	 * polishing, and the diagnostic counters. Implemented by delegation to {@code razornet.algorithms.strict} and
	 * {@code razornet.algorithms.relaxed}.
	 */
	private interface Core {
		Strategy strategy();

		int extensionCap(int groundSize);

		void extend(MutableD distances, Random random, ProgressListener progress, int extensionCap) throws CanceledException;

		boolean polishRound(IntGraph graph, MutableD distances, Random random, ProgressListener progress) throws CanceledException;

		void resetCounters();

		String summary();
	}

	private static final Core STRICT = new Core() {
		public Strategy strategy() {
			return Strategy.Strict;
		}

		public int extensionCap(int groundSize) {
			return razornet.algorithms.strict.DistanceMatrixExtension.extensionCap(groundSize);
		}

		public void extend(MutableD distances, Random random, ProgressListener progress, int extensionCap) throws CanceledException {
			razornet.algorithms.strict.DistanceMatrixExtension.extend(distances, random, progress, extensionCap);
		}

		public boolean polishRound(IntGraph graph, MutableD distances, Random random, ProgressListener progress) throws CanceledException {
			return razornet.algorithms.strict.PolishGraph.applyRound(graph, distances, random, progress);
		}

		public void resetCounters() {
			razornet.algorithms.strict.DistanceMatrixExtension.resetCounters();
		}

		public String summary() {
			return razornet.algorithms.strict.DistanceMatrixExtension.summary();
		}
	};

	private static final Core RELAXED = new Core() {
		public Strategy strategy() {
			return Strategy.Relaxed;
		}

		public int extensionCap(int groundSize) {
			return razornet.algorithms.relaxed.DistanceMatrixExtension.extensionCap(groundSize);
		}

		public void extend(MutableD distances, Random random, ProgressListener progress, int extensionCap) throws CanceledException {
			razornet.algorithms.relaxed.DistanceMatrixExtension.extend(distances, random, progress, extensionCap);
		}

		public boolean polishRound(IntGraph graph, MutableD distances, Random random, ProgressListener progress) throws CanceledException {
			return razornet.algorithms.relaxed.PolishGraph.applyRound(graph, distances, random, progress);
		}

		public void resetCounters() {
			razornet.algorithms.relaxed.DistanceMatrixExtension.resetCounters();
		}

		public String summary() {
			return razornet.algorithms.relaxed.DistanceMatrixExtension.summary();
		}
	};

	private static Core coreOf(Strategy strategy) {
		return strategy == Strategy.Strict ? STRICT : RELAXED;
	}

	/**
	 * The state of the search: the best graph so far, its length in input units and the strategy that produced it
	 * (null for the taxa-only realization), and each strategy's own best graph. Every graph passed to
	 * {@link #consider} realizes the input, so an interrupted run still returns a valid network.
	 */
	private static final class Search {
		final int[][] quantMatrix;
		final Quantization quantization;
		private double bestLength = Double.MAX_VALUE;
		private IntGraph bestGraph = null;
		private Strategy bestSource = null;
		private final EnumMap<Strategy, Double> lengthBySource = new EnumMap<>(Strategy.class);
		private final EnumMap<Strategy, IntGraph> graphBySource = new EnumMap<>(Strategy.class);

		Search(int[][] quantMatrix, Quantization quantization) {
			this.quantMatrix = quantMatrix;
			this.quantization = quantization;
		}

		/**
		 * thread-safe update: keeps a copy of the graph if it is strictly shorter than the best so far, overall and
		 * for its strategy
		 */
		synchronized void consider(IntGraph graph, Strategy source) {
			if (graph == null)
				return;
			var length = graph.getTotalLength(quantization.mapDistanceBack()); // summed in doubles to avoid int overflow
			if (length < bestLength) {
				bestLength = length;
				bestGraph = graph.copy();
				bestSource = source;
			}
			if (source != null && length < lengthBySource.getOrDefault(source, Double.MAX_VALUE)) {
				lengthBySource.put(source, length);
				graphBySource.put(source, graph.copy());
			}
		}

		synchronized double bestLength() {
			return bestLength;
		}

		synchronized IntGraph bestGraph() {
			return bestGraph;
		}

		synchronized Strategy bestSource() {
			return bestSource;
		}

		/** the shortest length this strategy reached so far, or MAX_VALUE if it produced nothing */
		synchronized double bestLength(Strategy source) {
			return lengthBySource.getOrDefault(source, Double.MAX_VALUE);
		}

		/** this strategy's own best graph (a private copy), or null if it produced nothing */
		synchronized IntGraph bestGraph(Strategy source) {
			var graph = graphBySource.get(source);
			return graph == null ? null : graph.copy();
		}
	}

	/** What one phase of the search did, for the report and the analysis log. */
	private static final class PhaseResult {
		final Strategy strategy;
		final boolean sharedBudget;          // running under a capped share of the search budget (inside Both)
		private int restarts = 0;            // restarts actually started (a deadline skips the rest of a batch)
		double bestLength = Double.MAX_VALUE; // the phase's own best, in input units, whether or not it beat the floor
		double seconds = 0;
		boolean stoppedEarly = false;
		boolean deadlineHit = false;

		PhaseResult(Strategy strategy, boolean sharedBudget) {
			this.strategy = strategy;
			this.sharedBudget = sharedBudget;
		}

		synchronized void started() {
			restarts++;
		}

		synchronized int restarts() {
			return restarts;
		}

		synchronized void record(double length) {
			bestLength = Math.min(bestLength, length);
		}

		String describe() {
			var buf = new StringBuilder();
			buf.append("%-12s %d restarts, best %s (%.2f s".formatted(strategy + ":", restarts(),
					bestLength == Double.MAX_VALUE ? "none" : StringUtils.trim((float) bestLength), seconds));
			if (stoppedEarly)
				buf.append("; stopped: a whole batch brought no improvement");
			if (deadlineHit)
				buf.append(sharedBudget ? "; its share of the time budget was used up" : "; time budget reached");
			buf.append(")");
			return buf.toString();
		}
	}

	@Override
	public void compute(ProgressListener progress, TaxaBlock taxaBlock, DistancesBlock distancesBlock, NetworkBlock networkBlock) {
		System.err.println("Running RazorNet (strategy: " + getOptionStrategy() + ")");
		progress.setTasks("RazorNet", "Initialization");

		// time budget: the search may use SEARCH_BUDGET_FRACTION of maxSeconds (0 = no limit); once the deadline is
		// reached, checkForCancel() throws and the best graph found so far is returned
		var timeoutProgress = new TimeoutProgress(progress, SEARCH_BUDGET_FRACTION * getOptionMaxSeconds());

		// optional preprocessing; the network realizes the preprocessed distances
		var inputDistances = getOptionPreprocessing().newPreprocessor(getOptionSplitFilterPercent() / 100.0).apply(distancesBlock.getDistances());

		// a non-metric input has no realization: warn, but continue
		MetricCheck.checkAndWarn(inputDistances, "RazorNet");

		// tracing for the analysis log has to start before the restarts
		final var logging = !getOptionLogFile().isBlank();
		final var logStart = System.currentTimeMillis();
		if (logging)
			ExtensionTrace.start();

		var quantization = Quantization.apply(inputDistances, getOptionSignificantDigits());
		var canonicalPermutation = MatrixPermutation.canonical(quantization.quantMatrix());
		var quantMatrix = canonicalPermutation.apply(quantization.quantMatrix());

		var search = new Search(quantMatrix, quantization);

		STRICT.resetCounters();
		RELAXED.resetCounters();

		// the taxa-only realization: exact, without auxiliary vertices; the search replaces it only by something shorter
		progress.setSubtask("Taxa-only realization");
		try {
			search.consider(GraphFromDistances.realizeTaxaOnly(quantMatrix, timeoutProgress), null);
		} catch (CanceledException ignored) {
			// time budget exhausted; the search still runs
		}

		// ----- the search: one strategy, or both, the likely winner first -----
		// Both: the relaxed strategy runs first if at least RELAXED_FIRST_INTERIOR_FRACTION of the taxa lie between
		// two others, otherwise the strict strategy. Under a time budget the first phase may use
		// FIRST_PHASE_BUDGET_FRACTION of the search budget and the second the rest; within a phase the restarts may
		// use RESTART_SHARE_OF_PHASE of the window and the finishing polish what remains.
		var phases = new ArrayList<PhaseResult>();
		var strategy = getOptionStrategy();
		var order = new ArrayList<Core>();
		String orderLine = null;
		if (strategy == Strategy.Both) {
			var n = quantMatrix.length;
			var interior = n <= INTERIOR_COUNT_MAX_TAXA ? SlackStructure.countInterior(quantMatrix) : -1;
			var relaxedFirst = interior >= RELAXED_FIRST_INTERIOR_FRACTION * n;
			order.add(relaxedFirst ? RELAXED : STRICT);
			order.add(relaxedFirst ? STRICT : RELAXED);
			orderLine = interior < 0
					? "Both: %,d taxa, too many to count the interior ones; running the strict strategy (Strict) first".formatted(n)
					: "Both: %d of %d taxa lie between two others (%.0f%%, threshold %.0f%%); running the %s strategy first".formatted(
					interior, n, 100.0 * interior / n, 100.0 * RELAXED_FIRST_INTERIOR_FRACTION, order.get(0).strategy());
			System.err.println(orderLine);
		} else
			order.add(coreOf(strategy));

		var beforeFinishing = new EnumMap<Strategy, Double>(Strategy.class);
		for (var i = 0; i < order.size(); i++) {
			var core = order.get(i);
			var first = strategy == Strategy.Both && i == 0;
			// the time window of this phase, and the share of it that the restarts may use
			TimeoutProgress windowProgress;
			TimeoutProgress restartProgress;
			if (getOptionMaxSeconds() > 0) {
				var remaining = timeoutProgress.secondsRemaining();
				var window = first ? Math.min(remaining, FIRST_PHASE_BUDGET_FRACTION * SEARCH_BUDGET_FRACTION * getOptionMaxSeconds()) : remaining;
				windowProgress = new TimeoutProgress(timeoutProgress, window);
				restartProgress = new TimeoutProgress(windowProgress, RESTART_SHARE_OF_PHASE * window);
			} else {
				windowProgress = timeoutProgress;
				restartProgress = timeoutProgress;
			}
			// the strict strategy stops after a batch of restarts without improvement; the relaxed strategy runs them all
			var adaptive = core == STRICT && getOptionMaxNonImprovingRestarts() > 0;
			var batch = adaptive ? getOptionMaxNonImprovingRestarts() : getOptionRandomRestarts();
			phases.add(runPhase(core, search, batch, adaptive, restartProgress, first, progress));

			// finishing polish: polish this strategy's own best graph with its own core until its length stops
			// decreasing; only if the graph has auxiliary vertices
			beforeFinishing.put(core.strategy(), search.bestLength(core.strategy()));
			if (isOptionFinishingPolish() && getOptionMaxPolishRounds() > 0) {
				var start = search.bestGraph(core.strategy());
				if (start != null && start.getNumberOfNodes() > quantMatrix.length && !windowProgress.isDeadlineReached()) {
					progress.setSubtask("Finishing polish (" + core.strategy() + ")");
					var finishProgress = getOptionMaxSeconds() > 0 ? windowProgress : new TimeoutProgress(progress, FINISH_FALLBACK_SECONDS);
					finishingPolish(core, start, search, finishProgress);
				}
			}
		}

		if (search.bestGraph() == null) {
			if (timeoutProgress.isDeadlineReached())
				throw new IllegalStateException(("RazorNet: no network produced within the time budget of %d s — increase the 'Max Seconds' option (0 = no limit).").formatted(getOptionMaxSeconds()));
			throw new IllegalStateException("Failed to compute graph");
		}
		if (timeoutProgress.isDeadlineReached())
			System.err.printf("RazorNet: time budget of %d s reached — returning best result found so far (length %s).%n",
					getOptionMaxSeconds(), StringUtils.trim((float) search.bestLength()));

		// ----- report what the strategies did -----
		var runSummary = new StringBuilder();
		if (orderLine != null)
			runSummary.append(orderLine).append(System.lineSeparator());
		for (var phase : phases)
			runSummary.append(phase.describe()).append(System.lineSeparator());
		if (!phases.isEmpty()) {
			var parts = new ArrayList<String>();
			for (var phase : phases)
				parts.add("%s %s -> %s".formatted(phase.strategy, lengthText(beforeFinishing.get(phase.strategy)),
						lengthText(search.bestLength(phase.strategy))));
			var finishingLine = "after finishing polish: " + String.join(", ", parts);
			System.err.println(finishingLine);
			runSummary.append(finishingLine).append(System.lineSeparator());
		}
		var winners = new ArrayList<String>();
		for (var phase : phases)
			if (search.bestLength(phase.strategy) == search.bestLength())
				winners.add(phase.strategy.toString());
		var winnerLine = winners.isEmpty()
				? "winner: the taxa-only realization, length %s (neither strategy improved on it)".formatted(lengthText(search.bestLength()))
				: winners.size() == 1
				? "winner: %s strategy, length %s".formatted(winners.get(0), lengthText(search.bestLength()))
				: "winner: tie between %s, length %s".formatted(String.join(" and ", winners), lengthText(search.bestLength()));
		System.err.println(winnerLine);
		runSummary.append(winnerLine);

		// optional lossy post-processing of the best graph (in permuted, quantized coordinates); if it is cancelled,
		// the edges removed so far stay removed
		if (getOptionMaxStretchPercent() > 0) {
			// the stretch filter runs on its own clock: without a limit for graphs of up to STRETCH_UNBOUNDED_MAX_NODES
			// nodes or when there is no time budget, otherwise with a budget of maxSeconds
			ProgressListener stretchProgress;
			if (getOptionMaxSeconds() <= 0 || search.bestGraph().getNumberOfNodes() <= STRETCH_UNBOUNDED_MAX_NODES)
				stretchProgress = progress;
			else
				stretchProgress = new TimeoutProgress(progress, getOptionMaxSeconds());
			try {
				GreedyStretchEdgeRemoval.apply(search.bestGraph(), quantization.isLabeled(),
						(u, v) -> quantMatrix[u][v], getOptionMaxStretchPercent() / 100.0, stretchProgress);
			} catch (CanceledException ignored) {
				// cancelled or out of time: keep the reduction achieved
			}
		}

		var finalGraph = canonicalPermutation.unpermuteGraph(search.bestGraph(), quantization.size());
		var network = networkBlock.getGraph();

		setupNetwork(taxaBlock, finalGraph, quantization, network);
		// a DistanceNetwork: the edge weights are lengths that realize the input distances
		networkBlock.setNetworkType(NetworkBlock.Type.DistanceNetwork);
		for (var e : network.edges()) {
			networkBlock.getEdgeData(e).put("weight", StringUtils.trim(network.getWeight(e)));
		}

		// final check of the pairwise distances; it throws on a disconnected graph, which must not abort the run
		try {
			CheckPairwiseDistances.apply(network, inputDistances, 0.000001);
		} catch (RuntimeException ex) {
			System.err.println("Note: pairwise-distance check skipped (" + ex.getMessage() + ")");
		}

		var counters = new StringBuilder();
		if (strategy != Strategy.Relaxed)
			counters.append("Strict extension: ").append(STRICT.summary()).append(System.lineSeparator());
		if (strategy != Strategy.Strict)
			counters.append("Relaxed extension:   ").append(RELAXED.summary()).append(System.lineSeparator());
		System.err.print(counters);

		if (logging) {
			ExtensionTrace.stop();
			AnalysisLog.write(getOptionLogFile(), getClass().getSimpleName(), optionsSummary(), taxaBlock, inputDistances,
					getOptionSignificantDigits(), network, (System.currentTimeMillis() - logStart) / 1000.0, Double.NaN, true,
					runSummary + System.lineSeparator() + counters);
		}

		// many more nodes than taxa: the distances have few ties at this precision, so advise the user
		var ntax = taxaBlock.getNtax();
		if (network.getNumberOfNodes() > EXPLOSION_NODES_PER_TAXON * ntax) {
			var lever = getOptionSignificantDigits() > 2
					? "try fewer significant digits (e.g. 2 or 3)"
					: "try a max stretch, or MDS preprocessing";
			System.err.printf("RazorNet: the network has %,d nodes for %,d taxa --- the distances look generic (few "
							+ "ties) at %d significant digits, so the extension over-expands. For a sparser network, %s.%n",
					network.getNumberOfNodes(), ntax, getOptionSignificantDigits(), lever);
		}
	}

	/**
	 * Runs one strategy's restarts, in parallel batches, recording every checkpoint in {@code search}. The restart
	 * seeds are drawn from {@code new Random(optionRandomSeed)} for each phase, so that a strategy executes the same
	 * restarts inside {@code Both} as alone.
	 *
	 * @param batchSize            restarts per parallel batch
	 * @param stopWhenNotImproving end the phase after a batch that does not shorten the best network
	 * @param timeoutProgress      the deadline of this phase
	 * @param sharedBudget         whether that deadline is a share of the search budget (inside Both)
	 */
	private PhaseResult runPhase(Core core, Search search, int batchSize, boolean stopWhenNotImproving,
								 TimeoutProgress timeoutProgress, boolean sharedBudget, ProgressListener progress) {
		var result = new PhaseResult(core.strategy(), sharedBudget);
		var start = System.currentTimeMillis();
		var extensionCap = core.extensionCap(search.quantMatrix.length);
		var numberOfRestarts = getOptionRandomRestarts();

		var seeds = new ArrayList<Long>();
		{
			var random = new Random(getOptionRandomSeed());
			for (var i = 0; i < numberOfRestarts; i++)
				seeds.add(random.nextLong());
		}

		var label = core.strategy() + (getOptionMaxPolishRounds() > 0 ? ": Extend+Polish" : ": Extend");
		progress.setSubtask(label);
		if (verbose)
			System.err.println("----------- " + label + " -----------");

		try {
			var next = 0;
			while (next < seeds.size() && !timeoutProgress.isDeadlineReached()) {
				var end = Math.min(seeds.size(), next + Math.max(1, batchSize));
				var batch = seeds.subList(next, end);
				var bestBefore = search.bestLength();
				ExecuteInParallel.apply(batch, seed -> runRestart(core, search, seed, extensionCap, timeoutProgress, result),
						ProgramExecutorService.getNumberOfCoresToUse(), timeoutProgress); // timeoutProgress drives the time bar; its job-count updates are swallowed
				next = end;
				if (stopWhenNotImproving && !(search.bestLength() < bestBefore) && next < seeds.size()) {
					result.stoppedEarly = true;
					break;
				}
			}
		} catch (Exception e) {
			Basic.caught(e);
			throw new RuntimeException(e);
		}
		result.deadlineHit = timeoutProgress.isDeadlineReached();
		result.seconds = (System.currentTimeMillis() - start) / 1000.0;
		if (verbose)
			System.err.println(result.describe());
		return result;
	}

	/**
	 * One restart: extend, realize, polish (up to {@code optionMaxPolishRounds} cycles), remove redundant edges,
	 * recording the graph in {@code search} after every step. A restart cut short by the deadline simply ends.
	 */
	private void runRestart(Core core, Search search, long seed, int extensionCap, TimeoutProgress timeoutProgress, PhaseResult phase) throws CanceledException {
		if (timeoutProgress.isDeadlineReached())
			return; // time budget already spent: skip any remaining restarts
		phase.started();
		var quantMatrix = search.quantMatrix;
		var quantization = search.quantization;
		try {
			final var distances = new MutableD(quantMatrix);
			IntGraph intGraph = null;

			var random = new Random(seed);

			var previousLength = Long.MAX_VALUE;
			for (var round = 0; round < Math.max(1, getOptionMaxPolishRounds()); round++) {
				var sizeBefore = distances.size();
				core.extend(distances, random, timeoutProgress, extensionCap);
				var extendChanged = distances.size() > sizeBefore;

				intGraph = GraphFromDistances.realizeFromDistances(distances::get, SetUtils.fullIndexSet(distances.size()), timeoutProgress);
				// checkpoint: the graph realizes the distances, so it is a valid result
				search.consider(intGraph, core.strategy());
				var polishChanged = getOptionMaxPolishRounds() > 0 && core.polishRound(intGraph, distances, random, timeoutProgress);

				var currentLength = intGraph.getTotalLength();
				if (currentLength >= previousLength)
					break;   // worse -> stop
				previousLength = currentLength;

				if (!extendChanged && !polishChanged) break;
				timeoutProgress.checkForCancel();
			}

			// checkpoint before the redundant-edge removal
			search.consider(intGraph, core.strategy());
			if (getOptionRemoveRedundant()) {
				GreedyGlobalEdgeRemoval.apply(intGraph, quantization.isLabeled(), (u, v) -> quantMatrix[u][v], timeoutProgress);
				search.consider(intGraph, core.strategy());
			}

			// the length this restart achieved, after edge removal
			var achieved = intGraph.getTotalLength(quantization.mapDistanceBack());
			ExtensionTrace.recordRestartLength(achieved);
			phase.record(achieved);

			if (verbose) {
				System.err.printf("Length achieved: %s (current best length: %s)%n",
						StringUtils.trim(achieved), StringUtils.trim((float) search.bestLength()));
			}
		} catch (TimeoutProgress.DeadlineReachedException timedOut) {
			// deadline reached: the best result so far is already recorded
		}
	}

	/**
	 * The finishing polish of one strategy: extend, polish and remove redundant edges on a copy of its best graph,
	 * with its own core, until a cycle does not shorten the graph or the deadline is reached.
	 */
	private void finishingPolish(Core core, IntGraph start, Search search, TimeoutProgress finishProgress) {
		var quantMatrix = search.quantMatrix;
		var quantization = search.quantization;
		try {
			var finishCap = Math.min(core.extensionCap(quantMatrix.length), start.getNumberOfNodes() + quantMatrix.length);
			var workM = matrixFromGraph(start, quantMatrix.length);
			var random = new Random(getOptionRandomSeed());
			var prevLength = (double) start.getTotalLength(); // incumbent (quantized units)
			for (var round = 0; round < MAX_FINISH_ROUNDS; round++) {
				if (finishProgress.isDeadlineReached())
					break;
				var sizeBefore = workM.size();
				core.extend(workM, random, finishProgress, finishCap);
				var extendChanged = workM.size() > sizeBefore;
				var work = GraphFromDistances.realizeFromDistances(workM::get, SetUtils.fullIndexSet(workM.size()), finishProgress);
				var polishChanged = core.polishRound(work, workM, random, finishProgress);
				if (getOptionRemoveRedundant())
					GreedyGlobalEdgeRemoval.apply(work, quantization.isLabeled(), (u, v) -> quantMatrix[u][v], finishProgress);
				search.consider(work, core.strategy());
				var curLength = (double) work.getTotalLength();
				if (curLength >= prevLength || (!extendChanged && !polishChanged))
					break; // did not beat the incumbent / converged
				prevLength = curLength;
			}
		} catch (CanceledException ignored) {
			// cancelled or out of time: the best so far is already recorded
		}
	}

	private static String lengthText(double length) {
		return length == Double.MAX_VALUE ? "none" : StringUtils.trim((float) length);
	}

	/**
	 * The distance matrix of a graph, as its all-pairs shortest paths, with the node ids made contiguous: taxa keep
	 * the ids 0..numTaxa-1, the auxiliary vertices get numTaxa, numTaxa+1, ...
	 */
	private static MutableD matrixFromGraph(IntGraph graph, int numTaxa) {
		var map = new HashMap<Integer, Integer>();
		var auxiliary = new ArrayList<Integer>();
		for (var v : graph.nodes()) {
			if (v < numTaxa) map.put(v, v);
			else auxiliary.add(v);
		}
		auxiliary.sort(null);
		var next = numTaxa;
		for (var v : auxiliary) map.put(v, next++);
		var matrix = new MutableD(next);
		var sp = DijkstraAllPairs.allShortestPathLengths(graph);
		for (var u : graph.nodes())
			for (var v : graph.nodes())
				matrix.set(map.get(u), map.get(v), sp.get(u).get(v));
		return matrix;
	}

	/** Copies the integer graph into the output network, re-attaching taxa (coincident taxa hang off their point by a zero edge). */
	private void setupNetwork(TaxaBlock taxaBlock, IntGraph intGraph, Quantization quantization, PhyloGraph network) {
		var nodeMap = new HashMap<Integer, Node>();
		for (var id : intGraph.nodes()) {
			if (!nodeMap.containsKey(id)) {
				var v = network.newNode();
				nodeMap.put(id, v);
				if (quantization.isLabeled().test(id)) {
					var taxa = quantization.mapNodeBack().apply(id);
					if (taxa.size() == 1) {
						var t = taxa.iterator().next() + 1;
						network.addTaxon(v, t);
						network.setLabel(v, taxaBlock.getLabel(t));
					} else if (taxa.size() > 1) {
						for (var t0 : taxa) {
							var t = t0 + 1;
							var w = network.newNode();
							network.addTaxon(w, t);
							network.setLabel(w, taxaBlock.getLabel(t));
							network.setWeight(network.newEdge(w, v), 0);
						}
					}
				}
			}
		}
		for (var e : intGraph.edges()) {
			var f = network.newEdge(nodeMap.get(e.u()), nodeMap.get(e.v()));
			// use getEdgeWeights() here: it sets up the network's edge weights
			network.getEdgeWeights().put(f, quantization.mapDistanceBack().applyAsDouble(intGraph.getWeight(e)));
		}
	}

	/** all options and their values, for the analysis log */
	private String optionsSummary() {
		var buf = new StringBuilder();
		for (var property : allOptionProperties())
			buf.append(buf.isEmpty() ? "" : ", ").append(property.getName().replaceAll("^option", ""))
					.append("=").append(property.getValue());
		return buf.toString();
	}

	public Strategy getOptionStrategy() {
		return optionStrategy.get();
	}

	public ObjectProperty<Strategy> optionStrategyProperty() {
		return optionStrategy;
	}

	public void setOptionStrategy(Strategy strategy) {
		optionStrategy.set(strategy);
	}

	public int getOptionMaxNonImprovingRestarts() {
		return optionMaxNonImprovingRestarts.get();
	}

	public IntegerProperty optionMaxNonImprovingRestartsProperty() {
		return optionMaxNonImprovingRestarts;
	}

	public int getOptionMaxPolishRounds() {
		return optionMaxPolishRounds.get();
	}

	public IntegerProperty optionMaxPolishRoundsProperty() {
		return optionMaxPolishRounds;
	}

	public int getOptionSignificantDigits() {
		return optionSignificantDigits.get();
	}

	public IntegerProperty optionSignificantDigitsProperty() {
		return optionSignificantDigits;
	}

	public PreprocessingMethod getOptionPreprocessing() {
		return optionPreprocessing.get();
	}

	public ObjectProperty<PreprocessingMethod> optionPreprocessingProperty() {
		return optionPreprocessing;
	}

	public double getOptionSplitFilterPercent() {
		return optionSplitFilterPercent.get();
	}

	public DoubleProperty optionSplitFilterPercentProperty() {
		return optionSplitFilterPercent;
	}

	public boolean getOptionRemoveRedundant() {
		return optionRemoveRedundant.get();
	}

	public BooleanProperty optionRemoveRedundantProperty() {
		return optionRemoveRedundant;
	}

	public boolean isOptionFinishingPolish() {
		return optionFinishingPolish.get();
	}

	public BooleanProperty optionFinishingPolishProperty() {
		return optionFinishingPolish;
	}

	public String getOptionLogFile() {
		return optionLogFile.get();
	}

	public StringProperty optionLogFileProperty() {
		return optionLogFile;
	}

	public void setOptionLogFile(String value) {
		optionLogFile.set(value);
	}

	public double getOptionMaxStretchPercent() {
		return optionMaxStretchPercent.get();
	}

	public DoubleProperty optionMaxStretchPercentProperty() {
		return optionMaxStretchPercent;
	}

	public int getOptionRandomRestarts() {
		return optionRandomRestarts.get();
	}

	public IntegerProperty optionRandomRestartsProperty() {
		return optionRandomRestarts;
	}

	public int getOptionRandomSeed() {
		return optionRandomSeed.get();
	}

	public IntegerProperty optionRandomSeedProperty() {
		return optionRandomSeed;
	}

	public int getOptionMaxSeconds() {
		return optionMaxSeconds.get();
	}

	public IntegerProperty optionMaxSecondsProperty() {
		return optionMaxSeconds;
	}
}
