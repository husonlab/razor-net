/*
 * RazorHaploNet.java Copyright (C) 2026 Daniel H. Huson
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
import jloda.graph.Edge;
import jloda.graph.Node;
import jloda.util.StringUtils;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.utils.GeodesicBetweenness;
import razornet.algorithms.utils.ParsimonyLabeler;
import splitstree6.algorithms.characters.characters2network.Characters2Network;
import splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances;
import splitstree6.algorithms.distances.distances2network.Distances2Network;
import splitstree6.data.*;
import splitstree6.view.alignment.AlignmentView;
import splitstree6.view.network.NetworkSequencesAnalyzer;

import java.io.IOException;
import java.util.*;
import java.util.function.IntFunction;
import java.util.function.ToDoubleFunction;

/**
 * RazorHaploNet: computes a haplotype network from an alignment. The Hamming distances of the sequences are
 * realized by {@link RazorNet}; sequences are then assigned to the auxiliary vertices by weighted parsimony
 * ({@link ParsimonyLabeler}), the weight of every edge is set to its number of mutations, superfluous edges are
 * removed and unlabeled edges are contracted. Columns with missing or ambiguous data are masked; a gap is a state,
 * unless the option Mask Gap Sites masks its column too, as PopART does. Sequences that are identical on the sites
 * kept are one haplotype: the engine, the labeling and the cleanup see one taxon per haplotype, and the other
 * sequences of a haplotype hang off its node by a zero-length edge in the end.
 * <p>
 * The engine that realizes the distances is obtained from {@link #getEngine()}, which
 * {@link TCSHaplotypeNetwork} overrides.
 * <p>
 * Source: Hayamizu et al., RazorNet, in preparation; see {@code docs/algorithm.md}.
 *
 * @author Daniel Huson and Momoko Hayamizu, 2026
 */
public class RazorHaploNet extends Characters2Network {
	private final BooleanProperty optionContractUnlabeled = new SimpleBooleanProperty(this, "optionContractUnlabeled", true);
	private final BooleanProperty optionRemoveSuperfluous = new SimpleBooleanProperty(this, "optionRemoveSuperfluous", true);

	private final IntegerProperty optionSignificantDigits = new SimpleIntegerProperty(this, "optionSignificantDigits", 6);

	private final IntegerProperty optionRandomRestarts = new SimpleIntegerProperty(this, "optionRandomRestarts", 1000);

	private final BooleanProperty optionRemoveRedundant = new SimpleBooleanProperty(this, "optionRemoveRedundant", true);

	// extend+polish cycles per restart, passed on to RazorNet
	private final IntegerProperty optionMaxPolishRounds = new SimpleIntegerProperty(this, "optionMaxPolishRounds", 1);

	private final IntegerProperty optionRandomSeed = new SimpleIntegerProperty(this, "optionRandomSeed", 0);

	private final IntegerProperty optionMaxSeconds = new SimpleIntegerProperty(this, "optionMaxSeconds", 60);

	// how sequence changes are placed on the edges: MinExcess weights each edge by its geodesic betweenness and so
	// minimizes the excess of the displayed distances; Parsimony minimizes the total number of changes
	public enum LabelingMethod {MinExcess, Parsimony}

	private final ObjectProperty<LabelingMethod> optionLabeling = new SimpleObjectProperty<>(this, "optionLabeling", LabelingMethod.MinExcess);

	// the engine's strategy, forwarded: Strict, Relaxed or Both (see RazorNet)
	private final ObjectProperty<RazorNet.Strategy> optionStrategy = new SimpleObjectProperty<>(this, "optionStrategy", RazorNet.Strategy.Both);

	// treat '?', 'N' and the IUPAC ambiguity codes as ordinary states, instead of masking their columns
	private final BooleanProperty optionMissingAsState = new SimpleBooleanProperty(this, "optionMissingAsState", false);

	// mask every column in which any sequence has a gap, as PopART does, instead of treating a gap as a state
	private final BooleanProperty optionMaskGapSites = new SimpleBooleanProperty(this, "optionMaskGapSites", false);

	private final RazorNet razorNet;

	/** prefix of the keys under which the options are saved in the program properties */
	private static final String KEY = "razornet.algorithm.RazorHaploNet";

	{
		ProgramProperties.track(KEY + optionSignificantDigits.getName(), optionSignificantDigits, 6);
		ProgramProperties.track(KEY + optionRandomRestarts.getName(), optionRandomRestarts, 1000);
		ProgramProperties.track(KEY + optionRandomSeed.getName(), optionRandomSeed, 0);
		ProgramProperties.track(KEY + optionMaxSeconds.getName(), optionMaxSeconds, 60);
		// not tracked, so they start at their defaults in every session: optionMaxPolishRounds, optionRemoveRedundant

		ProgramProperties.track(KEY + optionContractUnlabeled.getName(), optionContractUnlabeled, true);
		ProgramProperties.track(KEY + optionRemoveSuperfluous.getName(), optionRemoveSuperfluous, true);
		// not tracked: optionLabeling (an enum), optionMissingAsState and optionMaskGapSites (to be chosen anew in every session)
	}

	public RazorHaploNet() {
		// take the shared option defaults from the engine, so the two never disagree about them
		this.razorNet = new RazorNet();
		try {
			optionStrategy.set(razorNet.getOptionStrategy());
			optionRemoveRedundant.set(razorNet.getOptionRemoveRedundant());
			optionRandomRestarts.set(razorNet.getOptionRandomRestarts());
			optionMaxPolishRounds.set(razorNet.getOptionMaxPolishRounds());
			optionSignificantDigits.set(razorNet.getOptionSignificantDigits());
			optionRandomSeed.set(razorNet.getOptionRandomSeed());
			optionMaxSeconds.set(razorNet.getOptionMaxSeconds());
		} catch (Exception ignored) {
		}
	}

	@Override
	public String getShortDescription() {
		return "Computes a haplotype network from sequences: RazorNet on the Hamming distances, with mutations placed by parsimony.";
	}

	@Override
	public String getToolTip(String optionName) {
		if (optionRemoveSuperfluous.getName().equals(optionName))
			return "remove any superfluous edge for a which an alternative path of at most the same length exists";
		else if (optionContractUnlabeled.getName().equals(optionName))
			return "Contract any edge that is not labeled with any state change";
		else if (optionLabeling.getName().equals(optionName))
			return "how to place sequence changes on edges: MinExcess (faithful pairwise distances) or Parsimony (fewest changes)";
		else if (optionMissingAsState.getName().equals(optionName))
			return "treat '?', 'N' and ambiguity codes as ordinary states rather than masking out their columns; "
				   + "unknown positions then count as mutations - only for data where masking leaves too little";
		else if (optionMaskGapSites.getName().equals(optionName))
			return "mask out every column in which any sequence has a gap, as PopART does; "
				   + "otherwise a gap is a state of its own and an indel counts as a mutation";
		else return getEngine().getToolTip(optionName); // the shared options are the engine's, so are their tooltips
	}

	/**
	 * the engine that turns the Hamming distances into a network; a subclass overrides this and
	 * {@link #applyEngineOptions()} to use a different engine
	 */
	protected Distances2Network getEngine() {
		return razorNet;
	}

	/** copies the options of this algorithm to the engine, immediately before it runs */
	protected void applyEngineOptions() {
		try {
			razorNet.optionStrategyProperty().set(getOptionStrategy());
			razorNet.optionRemoveRedundantProperty().set(getOptionRemoveRedundant());
			razorNet.optionRandomRestartsProperty().set(getOptionRandomRestarts());
			razorNet.optionMaxPolishRoundsProperty().set(getOptionMaxPolishRounds());
			razorNet.optionSignificantDigitsProperty().set(getOptionSignificantDigits());
			razorNet.optionRandomSeedProperty().set(getOptionRandomSeed());
			razorNet.optionMaxSecondsProperty().set(getOptionMaxSeconds());
		} catch (Exception ignored) {
		}
	}

	@Override
	public List<String> listOptions() {
		// the basic options; in advanced mode, a separator and the advanced options follow (see AdvancedOptions)
		var list = new java.util.ArrayList<>(List.of(optionStrategy.getName(), optionLabeling.getName(), optionRemoveSuperfluous.getName(),
				optionContractUnlabeled.getName(), optionMissingAsState.getName(), optionMaskGapSites.getName(), optionMaxSeconds.getName()));
		if (razornet.algorithms.utils.AdvancedOptions.enabled()) {
			list.add(splitstree6.options.Option.SEPARATOR);
			list.addAll(List.of(optionRandomRestarts.getName(), optionRandomSeed.getName()));
		}
		return list;
	}

	public LabelingMethod getOptionLabeling() {
		return optionLabeling.get();
	}

	public ObjectProperty<LabelingMethod> optionLabelingProperty() {
		return optionLabeling;
	}

	public RazorNet.Strategy getOptionStrategy() {
		return optionStrategy.get();
	}

	public ObjectProperty<RazorNet.Strategy> optionStrategyProperty() {
		return optionStrategy;
	}

	@Override
	public void compute(ProgressListener progress, TaxaBlock taxaBlock, CharactersBlock charactersBlock, NetworkBlock networkBlock) throws IOException {
		// the active sites are taken from the alignment view of the workflow, if there is one
		compute(progress, taxaBlock, charactersBlock, networkBlock, findActiveSites());
	}

	/**
	 * runs the algorithm with the given active sites, which are used to report the positions of the mutations in
	 * the original alignment
	 *
	 * @param activeSites one-based positions of the active sites, or null to number the sites as given
	 */
	public void compute(ProgressListener progress, TaxaBlock taxaBlock, CharactersBlock charactersBlock, NetworkBlock networkBlock, BitSet activeSites) throws IOException {
		progress.setTasks("RazorHaploNet", "Initializing");
		var mapBack = mapBackIndices(activeSites);

		// Gaps and missing data: a gap is an ordinary state (an indel) and counts as a difference. A column in which
		// any sequence has '?', 'N' or an IUPAC ambiguity code is masked out for all sequences, so that all pairs are
		// compared on the same sites and the Hamming distances are a metric. With optionMissingAsState nothing is
		// masked and every character is a state of its own; with optionMaskGapSites a column in which any sequence
		// has a gap is masked out as well, as PopART does. The labeling follows the same convention.
		final int[] keptSites = keptSites(charactersBlock, !getOptionMissingAsState(), getOptionMaskGapSites());
		if (getOptionMissingAsState())
			warnMissingAsState(charactersBlock);
		var masked = (getOptionMissingAsState() ? "a gap" : getOptionMaskGapSites() ? "missing or ambiguous data or a gap" : "missing or ambiguous data");
		if (keptSites.length < charactersBlock.getNchar())
			System.err.printf("Masked %,d of %,d sites carrying %s; %,d remain%n",
					charactersBlock.getNchar() - keptSites.length, charactersBlock.getNchar(), masked, keptSites.length);
		if (keptSites.length == 0)
			throw new IOException(("All %,d sites carry %s: nothing left to compute on. "
								   + "Turn on 'Missing As State' to keep them, at the cost of counting unknown "
								   + "positions as mutations.").formatted(charactersBlock.getNchar(), masked));

		var sequences = maskedRows(charactersBlock, keptSites);   // one masked sequence per taxon, 0-based
		// site positions in the original alignment, for the edge labels
		final IntFunction<Integer> siteMapBack = k -> mapBack.apply(keptSites[k]);

		// ambiguityCharsAreStates = true: a difference is plain character inequality, so a gap counts as a state
		var utils = new NetworkSequencesAnalyzer(charactersBlock.getGapCharacter(), charactersBlock.getMissingCharacter(),
				charactersBlock.isHasAmbiguityCodes(), true);

		// sequences that are identical on the kept sites are one haplotype, represented by the first of them (1-based
		// taxon ids); the engine, the labeling and the cleanup run on the representatives, so that neither the
		// construction (TCS joins taxa in order of distance) nor the minimum-excess labeling (which weighs pairs of
		// taxa) sees a haplotype several times; the other members are attached at the end
		var groups = new LinkedHashMap<String, List<Integer>>();
		for (var t = 1; t <= taxaBlock.getNtax(); t++)
			groups.computeIfAbsent(sequences[t - 1], k -> new ArrayList<>()).add(t);
		var representatives = groups.values().stream().mapToInt(members -> members.get(0)).toArray();
		var engineTaxa = new TaxaBlock();
		engineTaxa.addTaxaByNames(Arrays.stream(representatives).mapToObj(taxaBlock::getLabel).toList());
		if (representatives.length < taxaBlock.getNtax())
			System.err.printf("Merged %,d sequences into %,d distinct haplotypes; the other sequences of a haplotype hang off its node by a zero-length edge%n",
					taxaBlock.getNtax(), representatives.length);

		var distancesBlock = new DistancesBlock();
		distancesBlock.setNtax(representatives.length);
		{
			var d = distancesBlock.getDistances();
			for (var i = 0; i < representatives.length; i++)
				for (var j = i + 1; j < representatives.length; j++)
					d[i][j] = d[j][i] = utils.differences(sequences[representatives[i] - 1], sequences[representatives[j] - 1]);
		}

		// Hamming distances on equal sites are a metric; this check is a safeguard
		razornet.algorithms.utils.MetricCheck.checkAndWarn(distancesBlock.getDistances(), "RazorHaploNet, Hamming distances");

		applyEngineOptions();

		getEngine().compute(progress, engineTaxa, distancesBlock, networkBlock);

		var network = networkBlock.getGraph();

		// the engine attached the representatives by their ids in engineTaxa: give every node its original taxon id
		for (var v : network.nodes()) {
			if (network.hasTaxa(v)) {
				var ids = new ArrayList<Integer>();
				for (var r : network.getTaxa(v))
					ids.add(representatives[r - 1]);
				network.clearTaxa(v);
				for (var t : ids)
					network.addTaxon(v, t);
			}
		}

		progress.setSubtask(getOptionLabeling() == LabelingMethod.MinExcess ? "Minimum-excess labeling" : "Parsimony labeling");

		// edge weights for the labeler: 1 for Parsimony; for MinExcess the geodesic betweenness of the edge (times
		// 10^6, plus 1 to break ties in favour of fewer changes), computed before edge removal and contraction
		ToDoubleFunction<Edge> edgeWeight;
		if (getOptionLabeling() == LabelingMethod.MinExcess) {
			var betweenness = GeodesicBetweenness.compute(network.nodes(), network.edges(),
					Edge::getSource, Edge::getTarget, network::getWeight, network::hasTaxa);
			edgeWeight = e -> 1.0e6 * betweenness.getOrDefault(e, 0) + 1.0;
		} else {
			edgeWeight = e -> 1.0;
		}

		var parsimonyLabeler = new ParsimonyLabeler<>(network.nodes(), network.edges(), Edge::getSource, Edge::getTarget, edgeWeight);
		// the labeler has to read the alphabet in the same way as the distances
		parsimonyLabeler.setMissingIsState(getOptionMissingAsState());
		var inputSequences = new HashMap<Node, String>();
		for (var v : network.nodes()) {
			if (network.hasTaxa(v)) {
				var row = network.getTaxon(v) - 1;
				// no '?', 'N' or IUPAC codes remain here unless optionMissingAsState kept their columns
				inputSequences.put(v, sequences[row]);
			}
		}

		var outputSequences = parsimonyLabeler.labelAllSites(inputSequences);
		for (var entry : outputSequences.entrySet()) {
			networkBlock.getNodeData(entry.getKey()).put(NetworkBlock.NODE_STATES_KEY, entry.getValue());
		}

		if (outputSequences.size() != network.getNumberOfNodes())
			System.err.println("Some nodes unlabeled");

		// from here on, the weight of an edge is its number of mutations
		setWeightsToMutationCounts(network, networkBlock, utils);

		if (getOptionRemoveSuperfluous()) {
			progress.setSubtask("Removing superfluous edges");
			System.err.println("----------- Superfluous Edge Removal  -----------");
			var initialEdges = network.getNumberOfEdges();

			// an edge is removed if an alternative path of the same length exists, which preserves all distances
			var edgesDeleted = 0;
			var edges = network.getEdgesAsList();
			edges.sort((a, b) -> -Double.compare(network.getWeight(a), network.getWeight(b)));
			progress.setMaximum(2L * network.getNumberOfEdges() + edges.size());
			for (var f : edges) {
				if (isSuperfluous(f.getSource(), f.getTarget(), f, network::getWeight)) {
					network.deleteEdge(f);
					edgesDeleted++;
				}
				progress.incrementProgress();
			}

			while (true) {
				var changed = false;
				var vOptional = network.nodeStream().filter(u -> u.getDegree() == 2 && !network.hasTaxa(u)).findAny();
				if (vOptional.isPresent()) {
					var v = vOptional.get();
					var e1 = v.getFirstAdjacentEdge();
					var e2 = v.getLastAdjacentEdge();
					var a = e1.getOpposite(v);
					var b = e2.getOpposite(v);
					var ab = a.getCommonEdge(b);
					if (ab == null) {
						var e = network.newEdge(a, b);
						network.setWeight(e, network.getWeight(e1) + network.getWeight(e2));
					} else {
						var w = 0.5 * (network.getWeight(ab) + network.getWeight(e1) + network.getWeight(e2));
						network.setWeight(ab, w);
					}
					network.deleteNode(v);
					changed = true;
				}
				var wOptional = network.nodeStream().filter(u -> u.getDegree() == 1 && !network.hasTaxa(u)).findAny();
				if (wOptional.isPresent()) {
					var w = wOptional.get();
					network.deleteNode(w);
					changed = true;
				}
				if (!changed)
					break;
				progress.incrementProgress();
			}
			System.err.printf("E removal: n=%4d e=%4d -> n=%4d e= %4d%n", network.getNumberOfNodes(), initialEdges, network.getNumberOfNodes(), network.getNumberOfEdges());
		}

		/* label all edges: */
		for (var e : network.edges()) {
			var sequence1 = networkBlock.getNodeData(e.getSource()).get(NetworkBlock.NODE_STATES_KEY);
			var sequence2 = networkBlock.getNodeData(e.getTarget()).get(NetworkBlock.NODE_STATES_KEY);
			if (sequence1 != null && sequence2 != null) {
				var label = utils.computeEdgeLabel(sequence1, sequence2, siteMapBack);
				networkBlock.getEdgeData(e).put(NetworkBlock.EDGE_SITES_KEY, label);
			} else {
				throw new IllegalStateException("Unlabeled node");
			}
		}

		if (getOptionContractUnlabeled()) {
			System.err.println("----------- Contracting Unlabeled Edges  -----------");
			var initialEdges = network.getNumberOfEdges();

			progress.setSubtask("Contracting unlabeled edges");
			progress.setMaximum(2L * network.getNumberOfEdges());
			progress.setProgress(0);

			for (var i = 0; i < 10000; i++) {
				var contracted = false;
				for (var e : network.edges()) {
					if (!network.hasTaxa(e.getSource()) || !network.hasTaxa(e.getTarget())) {
						var label = networkBlock.getEdgeData(e).get(NetworkBlock.EDGE_SITES_KEY);
						if (label == null || label.isBlank()) {
							var keep = (network.hasTaxa(e.getSource()) ? e.getSource() : e.getTarget());
							var other = e.getOpposite(keep);
							for (var f : other.adjacentEdges()) {
								if (f != e) {
									var u = f.getOpposite(other);
									var g = u.getCommonEdge(keep);
									if (g != null) {
										var w = 0.5 * (network.getWeight(f) + network.getWeight(g));
										network.setWeight(g, w);
									} else {
										var h = network.newEdge(u, keep);
										// the rewired edge keeps its own weight
										network.setWeight(h, network.getWeight(f));
										var sequence1 = networkBlock.getNodeData(h.getSource()).get(NetworkBlock.NODE_STATES_KEY);
										var sequence2 = networkBlock.getNodeData(h.getTarget()).get(NetworkBlock.NODE_STATES_KEY);
										if (sequence1 != null && sequence2 != null) {
											networkBlock.getEdgeData(h).put(NetworkBlock.EDGE_SITES_KEY, utils.computeEdgeLabel(sequence1, sequence2, siteMapBack));
										} else {
											throw new IllegalStateException("Unlabeled node");
										}

									}
								}
							}
							networkBlock.removeNodeData(other);
							network.deleteNode(other);
							contracted = true;
							break;
						}
					}
				}
				if (!contracted)
					break;
				progress.incrementProgress();
			}
			System.err.printf("Contracted:  n=%4d e=%4d -> n=%4d e= %4d%n", network.getNumberOfNodes(), initialEdges, network.getNumberOfNodes(), network.getNumberOfEdges());
		}

		networkBlock.setNetworkType(NetworkBlock.Type.HaplotypeNetwork);

		// contraction creates edges: set the weights again, and store them in the edge data for the viewer
		setWeightsToMutationCounts(network, networkBlock, utils);
		for (var e : network.edges())
			networkBlock.getEdgeData(e).put("weight", StringUtils.trim(network.getWeight(e)));

		warnIfExcessNegative(taxaBlock, distancesBlock, representatives, networkBlock);

		// the other sequences of each haplotype hang off its node by a zero-length edge, with the same states
		for (var members : groups.values()) {
			if (members.size() > 1) {
				var v = network.getTaxon2Node(members.get(0));
				var states = networkBlock.getNodeData(v).get(NetworkBlock.NODE_STATES_KEY);
				for (var t : members.subList(1, members.size())) {
					var w = network.newNode();
					network.addTaxon(w, t);
					network.setLabel(w, taxaBlock.getLabel(t));
					networkBlock.getNodeData(w).put(NetworkBlock.NODE_STATES_KEY, states);
					var e = network.newEdge(w, v);
					network.setWeight(e, 0);
					networkBlock.getEdgeData(e).put("weight", "0");
					networkBlock.getEdgeData(e).put(NetworkBlock.EDGE_SITES_KEY, "");
				}
			}
		}
	}

	/**
	 * Determine whether given method can be applied to given data.
	 *
	 * @param taxa  the taxa
	 * @param chars the characters matrix
	 * @return true, if method applies to given data
	 */
	@Override
	public boolean isApplicable(TaxaBlock taxa, CharactersBlock chars) {
		return taxa != null && chars != null && (chars.getDataType().isNucleotides() || chars.getNcolors() <= 8);
	}

	/**
	 * Is edge {@code e = (s,t)} superfluous, i.e. is there another path from s to t of the same length as e?
	 * Uses a Dijkstra search from s that does not use e and is bounded by the weight of e. A strictly shorter
	 * alternative path would mean that the graph does not realize its metric; the edge is then kept.
	 */
	public static boolean isSuperfluous(Node s, Node t, Edge e, ToDoubleFunction<Edge> getWeight) {
		final var EPS = 1e-9;
		final var W = getWeight.applyAsDouble(e);
		if (!e.nodes().containsAll(List.of(s, t)) || W < 0)
			throw new IllegalArgumentException("Edge e must be exactly (s,t) and have non-negative weight");

		// Dijkstra from s, bounded by W, forbidden to use e
		Map<Node, Double> dist = new HashMap<>();

		record NodeDist(Node n, double d) {
		}

		PriorityQueue<NodeDist> pq = new PriorityQueue<>(Comparator.comparingDouble(NodeDist::d));
		dist.put(s, 0.0);
		pq.add(new NodeDist(s, 0.0));

		while (!pq.isEmpty()) {
			NodeDist cur = pq.poll();
			var known = dist.get(cur.n);
			if (known == null || Math.abs(cur.d - known) > EPS) continue;   // stale
			if (cur.d > W + EPS) break;                                    // past the limit
			if (cur.n == t)
				return Math.abs(cur.d - W) <= EPS;
			for (var a : cur.n.adjacentEdges()) {
				if (a == e) continue;                                      // may not use e itself
				var nb = a.getOpposite(cur.n);
				if (nb == null) continue;
				var nd = cur.d + getWeight.applyAsDouble(a);
				if (nd > W + EPS) continue;                                // prune
				var best = dist.get(nb);
				if (best == null || nd < best - EPS) {
					dist.put(nb, nd);
					pq.add(new NodeDist(nb, nd));
				}
			}
		}
		return false;
	}


	public boolean getOptionContractUnlabeled() {
		return optionContractUnlabeled.get();
	}

	public BooleanProperty optionContractUnlabeledProperty() {
		return optionContractUnlabeled;
	}

	public boolean getOptionRemoveSuperfluous() {
		return optionRemoveSuperfluous.get();
	}

	public BooleanProperty optionRemoveSuperfluousProperty() {
		return optionRemoveSuperfluous;
	}

	public boolean getOptionMissingAsState() {
		return optionMissingAsState.get();
	}

	public BooleanProperty optionMissingAsStateProperty() {
		return optionMissingAsState;
	}

	public boolean getOptionMaskGapSites() {
		return optionMaskGapSites.get();
	}

	public BooleanProperty optionMaskGapSitesProperty() {
		return optionMaskGapSites;
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

	public boolean getOptionRemoveRedundant() {
		return optionRemoveRedundant.get();
	}

	public BooleanProperty optionRemoveRedundantProperty() {
		return optionRemoveRedundant;
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

	/** sets the weight of every edge to the number of sites at which the sequences of its two endpoints differ */
	private static void setWeightsToMutationCounts(jloda.phylo.PhyloGraph network, NetworkBlock networkBlock,
												  NetworkSequencesAnalyzer utils) {
		for (var e : network.edges()) {
			var a = networkBlock.getNodeData(e.getSource()).get(NetworkBlock.NODE_STATES_KEY);
			var b = networkBlock.getNodeData(e.getTarget()).get(NetworkBlock.NODE_STATES_KEY);
			if (a != null && b != null)
				network.setWeight(e, utils.differences(a, b));
		}
	}

	/**
	 * Reports the excess: the mutations displayed along the shortest paths minus the sequence differences, summed
	 * over all pairs of taxa. Warns if a pair of taxa is displayed with fewer mutations than its sequences differ
	 * by, which must not happen.
	 */
	private static void warnIfExcessNegative(TaxaBlock taxaBlock, DistancesBlock distancesBlock, int[] representatives, NetworkBlock networkBlock) {
		// the distances are those of the representatives, indexed by their position; the nodes carry original taxon ids
		var reducedIndex = new int[taxaBlock.getNtax() + 1];
		for (var i = 0; i < representatives.length; i++)
			reducedIndex[representatives[i]] = i + 1;
		var graph = networkBlock.getGraph();
		var nodes = new ArrayList<Node>();
		var taxonOf = new ArrayList<Integer>();
		graph.nodeStream().filter(graph::hasTaxa).forEach(v -> {
			nodes.add(v);
			taxonOf.add(graph.getTaxon(v));
		});
		var worstShortfall = 0.0;
		var worst = "";
		var totalExcess = 0.0;
		var pairs = 0;
		for (var a = 0; a < nodes.size(); a++) {
			// one shortest-path computation per taxon
			var distanceFromA = CheckPairwiseDistances.graphDistances(graph, nodes.get(a));
			for (var b = a + 1; b < nodes.size(); b++) {
				var expected = distancesBlock.get(reducedIndex[taxonOf.get(a)], reducedIndex[taxonOf.get(b)]);
				double shown = distanceFromA.get(nodes.get(b));
				totalExcess += shown - expected;
				pairs++;
				if (expected - shown > worstShortfall) {
					worstShortfall = expected - shown;
					worst = "%s vs %s: differ at %s sites, network shows only %s".formatted(
							taxaBlock.getLabel(taxonOf.get(a)), taxaBlock.getLabel(taxonOf.get(b)),
							StringUtils.trim(expected), StringUtils.trim(shown));
				}
			}
		}
		if (worstShortfall > 1e-6) {
			System.err.printf("WARNING (RazorHaploNet): the network UNDERSTATES a difference by %s (%s).%n",
					StringUtils.trim(worstShortfall), worst);
			System.err.printf("         Every differing site must change at least once on any path, so this is a bug.%n");
		} else if (pairs > 0)
			System.err.printf("Displayed mutations exceed the sequence differences by %s in total over %,d taxon pairs (the excess).%n",
					StringUtils.trim(totalExcess), pairs);
	}


	/**
	 * the 0-based columns that are kept: with {@code maskMissing}, a column in which any sequence has the missing
	 * character, 'N' or an IUPAC ambiguity code is dropped; with {@code maskGaps}, a column in which any sequence has
	 * a gap is dropped as well, as PopART does
	 */
	private static int[] keptSites(CharactersBlock charactersBlock, boolean maskMissing, boolean maskGaps) {
		var gap = Character.toUpperCase(charactersBlock.getGapCharacter());
		var missing = Character.toUpperCase(charactersBlock.getMissingCharacter());
		var kept = new ArrayList<Integer>();
		for (var j = 1; j <= charactersBlock.getNchar(); j++) {
			var ok = true;
			for (var i = 1; i <= charactersBlock.getNtax() && ok; i++) {
				var c = Character.toUpperCase(charactersBlock.get(i, j));
				if (c == gap ? maskGaps : maskMissing && (c == missing || c == 'N' || "RYSWKMBDHV".indexOf(c) >= 0))
					ok = false;
			}
			if (ok) kept.add(j - 1);
		}
		var result = new int[kept.size()];
		for (var i = 0; i < result.length; i++) result[i] = kept.get(i);
		return result;
	}

	/** warns that 'Missing As State' is on, and says what that means for the network */
	private static void warnMissingAsState(CharactersBlock charactersBlock) {
		var wouldMask = charactersBlock.getNchar() - keptSites(charactersBlock, true, false).length;
		System.err.printf("WARNING (RazorHaploNet): 'Missing As State' is on. All %,d sites are kept, including%n",
				charactersBlock.getNchar());
		System.err.printf("         the %,d that carry missing or ambiguous data and would otherwise be masked out (%.1f%%).%n",
				wouldMask, charactersBlock.getNchar() > 0 ? 100.0 * wouldMask / charactersBlock.getNchar() : 0.0);
		System.err.printf("         '?', 'N' and each ambiguity code are now ordinary states, distinct from one another%n");
		System.err.printf("         and from A/C/G/T, so an unknown position counts as a difference: 'N' against 'A' is%n");
		System.err.printf("         DRAWN AS A MUTATION although nothing is known to have changed. Sequences are pushed%n");
		System.err.printf("         apart in proportion to how much data they are missing, which on patchy alignments can%n");
		System.err.printf("         dominate the network; and many more columns now hold three or more states, so the%n");
		System.err.printf("         labeling falls to the approximate alpha-expansion solver rather than the exact min-cut.%n");
		System.err.printf("         The distances remain a metric and the network still realizes them exactly.%n");
	}

	/** One sequence per taxon (0-based), restricted to {@code keptSites}. */
	private static String[] maskedRows(CharactersBlock charactersBlock, int[] keptSites) {
		var rows = new String[charactersBlock.getNtax()];
		for (var t = 0; t < rows.length; t++) {
			var full = charactersBlock.getRow0(t);
			var buf = new char[keptSites.length];
			for (var k = 0; k < keptSites.length; k++) buf[k] = full[keptSites[k]];
			rows[t] = String.valueOf(buf);
		}
		return rows;
	}

	/**
	 * the active sites of the alignment view of the workflow that this algorithm is part of
	 *
	 * @return the active sites, or null if there is no workflow or no alignment view
	 */
	private BitSet findActiveSites() {
		if (getNode() == null || getNode().getOwner() == null)
			return null;
		var workflow = getNode().getOwner();
		var alignmentViewNode = workflow.getAlignmentViewNode();
		if (workflow.getInputDataBlock() instanceof CharactersBlock
			&& alignmentViewNode != null
			&& alignmentViewNode.getDataBlock() instanceof ViewBlock viewBlock
			&& viewBlock.getView() instanceof AlignmentView alignmentView)
			return alignmentView.getActiveSites();
		return null;
	}

	/**
	 * maps compressed (active-site) column indices back to original alignment positions
	 *
	 * @param activeSites one-based active positions, or null for the identity 1-based mapping
	 */
	private static IntFunction<Integer> mapBackIndices(BitSet activeSites) {
		if (activeSites == null)
			return i -> i + 1;
		var indexMap = new HashMap<Integer, Integer>();
		var compressed = 0;
		for (var original = activeSites.nextSetBit(1); original >= 0; original = activeSites.nextSetBit(original + 1)) {
			indexMap.put(compressed++, original);
		}
		return indexMap::get;
	}
}
