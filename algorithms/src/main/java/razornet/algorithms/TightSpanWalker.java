/*
 * TightSpanWalker.java Copyright (C) 2026 Daniel H. Huson
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

import javafx.beans.property.BooleanProperty;
import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;
import jloda.fx.util.ProgramProperties;
import jloda.graph.Node;
import jloda.phylo.PhyloGraph;
import jloda.util.CanceledException;
import jloda.util.StringUtils;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.tightspan.TightSpanConstruction;
import razornet.algorithms.utils.*;
import splitstree6.algorithms.distances.distances2network.CheckPairwiseDistances;
import splitstree6.algorithms.distances.distances2network.Distances2Network;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;

import java.util.HashMap;
import java.util.List;

/**
 * The Tight Span Walker: computes a network that realizes the distances exactly, by walking geodesics inside the
 * tight span of the metric ({@link TightSpanConstruction}). The network is then smoothed
 * ({@link NetworkSmoothing}), which preserves all distances.
 * <p>
 * Based on the implementation in PopART (Leigh &amp; Bryant 2015), modified to use exact integer arithmetic and
 * to join two points directly where no step is admissible, so that it does not fail. The tight span is due to
 * Isbell (1964) and Dress (1984).
 *
 * @author Daniel Huson, 8.2026
 */
public class TightSpanWalker extends Distances2Network {
	// number of significant digits kept when the input is scaled to integers
	private final IntegerProperty optionSignificantDigits = new SimpleIntegerProperty(this, "optionSignificantDigits", 6);

	// cap on the number of tight-span points; once it is reached, the remaining walks end in direct edges
	private final IntegerProperty optionMaxPoints = new SimpleIntegerProperty(this, "optionMaxPoints", 100000);

	// smoothing: suppress unlabeled nodes of degree at most two, and remove edges that duplicate a path that is no
	// longer; both preserve all distances
	private final BooleanProperty optionSmoothing = new SimpleBooleanProperty(this, "optionSmoothing", true);
	private final BooleanProperty optionRemoveSuperfluousEdges = new SimpleBooleanProperty(this, "optionRemoveSuperfluousEdges", true);

	private final IntegerProperty optionMaxSeconds = new SimpleIntegerProperty(this, "optionMaxSeconds", 60);

	{
		ProgramProperties.track(optionSignificantDigits, 6);
		ProgramProperties.track(optionMaxPoints, 100000);
		ProgramProperties.track(optionSmoothing, true);
		ProgramProperties.track(optionRemoveSuperfluousEdges, true);
		ProgramProperties.track(optionMaxSeconds, 60);
	}

	@Override
	public String getCitation() {
		return "Leigh & Bryant 2015; " +
			   "J.W. Leigh and D. Bryant. PopART: full-feature software for haplotype network construction. " +
			   "Methods in Ecology and Evolution, 6(9):1110-1116, 2015.;";
	}

	@Override
	public String getShortDescription() {
		return "Computes a network that realizes the given distances exactly, by walking geodesics in the tight span.";
	}

	@Override
	public List<String> listOptions() {
		return List.of(optionSignificantDigits.getName(), optionMaxPoints.getName(), optionSmoothing.getName(),
				optionRemoveSuperfluousEdges.getName(), optionMaxSeconds.getName());
	}

	@Override
	public String getToolTip(String optionName) {
		if (optionName.equals(optionSignificantDigits.getName()))
			return "number of significant digits to keep when scaling the distances to integers";
		else if (optionName.equals(optionMaxPoints.getName()))
			return "maximum number of tight-span points to create; beyond this, remaining geodesics are replaced by direct edges";
		else if (optionName.equals(optionSmoothing.getName()))
			return "smooth the network afterwards, suppressing unlabeled nodes of degree at most two (distances are preserved exactly)";
		else if (optionName.equals(optionRemoveSuperfluousEdges.getName()))
			return "remove edges whose endpoints are joined by another path of no greater length (distances are preserved exactly)";
		else if (optionName.equals(optionMaxSeconds.getName()))
			return "maximum time in seconds to run (0 = no limit)";
		else return super.getToolTip(optionName);
	}

	@Override
	public void compute(ProgressListener progress, TaxaBlock taxaBlock, DistancesBlock distancesBlock, NetworkBlock networkBlock) {
		System.err.println("Running TightSpanWalker");
		progress.setTasks("Tight Span Walker", "Initialization");

		var inputDistances = distancesBlock.getDistances();

		// a non-metric input has no realization: warn, but continue
		MetricCheck.checkAndWarn(inputDistances, "TightSpanWalker");

		// scale to even integers, so that all tight-span coordinates are integers; taxa at distance 0 are merged
		var quantization = Quantization.apply(inputDistances, getOptionSignificantDigits());

		var timeoutProgress = new TimeoutProgress(progress, getOptionMaxSeconds());

		TightSpanConstruction.Result result;
		try {
			result = TightSpanConstruction.apply(quantization.quantMatrix(), getOptionMaxPoints(), timeoutProgress);
		} catch (CanceledException ex) {
			// cancelled or out of time: a partial walk is not a realization, so there is nothing to return
			throw new RuntimeException("TightSpanWalker: cancelled or time budget exceeded before completion"
									   + " (increase maxSeconds, or reduce significantDigits)");
		}

		// smoothing preserves all distances between taxa
		progress.setSubtask("Smoothing");
		var smoothed = NetworkSmoothing.apply(result.graph(), quantization.isLabeled(), isOptionSmoothing(),
				isOptionRemoveSuperfluousEdges(), timeoutProgress);

		var network = networkBlock.getGraph();
		setupNetwork(taxaBlock, result.graph(), quantization, network);
		// a DistanceNetwork: the edge weights are lengths that realize the input distances
		networkBlock.setNetworkType(NetworkBlock.Type.DistanceNetwork);
		for (var e : network.edges()) {
			networkBlock.getEdgeData(e).put("weight", StringUtils.trim(network.getWeight(e)));
		}

		try {
			CheckPairwiseDistances.apply(network, inputDistances, 0.000001);
		} catch (RuntimeException ex) {
			System.err.println("Note: pairwise-distance check skipped (" + ex.getMessage() + ")");
		}

		if (!smoothed.isEmpty())
			System.err.printf("TightSpanWalker: smoothing removed %,d nodes and %,d edges%s%n",
					smoothed.nodesRemoved(), smoothed.edgesRemoved(),
					smoothed.stoppedEarly() ? " (stopped on the time budget --- the network still realizes the"
											  + " distances exactly, but is less tidy than it could be)" : "");
		System.err.printf("TightSpanWalker: %,d points walked, %,d nodes and %,d edges remain",
				result.numberOfPoints(), result.graph().getNumberOfNodes(), result.graph().getNumberOfEdges());
		if (result.directEdges() > 0)
			System.err.printf(", %,d of %,d geodesics ended in a direct edge", result.directEdges(),
					(quantization.size() * (quantization.size() - 1)) / 2);
		System.err.println();
		if (result.hitPointCap())
			System.err.printf("TightSpanWalker: hit the limit of %,d points --- the network realizes the distances but is "
							  + "coarser than it could be; raise maxPoints, or use fewer significant digits.%n", getOptionMaxPoints());
	}

	/**
	 * copies the internal integer graph into the network, mapping node ids back to taxa and weights back to
	 * doubles; taxa that were merged are attached to their node by edges of length zero
	 */
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

	public int getOptionSignificantDigits() {
		return optionSignificantDigits.get();
	}

	public IntegerProperty optionSignificantDigitsProperty() {
		return optionSignificantDigits;
	}

	public int getOptionMaxPoints() {
		return optionMaxPoints.get();
	}

	public IntegerProperty optionMaxPointsProperty() {
		return optionMaxPoints;
	}

	public boolean isOptionSmoothing() {
		return optionSmoothing.get();
	}

	public BooleanProperty optionSmoothingProperty() {
		return optionSmoothing;
	}

	public boolean isOptionRemoveSuperfluousEdges() {
		return optionRemoveSuperfluousEdges.get();
	}

	public BooleanProperty optionRemoveSuperfluousEdgesProperty() {
		return optionRemoveSuperfluousEdges;
	}

	public int getOptionMaxSeconds() {
		return optionMaxSeconds.get();
	}

	public IntegerProperty optionMaxSecondsProperty() {
		return optionMaxSeconds;
	}
}
