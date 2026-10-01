/*
 * TCS.java Copyright (C) 2026 Daniel H. Huson
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

import javafx.beans.property.IntegerProperty;
import javafx.beans.property.SimpleIntegerProperty;
import jloda.fx.util.ProgramProperties;
import jloda.graph.Node;
import jloda.phylo.PhyloGraph;
import jloda.util.CanceledException;
import jloda.util.StringUtils;
import jloda.util.progress.ProgressListener;
import razornet.algorithms.tcs.TCSConstruction;
import razornet.algorithms.utils.NetworkSmoothing;
import razornet.algorithms.utils.TimeoutProgress;
import splitstree6.algorithms.distances.distances2network.Distances2Network;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * TCS: the haplotype network method of Templeton, Crandall and Sing (1992). Components are joined in order of
 * increasing distance by paths of unit steps through intermediate vertices ({@link TCSConstruction}); unlabeled
 * vertices of degree at most two are then suppressed.
 * <p>
 * Based on the implementation in PopART (Leigh &amp; Bryant 2015), modified so that it does not fail when two
 * components are already closer than the distance being merged. The input is a distance matrix that counts
 * mutations; other distances are rounded to integers. The network need not realize the distances. The 95%
 * parsimony connection limit is not implemented.
 *
 * @author Daniel Huson, 8.2026
 */
public class TCS extends Distances2Network {
	private final IntegerProperty optionMaxSeconds = new SimpleIntegerProperty(this, "optionMaxSeconds", 60);

	{
		ProgramProperties.track(optionMaxSeconds, 60);
	}

	@Override
	public String getCitation() {
		return "Templeton et al 1992; " +
			   "A.R. Templeton, K.A. Crandall and C.F. Sing. A cladistic analysis of phenotypic associations with " +
			   "haplotypes inferred from restriction endonuclease mapping and DNA sequence data. III. Cladogram " +
			   "estimation. Genetics, 132(2):619-633, 1992.;" +
			   "Leigh & Bryant 2015; " +
			   "J.W. Leigh and D. Bryant. PopART: full-feature software for haplotype network construction. " +
			   "Methods in Ecology and Evolution, 6(9):1110-1116, 2015.;";
	}

	@Override
	public String getShortDescription() {
		return "Computes a haplotype network by joining haplotypes in order of increasing distance (statistical parsimony).";
	}

	@Override
	public List<String> listOptions() {
		return List.of(optionMaxSeconds.getName());
	}

	@Override
	public String getToolTip(String optionName) {
		if (optionName.equals(optionMaxSeconds.getName()))
			return "maximum time in seconds to run (0 = no limit)";
		else return super.getToolTip(optionName);
	}

	@Override
	public void compute(ProgressListener progress, TaxaBlock taxaBlock, DistancesBlock distancesBlock, NetworkBlock networkBlock) {
		System.err.println("Running TCS");
		progress.setTasks("TCS", "Initialization");

		var input = distancesBlock.getDistances();
		var ntax = input.length;

		// the construction lays paths of unit steps: round the distances to integers
		var rounded = new int[ntax][ntax];
		var maxDeviation = 0.0;
		for (var i = 0; i < ntax; i++) {
			for (var j = 0; j < ntax; j++) {
				rounded[i][j] = (int) Math.round(input[i][j]);
				maxDeviation = Math.max(maxDeviation, Math.abs(input[i][j] - rounded[i][j]));
			}
		}
		if (maxDeviation > 1e-9)
			System.err.printf("TCS: the distances are not whole numbers (largest rounding %.3g) --- they are rounded, "
							  + "but then an edge no longer counts mutations.%n", maxDeviation);

		// taxa at distance 0 are one haplotype: collapsed here, expanded again in setupNetwork
		var groups = collapseDuplicates(rounded);
		var reduced = new int[groups.size()][groups.size()];
		for (var a = 0; a < groups.size(); a++)
			for (var b = 0; b < groups.size(); b++)
				reduced[a][b] = rounded[groups.get(a).get(0)][groups.get(b).get(0)];
		if (groups.size() < ntax)
			System.err.printf("TCS: %,d taxa collapse to %,d distinct haplotypes%n", ntax, groups.size());

		var timeoutProgress = new TimeoutProgress(progress, getOptionMaxSeconds());

		TCSConstruction.Result result;
		try {
			result = TCSConstruction.apply(reduced, timeoutProgress);
		} catch (CanceledException ex) {
			throw new RuntimeException("TCS: cancelled or time budget exceeded before completion (increase maxSeconds)");
		}

		// only unlabeled vertices of degree at most two are suppressed; no edges are removed
		progress.setSubtask("Smoothing");
		var smoothed = NetworkSmoothing.apply(result.graph(), id -> id < groups.size(), true, false, timeoutProgress);

		var network = networkBlock.getGraph();
		setupNetwork(taxaBlock, result.graph(), groups, network);
		// a DistanceNetwork; when run on sequences, RazorHaploNet sets the type to HaplotypeNetwork
		networkBlock.setNetworkType(NetworkBlock.Type.DistanceNetwork);
		for (var e : network.edges()) {
			networkBlock.getEdgeData(e).put("weight", StringUtils.trim(network.getWeight(e)));
		}

		System.err.printf("TCS: %,d haplotypes, %,d intermediate vertices (%,d after smoothing), %,d edges%n",
				groups.size(), result.intermediates(), result.intermediates() - smoothed.nodesRemoved(),
				network.getNumberOfEdges());
		if (result.joinsSkipped() > 0)
			System.err.printf("TCS: %,d joins added nothing, the two components already being at least that close "
							  + "(PopART raises an error here)%n", result.joinsSkipped());
	}

	/**
	 * group taxa that are at distance 0 from one another; the first member of each group represents it
	 */
	private static List<List<Integer>> collapseDuplicates(int[][] distances) {
		var groups = new ArrayList<List<Integer>>();
		var assigned = new boolean[distances.length];
		for (var i = 0; i < distances.length; i++) {
			if (assigned[i])
				continue;
			var group = new ArrayList<Integer>();
			group.add(i);
			assigned[i] = true;
			for (var j = i + 1; j < distances.length; j++) {
				if (!assigned[j] && distances[i][j] == 0) {
					group.add(j);
					assigned[j] = true;
				}
			}
			groups.add(group);
		}
		return groups;
	}

	/**
	 * copy the internal graph into the network. Node ids below {@code groups.size()} are haplotypes; a
	 * haplotype shared by several taxa is expanded into one zero-length pendant per taxon, as RazorNet does.
	 */
	private void setupNetwork(TaxaBlock taxaBlock, razornet.algorithms.utils.IntGraph intGraph,
							  List<List<Integer>> groups, PhyloGraph network) {
		var nodeMap = new HashMap<Integer, Node>();
		for (var id : intGraph.nodes()) {
			var v = network.newNode();
			nodeMap.put(id, v);
			if (id < groups.size()) {
				var taxa = groups.get(id);
				if (taxa.size() == 1) {
					var t = taxa.get(0) + 1;
					network.addTaxon(v, t);
					network.setLabel(v, taxaBlock.getLabel(t));
				} else {
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
		for (var e : intGraph.edges()) {
			var f = network.newEdge(nodeMap.get(e.u()), nodeMap.get(e.v()));
			// use getEdgeWeights() here: it sets up the network's edge weights
			network.getEdgeWeights().put(f, (double) intGraph.getWeight(e));
		}
	}

	public int getOptionMaxSeconds() {
		return optionMaxSeconds.get();
	}

	public IntegerProperty optionMaxSecondsProperty() {
		return optionMaxSeconds;
	}
}
