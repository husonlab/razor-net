/*
 * TCSHaplotypeNetwork.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.progress.ProgressListener;
import splitstree6.algorithms.distances.distances2network.Distances2Network;
import splitstree6.data.CharactersBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.data.TaxaBlock;

import java.io.IOException;
import java.util.BitSet;
import java.util.List;

/**
 * Computes a haplotype network from an alignment using {@link TCS}: {@link RazorHaploNet} with TCS as the
 * engine. TCS does not realize the Hamming distances, so the mutations displayed between two taxa may be fewer
 * than the number of sites at which their sequences differ.
 *
 * @author Daniel Huson, 8.2026
 */
public class TCSHaplotypeNetwork extends RazorHaploNet {
	private final TCS tcs;

	public TCSHaplotypeNetwork() {
		// the shared option defaults are taken from the TCS engine
		this.tcs = new TCS();
		try {
			optionMaxSecondsProperty().set(tcs.getOptionMaxSeconds());
		} catch (Exception ignored) {
		}
	}

	@Override
	protected Distances2Network getEngine() {
		return tcs;
	}

	@Override
	protected void applyEngineOptions() {
		try {
			tcs.optionMaxSecondsProperty().set(getOptionMaxSeconds());
		} catch (Exception ignored) {
		}
	}

	@Override
	public String getCitation() {
		return tcs.getCitation();
	}

	@Override
	public String getShortDescription() {
		return "Computes a haplotype network from sequences using TCS (statistical parsimony).";
	}

	@Override
	public void compute(ProgressListener progress, TaxaBlock taxaBlock, CharactersBlock charactersBlock,
						NetworkBlock networkBlock, BitSet activeSites) throws IOException {
		warnNotARealization();
		super.compute(progress, taxaBlock, charactersBlock, networkBlock, activeSites);
	}

	/** notes that TCS does not realize the Hamming distances, so that the reported excess can be negative */
	private static void warnNotARealization() {
		System.err.printf("NOTE (TCSHaplotypeNetwork): TCS does not realize the Hamming distances --- it joins%n");
		System.err.printf("      haplotypes in order of increasing distance and only forbids the worst shortcuts. So the%n");
		System.err.printf("      mutations drawn along a geodesic need not number at least the number of sites at which its%n");
		System.err.printf("      two taxa differ: the reported excess can go negative and the shortfall need not be zero.%n");
		System.err.printf("      That is a fact about the network, not a fault in the labeling.%n");
	}

	/** the options of the parent that apply to TCS */
	@Override
	public List<String> listOptions() {
		return List.of(optionRemoveSuperfluousProperty().getName(), optionContractUnlabeledProperty().getName(),
				optionLabelingProperty().getName(), optionMissingAsStateProperty().getName(),
				optionMaxSecondsProperty().getName());
	}
}
