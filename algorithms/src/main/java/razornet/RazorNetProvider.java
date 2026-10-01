/*
 * RazorNetProvider.java Copyright (C) 2026 Daniel H. Huson
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

package razornet;

import razornet.algorithms.*;
import splitstree6.algorithms.AlgorithmProvider;
import splitstree6.data.CharactersBlock;
import splitstree6.data.DistancesBlock;
import splitstree6.data.NetworkBlock;
import splitstree6.workflow.Algorithm;

import java.util.List;
import java.util.Map;

/**
 * Provides the RazorNet algorithms to SplitsTree6 and to the RazorNet program: {@link RazorNet} and
 * {@link RazorHaploNet} and, for comparison, the Tight Span Walker and TCS with their haplotype variants.
 *
 * @author Daniel Huson and Momoko Hayamizu, 2026
 */
public class RazorNetProvider implements AlgorithmProvider {
	@Override
	public List<Algorithm<?, ?>> getAlgorithms() {
		return List.of(createRazorNet(Map.of()), createRazorHaploNet(Map.of()),
				createTightSpanWalker(Map.of()), createTightSpanHaplotypeNetwork(Map.of()),
				createTCS(Map.of()), createTCSHaplotypeNetwork(Map.of()));
	}

	@Override
	public boolean provides(String simpleClassName) {
		return simpleClassName.equals("RazorNet") || simpleClassName.equals("RazorHaploNet")
			   || simpleClassName.equals("TightSpanWalker") || simpleClassName.equals("TightSpanHaplotypeNetwork")
			   || simpleClassName.equals("TCS") || simpleClassName.equals("TCSHaplotypeNetwork");
	}

	@Override
	public Algorithm<?, ?> createAlgorithm(String simpleClassName, Map<String, String> parameters) {
		return switch (simpleClassName) {
			case "RazorNet" -> createRazorNet(parameters);
			case "RazorHaploNet" -> createRazorHaploNet(parameters);
			case "TightSpanWalker" -> createTightSpanWalker(parameters);
			case "TightSpanHaplotypeNetwork" -> createTightSpanHaplotypeNetwork(parameters);
			case "TCS" -> createTCS(parameters);
			case "TCSHaplotypeNetwork" -> createTCSHaplotypeNetwork(parameters);
			default -> null;
		};
	}

	/** RazorNet: distances to a minimum-length realization, with the Strict, Relaxed or both strategies */
	public Algorithm<DistancesBlock, NetworkBlock> createRazorNet(Map<String, String> parameters) {
		return new RazorNet();
	}

	/** the haplotype network on RazorNet */
	public Algorithm<CharactersBlock, NetworkBlock> createRazorHaploNet(Map<String, String> parameters) {
		return new RazorHaploNet();
	}

	/** the Tight Span Walker, after PopART (Leigh &amp; Bryant 2015) */
	public Algorithm<DistancesBlock, NetworkBlock> createTightSpanWalker(Map<String, String> parameters) {
		return new TightSpanWalker();
	}

	/** the haplotype network on the Tight Span Walker */
	public Algorithm<CharactersBlock, NetworkBlock> createTightSpanHaplotypeNetwork(Map<String, String> parameters) {
		return new TightSpanHaplotypeNetwork();
	}

	/** TCS (Templeton, Crandall &amp; Sing 1992), after PopART (Leigh &amp; Bryant 2015) */
	public Algorithm<DistancesBlock, NetworkBlock> createTCS(Map<String, String> parameters) {
		return new TCS();
	}

	/** the haplotype network on TCS */
	public Algorithm<CharactersBlock, NetworkBlock> createTCSHaplotypeNetwork(Map<String, String> parameters) {
		return new TCSHaplotypeNetwork();
	}
}
