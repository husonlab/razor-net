/*
 * TightSpanHaplotypeNetwork.java Copyright (C) 2026 Daniel H. Huson
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
import splitstree6.algorithms.distances.distances2network.Distances2Network;

import java.util.List;

/**
 * Computes a haplotype network from an alignment using the {@link TightSpanWalker}: {@link RazorHaploNet} with
 * the Tight Span Walker as the engine, which realizes the Hamming distances exactly.
 *
 * @author Daniel Huson, 8.2026
 */
public class TightSpanHaplotypeNetwork extends RazorHaploNet {
	private final TightSpanWalker tightSpanWalker;

	// the engine's cap on the number of tight-span points, offered as an option of this algorithm
	private final IntegerProperty optionMaxPoints = new SimpleIntegerProperty(this, "optionMaxPoints", 100000);

	{
		ProgramProperties.track(optionMaxPoints, 100000);
	}

	public TightSpanHaplotypeNetwork() {
		// the shared option defaults are taken from the Tight Span Walker engine
		this.tightSpanWalker = new TightSpanWalker();
		try {
			optionSignificantDigitsProperty().set(tightSpanWalker.getOptionSignificantDigits());
			optionMaxSecondsProperty().set(tightSpanWalker.getOptionMaxSeconds());
			optionMaxPoints.set(tightSpanWalker.getOptionMaxPoints());
		} catch (Exception ignored) {
		}
	}

	@Override
	protected Distances2Network getEngine() {
		return tightSpanWalker;
	}

	@Override
	protected void applyEngineOptions() {
		try {
			tightSpanWalker.optionSignificantDigitsProperty().set(getOptionSignificantDigits());
			tightSpanWalker.optionMaxSecondsProperty().set(getOptionMaxSeconds());
			tightSpanWalker.optionMaxPointsProperty().set(getOptionMaxPoints());
		} catch (Exception ignored) {
		}
	}

	@Override
	public String getCitation() {
		return tightSpanWalker.getCitation();
	}

	@Override
	public String getShortDescription() {
		return "Computes a haplotype network that realizes the Hamming distances exactly, using the Tight Span Walker.";
	}

	/** the options of the parent that apply to the Tight Span Walker, and its cap on the number of points */
	@Override
	public List<String> listOptions() {
		return List.of(optionRemoveSuperfluousProperty().getName(), optionContractUnlabeledProperty().getName(),
				optionLabelingProperty().getName(), optionMissingAsStateProperty().getName(),
				optionMaxPoints.getName(), optionMaxSecondsProperty().getName());
	}

	@Override
	public String getToolTip(String optionName) {
		if (optionMaxPoints.getName().equals(optionName))
			return "maximum number of tight-span points to create; beyond this, remaining geodesics are replaced by direct edges";
		else return super.getToolTip(optionName);
	}

	public int getOptionMaxPoints() {
		return optionMaxPoints.get();
	}

	public IntegerProperty optionMaxPointsProperty() {
		return optionMaxPoints;
	}
}
