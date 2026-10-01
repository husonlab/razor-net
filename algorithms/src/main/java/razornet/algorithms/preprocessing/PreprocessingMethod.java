/*
 * PreprocessingMethod.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.preprocessing;

/**
 * The distance preprocessing methods that RazorNet offers. To add one, add a constant here and return its
 * {@link DistancesPreprocessor} from {@link #newPreprocessor(double)}.
 */
public enum PreprocessingMethod {
	/** no preprocessing: use the input distances as given (the default) */
	None,
	/** classical MDS (Torgerson/Gower) to 2 dimensions, then use the resulting planar Euclidean distances */
	ClassicalMDS2D,
	/** classical MDS (Torgerson/Gower) to 3 dimensions, then use the resulting Euclidean distances */
	ClassicalMDS3D,
	/** split pruning: remove the splits of smallest weight from the split decomposition, keeping the residue */
	PruneSmallSplits;

	/**
	 * @param splitRemoveFraction fraction of the total split weight to remove, for {@link #PruneSmallSplits}
	 *                            (ignored by the other methods)
	 * @return a preprocessor implementing this method (the identity transform for {@link #None})
	 */
	public DistancesPreprocessor newPreprocessor(double splitRemoveFraction) {
		return switch (this) {
			case None -> distances -> distances;
			case ClassicalMDS2D -> new ClassicalMdsPreprocessor(2);
			case ClassicalMDS3D -> new ClassicalMdsPreprocessor(3);
			case PruneSmallSplits -> new SplitDecompositionPreprocessor(splitRemoveFraction);
		};
	}
}
