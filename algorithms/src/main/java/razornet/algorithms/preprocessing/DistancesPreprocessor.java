/*
 * DistancesPreprocessor.java Copyright (C) 2026 Daniel H. Huson
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
 * A transformation of the input distances that is applied before the RazorNet search, such as classical MDS. The
 * taxa and their order are not changed, only the distances. The implementations are listed in
 * {@link PreprocessingMethod}.
 */
public interface DistancesPreprocessor {
	/**
	 * transform the given distances
	 *
	 * @param distances symmetric input distance matrix (0-based, zero diagonal); must not be mutated
	 * @return the transformed matrix, or the given {@code distances} array itself if the transform does not apply
	 */
	double[][] apply(double[][] distances);
}
