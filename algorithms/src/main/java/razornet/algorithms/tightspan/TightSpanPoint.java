/*
 * TightSpanPoint.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.tightspan;

import java.util.Arrays;

/**
 * A point of the tight span of a metric on n taxa, given by the vector of its distances to the taxa. The point
 * of taxon i is row i of the distance matrix. The coordinates are integers, in the quantized units of
 * {@link razornet.algorithms.utils.Quantization}. Points are compared by value, so that they can be used as the
 * keys of a map; the array of a point must not be modified.
 *
 * @author Daniel Huson, 8.2026
 */
public record TightSpanPoint(int[] coordinates) {
	/**
	 * the L-infinity distance between two points; for the points of two taxa it equals their input distance
	 *
	 * @return max over i of |p[i] - q[i]|, in quantized units
	 */
	public static int distance(TightSpanPoint p, TightSpanPoint q) {
		var result = 0;
		for (var i = 0; i < p.coordinates.length; i++) {
			var diff = Math.abs(p.coordinates[i] - q.coordinates[i]);
			if (diff > result)
				result = diff;
		}
		return result;
	}

	/**
	 * the slack of the pair (i,j) at this point, {@code f[i] + f[j] - d[i][j]}: non-negative for a point of the
	 * tight span, and 0 if the pair is tight. For {@code i == j} this is {@code 2*f[i]}. Computed in long, as the sum
	 * of two coordinates need not fit into an int.
	 */
	public long slack(int i, int j, int[][] distances) {
		return (long) coordinates[i] + coordinates[j] - distances[i][j];
	}

	public int size() {
		return coordinates.length;
	}

	public int get(int i) {
		return coordinates[i];
	}

	@Override
	public boolean equals(Object other) {
		return this == other || (other instanceof TightSpanPoint that && Arrays.equals(coordinates, that.coordinates));
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(coordinates);
	}

	@Override
	public String toString() {
		return Arrays.toString(coordinates);
	}
}
