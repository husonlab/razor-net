/*
 * SlackStructure.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.utils;

/**
 * Counts the interior taxa of a metric: a taxon is interior if it lies between two other taxa,
 * {@code d(x,y) + d(x,z) = d(y,z)}. Such a taxon has slack zero. RazorNet uses the count to decide which strategy
 * runs first.
 *
 * @author Daniel Huson, 2026
 */
public class SlackStructure {
	/**
	 * counts the points that lie between two other points
	 *
	 * @param d a metric as a symmetric matrix with zero diagonal (distinct points, as after quantization)
	 * @return the number of interior points; O(n^3), with an early exit per point
	 */
	public static int countInterior(int[][] d) {
		var n = d.length;
		var count = 0;
		for (var x = 0; x < n; x++) {
			var row = d[x];
			search:
			for (var y = 0; y < n; y++) {
				if (y == x) continue;
				var dxy = row[y];
				var rowY = d[y];
				for (var z = y + 1; z < n; z++) {
					if (z != x && dxy + row[z] == rowY[z]) {
						count++;
						break search;
					}
				}
			}
		}
		return count;
	}
}
