/*
 * GraphCleanAndSmooth.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.CanceledException;
import jloda.util.progress.ProgressListener;

import java.util.ArrayList;
import java.util.function.IntPredicate;

/**
 * remove superfluous nodes from graph (unlabeled nodes of degree 1 or two)
 */
public class GraphCleanAndSmooth {
	/**
	 * remove superfluous nodes from the graph
	 *
	 * @param G         graph
	 * @param isLabeled is node labeled
	 */
	public static void apply(IntGraph G, IntPredicate isLabeled, ProgressListener progress) throws CanceledException {
		var originalSize = G.getNumberOfNodes();
		while (true) {
			var size = G.getNumberOfNodes();
			var nodes = new ArrayList<>(G.nodes());
			for (var u : nodes) {
				progress.checkForCancel();
				if (!isLabeled.test(u)) {
					var d = G.getDegree(u);
					if (d <= 1) {
						G.removeNode(u);
					} else if (d == 2) {
						var edges = G.getAdjacentEdges(u);
						var e1 = edges.get(0);
						var v1 = e1.other(u);
						var e2 = edges.get(1);
						var v2 = e2.other(u);
						G.addEdge(v1, v2, G.getWeight(e1) + G.getWeight(e2));
						G.removeNode(u);
					}
				}
			}
			if (G.getNumberOfNodes() == size)
				break;
		}
		if (G.getNumberOfNodes() < originalSize) {
			G.compact();
		}
	}
}
