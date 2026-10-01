/*
 * GeodesicBetweenness.java Copyright (C) 2026 Daniel H. Huson
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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;

/**
 * Geodesic betweenness of the edges of a weighted graph: for each edge, the number of pairs of taxa whose
 * shortest path runs through it. One shortest path is used per pair. Used as the edge weight of the
 * {@link ParsimonyLabeler} for the minimum-excess labeling.
 *
 * @author Daniel Huson, 2026
 */
public class GeodesicBetweenness {

	public static <N, E> Map<E, Integer> compute(Iterable<N> nodes, Iterable<E> edges,
												 Function<E, N> source, Function<E, N> target,
												 ToDoubleFunction<E> weight, Predicate<N> isTaxon) {
		// incident edges per node
		var incident = new HashMap<N, List<E>>();
		for (var n : nodes) incident.put(n, new ArrayList<>());
		for (var e : edges) {
			incident.get(source.apply(e)).add(e);
			incident.get(target.apply(e)).add(e);
		}
		var taxa = new ArrayList<N>();
		for (var n : nodes) if (isTaxon.test(n)) taxa.add(n);

		var betweenness = new HashMap<E, Integer>();
		for (var e : edges) betweenness.put(e, 0);

		for (var i = 0; i < taxa.size(); i++) {
			var s = taxa.get(i);
			var prevEdge = new HashMap<N, E>();
			dijkstra(s, incident, source, target, weight, prevEdge);
			for (var j = i + 1; j < taxa.size(); j++) {
				var t = taxa.get(j);
				// walk the geodesic t -> s via predecessor edges, incrementing betweenness
				var cur = t;
				while (!cur.equals(s)) {
					var e = prevEdge.get(cur);
					if (e == null) break; // t unreachable from s (disconnected); skip
					betweenness.merge(e, 1, Integer::sum);
					cur = source.apply(e).equals(cur) ? target.apply(e) : source.apply(e);
				}
			}
		}
		return betweenness;
	}

	private static <N, E> void dijkstra(N src, Map<N, List<E>> incident, Function<E, N> source, Function<E, N> target,
										ToDoubleFunction<E> weight, Map<N, E> prevEdge) {
		var dist = new HashMap<N, Double>();
		dist.put(src, 0.0);
		var pq = new PriorityQueue<N>((a, b) -> Double.compare(dist.getOrDefault(a, Double.MAX_VALUE), dist.getOrDefault(b, Double.MAX_VALUE)));
		pq.add(src);
		var done = new HashMap<N, Boolean>();
		while (!pq.isEmpty()) {
			var u = pq.poll();
			if (done.put(u, Boolean.TRUE) != null) continue;
			var du = dist.get(u);
			for (var e : incident.get(u)) {
				var w = source.apply(e).equals(u) ? target.apply(e) : source.apply(e);
				var nd = du + weight.applyAsDouble(e);
				if (nd < dist.getOrDefault(w, Double.MAX_VALUE)) {
					dist.put(w, nd);
					prevEdge.put(w, e);
					pq.add(w);
				}
			}
		}
	}
}
