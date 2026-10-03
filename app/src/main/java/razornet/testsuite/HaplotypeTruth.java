/*
 * HaplotypeTruth.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.testsuite;

import java.io.IOException;
import java.util.*;
import java.util.function.UnaryOperator;

/**
 * The true haplotype tree of a sample: the genealogy with the sequences of all its nodes, in which every edge
 * without a mutation has been contracted and every unsampled node with two neighbours has been removed. Its nodes
 * are the sampled haplotypes and the unsampled ancestors at which lineages part; the length of an edge is the
 * number of sites at which the sequences along it differ.
 *
 * @author Daniel Huson, 2026
 */
public final class HaplotypeTruth {
	private final TaxonGraph graph;
	private final Set<String> edges;
	private final Set<String> unplaced = new TreeSet<>();

	/**
	 * @param newick     the genealogy, all nodes named
	 * @param sequences  the sequence of every node of the genealogy, by name
	 * @param mask       restricts a sequence to the sites in use
	 * @param haplotypes the sampled haplotypes: label and sequence, restricted to the sites in use
	 */
	public HaplotypeTruth(String newick, Map<String, String> sequences, UnaryOperator<String> mask, Map<String, String> haplotypes) throws IOException {
		var tree = NewickTree.parse(newick);
		var n = tree.size();
		var sequence = new String[n];
		for (var v = 0; v < n; v++) {
			var full = sequences.get(tree.label(v));
			if (full == null)
				throw new IOException("No sequence for node " + tree.label(v));
			sequence[v] = mask.apply(full);
		}
		// contract the edges without a mutation
		var group = new int[n];
		for (var v = 0; v < n; v++)
			group[v] = v;
		for (var v = 1; v < n; v++)
			if (sequence[v].equals(sequence[tree.parent(v)]))
				group[find(group, v)] = find(group, tree.parent(v));
		var sampled = new boolean[n];
		for (var v = 0; v < n; v++)
			if (tree.isLeaf(v))
				sampled[find(group, v)] = true;
		var neighbors = new HashMap<Integer, Map<Integer, Double>>();
		for (var v = 0; v < n; v++)
			if (find(group, v) == v)
				neighbors.put(v, new HashMap<>());
		for (var v = 1; v < n; v++) {
			var a = find(group, v);
			var b = find(group, tree.parent(v));
			if (a != b) {
				var length = (double) TaxonGraph.differences(sequence[a], sequence[b]);
				neighbors.get(a).put(b, length);
				neighbors.get(b).put(a, length);
			}
		}
		removeUnsampledOfDegreeTwo(neighbors, v -> sampled[v]);

		var builder = new TaxonGraph.Builder();
		var index = new HashMap<Integer, Integer>();
		for (var v : neighbors.keySet()) {
			index.put(v, builder.addNode());
			builder.setSequence(index.get(v), sequence[v]);
		}
		for (var v : neighbors.keySet())
			for (var entry : neighbors.get(v).entrySet())
				if (v < entry.getKey())
					builder.addEdge(index.get(v), index.get(entry.getKey()), entry.getValue());
		// a haplotype is placed at the node that carries its sequence, if there is exactly one such sampled node
		for (var haplotype : haplotypes.entrySet()) {
			var nodes = new ArrayList<Integer>();
			for (var v : neighbors.keySet())
				if (sampled[v] && sequence[v].equals(haplotype.getValue()))
					nodes.add(v);
			if (nodes.size() == 1)
				builder.setTaxon(haplotype.getKey(), index.get(nodes.get(0)));
			else
				unplaced.add(haplotype.getKey());
		}
		graph = builder.build();
		edges = edgesBySequence(graph);
	}

	private static int find(int[] group, int v) {
		while (group[v] != v)
			v = group[v] = group[group[v]];
		return v;
	}

	private static void removeUnsampledOfDegreeTwo(Map<Integer, Map<Integer, Double>> neighbors, java.util.function.IntPredicate sampled) {
		var changed = true;
		while (changed) {
			changed = false;
			for (var v : new ArrayList<>(neighbors.keySet())) {
				var adjacent = neighbors.get(v);
				if (!sampled.test(v) && adjacent.size() == 2) {
					var iterator = adjacent.entrySet().iterator();
					var first = iterator.next();
					var second = iterator.next();
					int a = first.getKey(), b = second.getKey();
					if (neighbors.get(a).containsKey(b))
						continue; // its two neighbours are joined already: the node stays
					var length = first.getValue() + second.getValue();
					neighbors.get(a).remove(v);
					neighbors.get(b).remove(v);
					neighbors.get(a).put(b, length);
					neighbors.get(b).put(a, length);
					neighbors.remove(v);
					changed = true;
				}
			}
		}
	}

	/**
	 * the edges of a haplotype network as pairs of sequences, after the nodes without a taxon that have two
	 * neighbours have been removed; an edge at a node without a sequence is given a key that matches nothing
	 */
	public static Set<String> edgesBySequence(TaxonGraph graph) {
		var neighbors = new HashMap<Integer, Map<Integer, Double>>();
		for (var v = 0; v < graph.getNumberOfNodes(); v++)
			neighbors.put(v, new HashMap<>());
		for (var e = 0; e < graph.getNumberOfEdges(); e++) {
			neighbors.get(graph.getSource(e)).put(graph.getTarget(e), graph.getWeight(e));
			neighbors.get(graph.getTarget(e)).put(graph.getSource(e), graph.getWeight(e));
		}
		removeUnsampledOfDegreeTwo(neighbors, graph::hasTaxon);
		var result = new TreeSet<String>();
		for (var v : neighbors.keySet()) {
			for (var w : neighbors.get(v).keySet()) {
				if (v < w) {
					var a = graph.getSequence(v) != null ? graph.getSequence(v).toUpperCase() : "?" + v;
					var b = graph.getSequence(w) != null ? graph.getSequence(w).toUpperCase() : "?" + w;
					result.add(a.compareTo(b) <= 0 ? a + "|" + b : b + "|" + a);
				}
			}
		}
		return result;
	}

	/** the tree; its taxa are the haplotypes that could be placed */
	public TaxonGraph getGraph() {
		return graph;
	}

	/** the edges of the tree, as pairs of sequences */
	public Set<String> getEdges() {
		return edges;
	}

	/** the haplotypes whose sequence arose more than once in the genealogy, so that they have no single node */
	public Set<String> getUnplaced() {
		return unplaced;
	}
}
