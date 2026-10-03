/*
 * TaxonGraph.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.graph.Node;
import jloda.phylo.PhyloGraph;

import java.io.IOException;
import java.util.*;
import java.util.function.IntFunction;
import java.util.function.ToDoubleFunction;

/**
 * A weighted undirected graph with taxa on some of its nodes, held in arrays: the form in which the comparison
 * tools measure a network, a split network or a tree. Taxa are identified by their labels.
 *
 * @author Daniel Huson, 2026
 */
public final class TaxonGraph {
	private final int numberOfNodes;
	private final int[] source;
	private final int[] target;
	private final double[] weight;
	private final Map<String, Integer> taxon2node;
	private final String[] sequence;     // node -> its sequence, or null
	private final int[][] adjacentEdges; // node -> indices of its edges

	private TaxonGraph(int numberOfNodes, int[] source, int[] target, double[] weight, Map<String, Integer> taxon2node, String[] sequence) {
		this.numberOfNodes = numberOfNodes;
		this.source = source;
		this.target = target;
		this.weight = weight;
		this.taxon2node = taxon2node;
		this.sequence = sequence;
		var degree = new int[numberOfNodes];
		for (var e = 0; e < source.length; e++) {
			degree[source[e]]++;
			degree[target[e]]++;
		}
		adjacentEdges = new int[numberOfNodes][];
		for (var v = 0; v < numberOfNodes; v++)
			adjacentEdges[v] = new int[degree[v]];
		Arrays.fill(degree, 0);
		for (var e = 0; e < source.length; e++) {
			adjacentEdges[source[e]][degree[source[e]]++] = e;
			adjacentEdges[target[e]][degree[target[e]]++] = e;
		}
	}

	/** collects nodes, edges and taxa, then builds the graph */
	public static final class Builder {
		private int nodes = 0;
		private final ArrayList<int[]> edges = new ArrayList<>();
		private final ArrayList<Double> weights = new ArrayList<>();
		private final Map<String, Integer> taxon2node = new LinkedHashMap<>();
		private final Map<Integer, String> sequences = new HashMap<>();

		public int addNode() {
			return nodes++;
		}

		public void setSequence(int node, String sequence) {
			sequences.put(node, sequence);
		}

		public void addEdge(int u, int v, double weight) {
			edges.add(new int[]{u, v});
			weights.add(weight);
		}

		public void setTaxon(String label, int node) {
			taxon2node.put(label, node);
		}

		public TaxonGraph build() {
			var source = new int[edges.size()];
			var target = new int[edges.size()];
			var weight = new double[edges.size()];
			for (var e = 0; e < edges.size(); e++) {
				source[e] = edges.get(e)[0];
				target[e] = edges.get(e)[1];
				weight[e] = weights.get(e);
			}
			var sequence = new String[nodes];
			for (var entry : sequences.entrySet())
				sequence[entry.getKey()] = entry.getValue();
			return new TaxonGraph(nodes, source, target, weight, taxon2node, sequence);
		}
	}

	/**
	 * the graph of a haplotype network: a node has the sequence given for it, and the length of an edge is the
	 * number of sites at which the sequences of its two nodes differ, or its weight if a sequence is missing
	 */
	public static TaxonGraph ofHaplotypeNetwork(PhyloGraph graph, IntFunction<String> taxonLabel, java.util.function.Function<Node, String> nodeSequence) {
		var builder = new Builder();
		var node2index = new HashMap<Node, Integer>();
		for (var v : graph.nodes()) {
			var index = builder.addNode();
			node2index.put(v, index);
			if (nodeSequence.apply(v) != null)
				builder.setSequence(index, nodeSequence.apply(v));
			if (graph.hasTaxa(v))
				for (var t : graph.getTaxa(v))
					builder.setTaxon(taxonLabel.apply(t), index);
		}
		for (var e : graph.edges()) {
			var a = nodeSequence.apply(e.getSource());
			var b = nodeSequence.apply(e.getTarget());
			builder.addEdge(node2index.get(e.getSource()), node2index.get(e.getTarget()),
					a != null && b != null && a.length() == b.length() ? differences(a, b) : graph.getWeight(e));
		}
		return builder.build();
	}

	/** the number of positions at which two sequences of the same length differ, ignoring case */
	public static int differences(String a, String b) {
		var count = 0;
		for (var i = 0; i < a.length(); i++)
			if (Character.toUpperCase(a.charAt(i)) != Character.toUpperCase(b.charAt(i)))
				count++;
		return count;
	}

	/**
	 * the graph of a network or split network
	 *
	 * @param taxonLabel the label of a taxon id
	 * @param edgeWeight the length of an edge
	 */
	public static TaxonGraph of(PhyloGraph graph, IntFunction<String> taxonLabel, ToDoubleFunction<jloda.graph.Edge> edgeWeight) {
		var builder = new Builder();
		var node2index = new HashMap<Node, Integer>();
		for (var v : graph.nodes()) {
			var index = builder.addNode();
			node2index.put(v, index);
			if (graph.hasTaxa(v))
				for (var t : graph.getTaxa(v))
					builder.setTaxon(taxonLabel.apply(t), index);
		}
		for (var e : graph.edges())
			builder.addEdge(node2index.get(e.getSource()), node2index.get(e.getTarget()), edgeWeight.applyAsDouble(e));
		return builder.build();
	}

	/** a tree in Newick format; its labeled leaves are the taxa, and edges without a length have length 1 */
	public static TaxonGraph ofNewick(String newick) throws IOException {
		var tree = NewickTree.parse(newick);
		var builder = new Builder();
		for (var v = 0; v < tree.size(); v++) {
			builder.addNode();
			if (tree.isLeaf(v) && tree.label(v) != null)
				builder.setTaxon(tree.label(v), v);
		}
		for (var v = 0; v < tree.size(); v++)
			if (tree.parent(v) >= 0)
				builder.addEdge(tree.parent(v), v, Double.isNaN(tree.length(v)) ? 1.0 : tree.length(v));
		return builder.build();
	}

	public int getNumberOfNodes() {
		return numberOfNodes;
	}

	public int getNumberOfEdges() {
		return source.length;
	}

	public int getSource(int edge) {
		return source[edge];
	}

	public int getTarget(int edge) {
		return target[edge];
	}

	public double getWeight(int edge) {
		return weight[edge];
	}

	public int getDegree(int node) {
		return adjacentEdges[node].length;
	}

	/** the sequence of a node, or null */
	public String getSequence(int node) {
		return sequence[node];
	}

	/** is some taxon on this node? */
	public boolean hasTaxon(int node) {
		return taxon2node.containsValue(node);
	}

	public Set<String> taxa() {
		return taxon2node.keySet();
	}

	/** the node of a taxon, or -1 */
	public int nodeOf(String taxon) {
		return taxon2node.getOrDefault(taxon, -1);
	}

	public double length() {
		return Arrays.stream(weight).sum();
	}

	/** the number of independent cycles: edges - nodes + connected components */
	public int cycles() {
		var seen = new boolean[numberOfNodes];
		var components = 0;
		var stack = new ArrayDeque<Integer>();
		for (var s = 0; s < numberOfNodes; s++) {
			if (!seen[s]) {
				components++;
				seen[s] = true;
				stack.push(s);
				while (!stack.isEmpty()) {
					var v = stack.pop();
					for (var e : adjacentEdges[v]) {
						var w = source[e] == v ? target[e] : source[e];
						if (!seen[w]) {
							seen[w] = true;
							stack.push(w);
						}
					}
				}
			}
		}
		return source.length - numberOfNodes + components;
	}

	/** the distance of every node from the nearest of the given nodes (Dijkstra with several sources) */
	public float[] distancesFrom(int[] starts) {
		var distance = new double[numberOfNodes];
		Arrays.fill(distance, Double.POSITIVE_INFINITY);
		var queue = new PriorityQueue<double[]>(Comparator.comparingDouble(a -> a[0]));
		for (var start : starts) {
			distance[start] = 0;
			queue.add(new double[]{0, start});
		}
		relax(distance, queue);
		var result = new float[numberOfNodes];
		for (var v = 0; v < numberOfNodes; v++)
			result[v] = (float) distance[v];
		return result;
	}

	/** shortest-path distances from a node to all nodes (Dijkstra); unreachable nodes are at infinity */
	public double[] distancesFrom(int start) {
		var distance = new double[numberOfNodes];
		Arrays.fill(distance, Double.POSITIVE_INFINITY);
		distance[start] = 0;
		var queue = new PriorityQueue<double[]>(Comparator.comparingDouble(a -> a[0]));
		queue.add(new double[]{0, start});
		relax(distance, queue);
		return distance;
	}

	private void relax(double[] distance, PriorityQueue<double[]> queue) {
		while (!queue.isEmpty()) {
			var top = queue.poll();
			var v = (int) top[1];
			if (top[0] > distance[v])
				continue;
			for (var e : adjacentEdges[v]) {
				var w = source[e] == v ? target[e] : source[e];
				var d = top[0] + weight[e];
				if (d < distance[w]) {
					distance[w] = d;
					queue.add(new double[]{d, w});
				}
			}
		}
	}

	/** the distances between the given taxa, as displayed by the graph; NaN for a taxon that is not in the graph */
	public double[][] taxonDistances(List<String> taxa) {
		var n = taxa.size();
		var result = new double[n][n];
		for (var i = 0; i < n; i++) {
			var u = nodeOf(taxa.get(i));
			var from = u >= 0 ? distancesFrom(u) : null;
			for (var j = 0; j < n; j++) {
				var v = nodeOf(taxa.get(j));
				result[i][j] = (from == null || v < 0) ? Double.NaN : from[v];
			}
		}
		return result;
	}
}
