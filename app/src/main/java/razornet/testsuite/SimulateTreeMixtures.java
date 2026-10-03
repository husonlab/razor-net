/*
 * SimulateTreeMixtures.java Copyright (C) 2026 Daniel H. Huson
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

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Random;

/**
 * Generates data sets with a known conflicting signal: sequences that evolved down several trees which each
 * differ from one species tree by a single rooted subtree-prune-and-regraft move, concatenated.
 * <p>
 * The species tree is a Yule tree, scaled to a given height in expected substitutions per site. A gene tree is
 * obtained from it by pruning a random subtree and regrafting it onto a lineage that exists at the same time, as in
 * a horizontal transfer. Sequences are evolved down each gene tree with AliSim (Ly-Trong et al. 2022, IQ-TREE),
 * or under the infinite-sites model, in which no site mutates twice, and the alignments of the gene trees are
 * concatenated; the gene trees share the alignment length equally.
 * <p>
 * For each replicate the tool writes a directory with {@code species.nwk}, {@code genes.nwk} (one tree per line),
 * {@code genes.mutations.nwk} (the gene trees with the number of mutations on each edge as its length, counted as
 * the sites at which the sequences at its two ends differ), {@code alignment.fasta}, {@code hamming.dist} (the
 * number of differing sites for each pair of taxa) and {@code info.txt}.
 *
 * <pre>
 * SimulateTreeMixtures --output directory [--taxa N] [--trees K] [--length L] [--height H] [--replicates R]
 *                      [--seed S] [--model M | --infiniteSites] [--iqtree program]
 * </pre>
 *
 * @author Daniel Huson, 2026
 */
public class SimulateTreeMixtures {
	/** a rooted binary tree with node heights; leaves are 0..n-1 */
	static final class Tree {
		final int leaves;
		final int[] parent;
		final int[][] children;
		final double[] height;
		int root;

		Tree(int leaves) {
			this.leaves = leaves;
			parent = new int[2 * leaves - 1];
			children = new int[2 * leaves - 1][];
			height = new double[2 * leaves - 1];
			java.util.Arrays.fill(parent, -1);
		}

		Tree copy() {
			var that = new Tree(leaves);
			System.arraycopy(parent, 0, that.parent, 0, parent.length);
			System.arraycopy(height, 0, that.height, 0, height.length);
			for (var v = 0; v < children.length; v++)
				that.children[v] = children[v] == null ? null : children[v].clone();
			that.root = root;
			return that;
		}

		boolean isBelow(int v, int ancestor) {
			for (; v >= 0; v = parent[v])
				if (v == ancestor)
					return true;
			return false;
		}

		int sibling(int v) {
			var p = parent[v];
			return children[p][0] == v ? children[p][1] : children[p][0];
		}

		void replaceChild(int p, int oldChild, int newChild) {
			if (children[p][0] == oldChild)
				children[p][0] = newChild;
			else
				children[p][1] = newChild;
			parent[newChild] = p;
		}

		/** the name of a node: that of its leaf, or i followed by its number */
		String name(int v, String[] leafNames) {
			return v < leaves ? leafNames[v] : "i" + v;
		}

		/** Newick, with the difference of the heights as the length of an edge */
		String toNewick(String[] leafNames) {
			return toNewick(leafNames, v -> "%.8f".formatted(height[parent[v]] - height[v]));
		}

		/**
		 * Newick, with all nodes named
		 *
		 * @param edgeLength the length of the edge from a node to its parent, as text
		 */
		String toNewick(String[] leafNames, java.util.function.IntFunction<String> edgeLength) {
			var buf = new StringBuilder();
			write(root, leafNames, edgeLength, buf);
			return buf.append(";").toString();
		}

		private void write(int v, String[] leafNames, java.util.function.IntFunction<String> edgeLength, StringBuilder buf) {
			if (children[v] != null) {
				buf.append("(");
				write(children[v][0], leafNames, edgeLength, buf);
				buf.append(",");
				write(children[v][1], leafNames, edgeLength, buf);
				buf.append(")");
			}
			buf.append(name(v, leafNames));
			if (parent[v] >= 0)
				buf.append(":").append(edgeLength.apply(v));
		}
	}

	/** a Yule tree: going back in time, two random lineages join after a waiting time with rate k for k lineages */
	static Tree yuleTree(int leaves, double rootHeight, Random random) {
		var tree = new Tree(leaves);
		var active = new ArrayList<Integer>();
		for (var v = 0; v < leaves; v++)
			active.add(v);
		var time = 0.0;
		var next = leaves;
		while (active.size() > 1) {
			time += -Math.log(1.0 - random.nextDouble()) / active.size();
			var a = active.remove(random.nextInt(active.size()));
			var b = active.remove(random.nextInt(active.size()));
			tree.children[next] = new int[]{a, b};
			tree.parent[a] = tree.parent[b] = next;
			tree.height[next] = time;
			active.add(next++);
		}
		tree.root = active.get(0);
		var factor = rootHeight / tree.height[tree.root];
		for (var v = 0; v < tree.height.length; v++)
			tree.height[v] *= factor;
		return tree;
	}

	/**
	 * a copy of the tree in which a random subtree has been pruned and regrafted onto a lineage that exists at the
	 * time of the cut, and that is not the lineage of its sibling; null if no such move was found
	 */
	static Tree subtreePruneAndRegraft(Tree species, Random random) {
		for (var attempt = 0; attempt < 1000; attempt++) {
			var tree = species.copy();
			var v = random.nextInt(tree.parent.length);
			if (v == tree.root)
				continue;
			var p = tree.parent[v];
			var s = tree.sibling(v);
			var time = tree.height[v] + random.nextDouble() * (tree.height[p] - tree.height[v]);
			var targets = new ArrayList<Integer>();
			for (var c = 0; c < tree.parent.length; c++)
				if (c != tree.root && c != s && c != p && !tree.isBelow(c, v)
					&& tree.height[c] <= time && time < tree.height[tree.parent[c]])
					targets.add(c);
			if (targets.isEmpty())
				continue;
			var c = targets.get(random.nextInt(targets.size()));
			// remove p: its other child takes its place
			if (p == tree.root) {
				tree.root = s;
				tree.parent[s] = -1;
			} else
				tree.replaceChild(tree.parent[p], p, s);
			// p is reused as the new parent of v and c, at the time of the cut
			var above = tree.parent[c];
			tree.replaceChild(above, c, p);
			tree.children[p] = new int[]{c, v};
			tree.parent[c] = tree.parent[v] = p;
			tree.height[p] = time;
			return tree;
		}
		return null;
	}

	/** runs AliSim on a tree and returns the sequences of all its nodes by name */
	static LinkedHashMap<String, String> evolve(String iqtree, String newick, String model, int length, long seed,
												 File directory) throws IOException, InterruptedException {
		var treeFile = new File(directory, "alisim-tree.nwk");
		Files.writeString(treeFile.toPath(), newick + "\n");
		var command = new ArrayList<>(List.of(iqtree, "--alisim", "alisim-out", "-t", treeFile.getName(), "-m", model,
				"--length", String.valueOf(length), "--seed", String.valueOf(seed), "-af", "fasta", "-redo", "-nt", "1",
				"--write-all"));
		var process = new ProcessBuilder(command).directory(directory).redirectErrorStream(true)
				.redirectOutput(new File(directory, "alisim-run.txt")).start();
		if (process.waitFor() != 0)
			throw new IOException("AliSim failed, see " + new File(directory, "alisim-run.txt"));
		var sequences = readFasta(new File(directory, "alisim-out.fa"));
		for (var file : directory.listFiles((dir, name) -> name.startsWith("alisim-")))
			file.delete();
		return sequences;
	}

	/**
	 * evolves sequences down a tree under the infinite-sites model: the number of mutations on an edge is Poisson
	 * with the length of the edge times the number of sites as its mean, and every mutation is at a site of its
	 * own, at which the nodes below the edge have C and all others A. Returns the sequences of all nodes by name.
	 */
	static LinkedHashMap<String, String> evolveInfiniteSites(Tree tree, String[] leafNames, java.util.function.IntToDoubleFunction edgeLength,
															  int sites, Random random) {
		var nodes = tree.parent.length;
		var buffers = new StringBuilder[nodes];
		for (var v = 0; v < nodes; v++)
			buffers[v] = new StringBuilder();
		for (var v = 0; v < nodes; v++) {
			if (v == tree.root)
				continue;
			var mean = edgeLength.applyAsDouble(v) * sites;
			// the number of events of a Poisson process of rate 1 before time mean
			var mutations = 0;
			for (var time = -Math.log(1.0 - random.nextDouble()); time < mean; time += -Math.log(1.0 - random.nextDouble()))
				mutations++;
			for (var m = 0; m < mutations; m++)
				for (var w = 0; w < nodes; w++)
					buffers[w].append(tree.isBelow(w, v) ? 'C' : 'A');
		}
		var sequences = new LinkedHashMap<String, String>();
		for (var v = 0; v < nodes; v++)
			sequences.put(tree.name(v, leafNames), buffers[v].toString());
		return sequences;
	}

	static LinkedHashMap<String, String> readFasta(File file) throws IOException {
		var sequences = new LinkedHashMap<String, String>();
		String name = null;
		var buf = new StringBuilder();
		for (var line : Files.readAllLines(file.toPath())) {
			if (line.startsWith(">")) {
				if (name != null)
					sequences.put(name, buf.toString());
				name = line.substring(1).trim();
				buf.setLength(0);
			} else
				buf.append(line.trim());
		}
		if (name != null)
			sequences.put(name, buf.toString());
		return sequences;
	}

	static int differences(String x, String y) {
		var differences = 0;
		for (var j = 0; j < x.length(); j++)
			if (x.charAt(j) != y.charAt(j))
				differences++;
		return differences;
	}

	static void writeFasta(File file, LinkedHashMap<String, String> sequences) throws IOException {
		try (var w = new PrintWriter(file)) {
			for (var entry : sequences.entrySet())
				w.printf(">%s%n%s%n", entry.getKey(), entry.getValue());
		}
	}

	public static void main(String[] args) throws Exception {
		var output = "";
		var taxa = 20;
		var trees = 3;
		var length = 6000;
		var height = 0.1;
		var replicates = 10;
		var seed = 1L;
		var model = "HKY{2.0}+F{0.3/0.2/0.2/0.3}+G4{0.5}";
		var infiniteSites = false;
		var iqtree = "iqtree2";
		for (var i = 0; i < args.length; i++) {
			switch (args[i]) {
				case "--output", "-o" -> output = args[++i];
				case "--infiniteSites" -> infiniteSites = true;
				case "--taxa", "-n" -> taxa = Integer.parseInt(args[++i]);
				case "--trees", "-k" -> trees = Integer.parseInt(args[++i]);
				case "--length", "-l" -> length = Integer.parseInt(args[++i]);
				case "--height" -> height = Double.parseDouble(args[++i]);
				case "--replicates", "-r" -> replicates = Integer.parseInt(args[++i]);
				case "--seed" -> seed = Long.parseLong(args[++i]);
				case "--model", "-m" -> model = args[++i];
				case "--iqtree" -> iqtree = args[++i];
				default -> {
					System.err.println("Unknown argument: " + args[i]);
					System.exit(2);
				}
			}
		}
		if (output.isEmpty() || taxa < 4 || trees < 1) {
			System.err.println("Usage: SimulateTreeMixtures --output directory [--taxa N] [--trees K] [--length L] [--height H] [--replicates R] [--seed S] [--model M | --infiniteSites] [--iqtree program]");
			System.exit(2);
		}
		var names = new String[taxa];
		for (var t = 0; t < taxa; t++)
			names[t] = "t%03d".formatted(t + 1);

		var random = new Random(seed);
		for (var replicate = 1; replicate <= replicates; replicate++) {
			var directory = new File(output, "rep%03d".formatted(replicate));
			Files.createDirectories(directory.toPath());
			var species = yuleTree(taxa, height, random);
			var genes = new ArrayList<Tree>();
			while (genes.size() < trees) {
				var gene = subtreePruneAndRegraft(species, random);
				if (gene == null)
					throw new IOException("No subtree-prune-and-regraft move found");
				genes.add(gene);
			}

			var alignment = new LinkedHashMap<String, String>();
			for (var name : names)
				alignment.put(name, "");
			var lines = new ArrayList<String>();
			var mutationLines = new ArrayList<String>();
			for (var g = 0; g < trees; g++) {
				var gene = genes.get(g);
				var newick = gene.toNewick(names);
				lines.add(newick);
				var block = length / trees + (g < length % trees ? 1 : 0);
				var sequences = infiniteSites
						? evolveInfiniteSites(gene, names, v -> gene.height[gene.parent[v]] - gene.height[v], block, random)
						: evolve(iqtree, newick, model, block, random.nextInt(Integer.MAX_VALUE), directory);
				for (var name : names)
					alignment.put(name, alignment.get(name) + sequences.get(name));
				// the same tree with the number of sites that differ between the two ends of an edge as its length
				mutationLines.add(gene.toNewick(names, v -> String.valueOf(
						differences(sequences.get(gene.name(v, names)), sequences.get(gene.name(gene.parent[v], names))))));
			}
			Files.writeString(new File(directory, "species.nwk").toPath(), species.toNewick(names) + "\n");
			Files.write(new File(directory, "genes.nwk").toPath(), lines);
			Files.write(new File(directory, "genes.mutations.nwk").toPath(), mutationLines);
			writeFasta(new File(directory, "alignment.fasta"), alignment);
			try (var w = new PrintWriter(new File(directory, "hamming.dist"))) {
				w.println(taxa);
				for (var a : names) {
					w.print(a);
					for (var b : names)
						w.print(" " + differences(alignment.get(a), alignment.get(b)));
					w.println();
				}
			}
			Files.writeString(new File(directory, "info.txt").toPath(),
					"taxa=%d%ntrees=%d%nlength=%d%nheight=%s%nmodel=%s%nseed=%d%nreplicate=%d%n".formatted(taxa, trees, length, height,
							infiniteSites ? "infinite sites" : model, seed, replicate));
		}
		System.err.printf("Wrote %d replicates to %s%n", replicates, output);
	}
}
