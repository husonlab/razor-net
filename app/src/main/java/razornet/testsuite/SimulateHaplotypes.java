/*
 * SimulateHaplotypes.java Copyright (C) 2026 Daniel H. Huson
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
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Random;

/**
 * Generates samples of sequences from a population with a known genealogy. The genealogy of the sample is drawn
 * from the standard neutral coalescent (Kingman 1982), and sequences are evolved down it with AliSim (Ly-Trong et
 * al. 2022, IQ-TREE), or under the infinite-sites model, in which no site mutates twice, so that the sequence of
 * every ancestor is known. The mutation rate is given as theta, the expected number of sites at which two sampled
 * sequences differ if no site mutates twice.
 * <p>
 * For each replicate the tool writes a directory with {@code samples.fasta} (the sampled sequences),
 * {@code genealogy.nwk} (all nodes named, edge lengths in expected substitutions per site), {@code nodes.fasta}
 * (the sequences of all nodes of the genealogy) and {@code info.txt}.
 *
 * <pre>
 * SimulateHaplotypes --output directory [--samples N] [--sites L] [--theta T] [--replicates R] [--seed S]
 *                    [--model M | --infiniteSites] [--iqtree program]
 * </pre>
 *
 * @author Daniel Huson, 2026
 */
public class SimulateHaplotypes {
	/** a coalescent genealogy: going back in time, two random lineages join at rate k(k-1)/2 for k lineages */
	static SimulateTreeMixtures.Tree coalescentTree(int leaves, Random random) {
		var tree = new SimulateTreeMixtures.Tree(leaves);
		var active = new ArrayList<Integer>();
		for (var v = 0; v < leaves; v++)
			active.add(v);
		var time = 0.0;
		var next = leaves;
		while (active.size() > 1) {
			var k = active.size();
			time += -Math.log(1.0 - random.nextDouble()) / (k * (k - 1) / 2.0);
			var a = active.remove(random.nextInt(active.size()));
			var b = active.remove(random.nextInt(active.size()));
			tree.children[next] = new int[]{a, b};
			tree.parent[a] = tree.parent[b] = next;
			tree.height[next] = time;
			active.add(next++);
		}
		tree.root = active.get(0);
		return tree;
	}

	public static void main(String[] args) throws Exception {
		var output = "";
		var samples = 50;
		var sites = 1000;
		var theta = 10.0;
		var replicates = 10;
		var seed = 1L;
		var model = "HKY{2.0}+F{0.3/0.2/0.2/0.3}+G4{0.5}";
		var infiniteSites = false;
		var iqtree = "iqtree2";
		for (var i = 0; i < args.length; i++) {
			switch (args[i]) {
				case "--output", "-o" -> output = args[++i];
				case "--infiniteSites" -> infiniteSites = true;
				case "--samples", "-n" -> samples = Integer.parseInt(args[++i]);
				case "--sites", "-l" -> sites = Integer.parseInt(args[++i]);
				case "--theta" -> theta = Double.parseDouble(args[++i]);
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
		if (output.isEmpty() || samples < 2) {
			System.err.println("Usage: SimulateHaplotypes --output directory [--samples N] [--sites L] [--theta T] [--replicates R] [--seed S] [--model M | --infiniteSites] [--iqtree program]");
			System.exit(2);
		}
		var names = new String[samples];
		for (var t = 0; t < samples; t++)
			names[t] = "s%03d".formatted(t + 1);

		var random = new Random(seed);
		for (var replicate = 1; replicate <= replicates; replicate++) {
			var directory = new File(output, "rep%03d".formatted(replicate));
			Files.createDirectories(directory.toPath());
			var tree = coalescentTree(samples, random);
			// two lineages coalesce after time 1 on average and then differ at theta sites
			var perSite = theta / (2.0 * sites);
			var newick = tree.toNewick(names, v -> "%.10f".formatted((tree.height[tree.parent[v]] - tree.height[v]) * perSite));
			var sequences = infiniteSites
					? SimulateTreeMixtures.evolveInfiniteSites(tree, names, v -> (tree.height[tree.parent[v]] - tree.height[v]) * perSite, sites, random)
					: SimulateTreeMixtures.evolve(iqtree, newick, model, sites, random.nextInt(Integer.MAX_VALUE), directory);
			var sampled = new LinkedHashMap<String, String>();
			for (var name : names)
				sampled.put(name, sequences.get(name));
			Files.writeString(new File(directory, "genealogy.nwk").toPath(), newick + "\n");
			SimulateTreeMixtures.writeFasta(new File(directory, "nodes.fasta"), sequences);
			SimulateTreeMixtures.writeFasta(new File(directory, "samples.fasta"), sampled);
			Files.writeString(new File(directory, "info.txt").toPath(),
					"samples=%d%nsites=%d%ntheta=%s%nmodel=%s%nseed=%d%nreplicate=%d%n".formatted(samples, sites, theta,
							infiniteSites ? "infinite sites" : model, seed, replicate));
		}
		System.err.printf("Wrote %d replicates to %s%n", replicates, output);
	}
}
