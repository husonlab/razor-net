/*
 *  PopartCheck.java Copyright (C) 2026 Daniel H. Huson
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
 */

package razornet.testsuite;

import jloda.util.FileUtils;
import jloda.util.progress.ProgressSilent;
import razornet.algorithms.tcs.TCSConstruction;
import razornet.algorithms.utils.NetworkSmoothing;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

/**
 * Writes the common input of the haplotype comparison ({@link CommonInput}) in the NEXUS dialect that PopART
 * reads, and reports the size of the TCS network that our construction builds on it, for a check against the
 * TCS network that PopART draws on the same file: the same numbers of inferred intermediates, steps (hatch
 * marks) and cycles are expected.
 * <p>
 * PopartCheck --out directory alignment...
 * <p>
 * Writes directory/name-popart.nex for each alignment and prints one table row per alignment. The numbers are
 * those of the construction before sequences are assigned to the intermediate nodes, which is what PopART
 * draws; the tables of the paper report the labeled network, whose number of mutations can differ.
 *
 * @author Daniel Huson, 10.2026
 */
public class PopartCheck {
	public static void main(String[] args) throws Exception {
		var out = "";
		var inputs = new ArrayList<String>();
		for (var i = 0; i < args.length; i++) {
			if (args[i].equals("--out"))
				out = args[++i];
			else
				inputs.add(args[i]);
		}
		if (out.isEmpty() || inputs.isEmpty()) {
			System.err.println("Usage: PopartCheck --out directory alignment...");
			System.exit(1);
		}
		if (!new File(out).isDirectory() && !new File(out).mkdirs())
			throw new IOException("Cannot create directory: " + out);

		System.out.println("| data set | sequences | haplotypes | sites | intermediates | steps | cycles | file |");
		System.out.println("|---|---:|---:|---:|---:|---:|---:|---|");
		for (var input : inputs) {
			var base = FileUtils.replaceFileSuffix(new File(input).getName(), "").replace(' ', '_');
			var common = CommonInput.prepare(input);
			var file = new File(out, base + "-popart.nex");
			writePopart(common, file);

			var labels = new ArrayList<>(common.haplotypes().keySet());
			var n = labels.size();
			var distances = new int[n][n];
			for (var i = 0; i < n; i++)
				for (var j = 0; j < n; j++)
					distances[i][j] = (int) Math.round(TaxonGraph.differences(common.haplotypes().get(labels.get(i)),
							common.haplotypes().get(labels.get(j))));
			var result = TCSConstruction.apply(distances, new ProgressSilent());
			var graph = result.graph();
			var steps = 0;
			for (var e : graph.edges())
				steps += graph.getWeight(e);
			var cycles = graph.countNumberOfEdges() - graph.getNumberOfNodes() + 1;
			// PopART draws the intermediates of degree three or more; the others are hatch marks on an edge
			var smoothed = NetworkSmoothing.apply(graph.copy(), id -> id < n, true, false, new ProgressSilent());
			var intermediates = result.intermediates() - smoothed.nodesRemoved();
			System.out.printf("| %s | %d | %d | %d | %d | %d | %d | %s |%n", base, common.sequences(), n,
					common.kept().size(), intermediates, steps, cycles, file.getName());
		}
	}

	/**
	 * every sequence on the kept columns, in the plainest NEXUS: labels of letters, digits and underscores,
	 * no ntax in the CHARACTERS block, no quotes
	 */
	private static void writePopart(CommonInput common, File file) throws IOException {
		var names = new ArrayList<String>();
		var used = new HashSet<String>();
		for (var label : common.masked().keySet()) {
			var name = label.replaceAll("[^A-Za-z0-9_]", "_");
			while (!used.add(name))
				name += "_";
			names.add(name);
		}
		try (var w = new PrintWriter(file)) {
			w.println("#NEXUS");
			w.println();
			w.println("BEGIN TAXA;");
			w.printf("DIMENSIONS NTAX=%d;%n", names.size());
			w.println("TAXLABELS");
			for (var name : names)
				w.println(name);
			w.println(";");
			w.println("END;");
			w.println();
			w.println("BEGIN CHARACTERS;");
			w.printf("DIMENSIONS NCHAR=%d;%n", common.kept().size());
			w.println("FORMAT DATATYPE=DNA MISSING=? GAP=-;");
			w.println("MATRIX");
			var i = 0;
			for (var sequence : common.masked().values())
				w.printf("%s %s%n", names.get(i++), sequence);
			w.println(";");
			w.println("END;");
		}
	}
}
