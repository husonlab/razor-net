/*
 * CharactersToDistances.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.progress.ProgressSilent;
import splitstree6.data.CharactersBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.ImportManager;
import splitstree6.io.utils.DataReaderBase;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;

/**
 * Writes the Hamming distances of an alignment as a Nexus distances file, prepared as the comparison of the
 * haplotype network methods prepares its input ({@link CompareHaplotypeMethods}): every column in which some
 * sequence has a gap, a missing or an ambiguous character is dropped, and identical sequences are merged into one
 * taxon, named after the first of them. For {@link CompareDistanceMethods}, which reads distances.
 *
 * <pre>
 * CharactersToDistances input output.nex
 * </pre>
 *
 * @author Daniel Huson, 2026
 */
public class CharactersToDistances {
	public static void main(String[] args) throws Exception {
		if (args.length != 2) {
			System.err.println("Usage: CharactersToDistances input output.nex");
			System.exit(1);
		}
		var input = args[0];
		var output = args[1];
		var taxa = new TaxaBlock();
		var characters = new CharactersBlock();
		var reader = ImportManager.getInstance().getReader(input);
		if (reader == null || !CharactersBlock.class.isAssignableFrom(reader.getToClass()))
			throw new Exception("Not an alignment: " + input);
		((DataReaderBase<CharactersBlock>) reader).read(new ProgressSilent(), input, taxa, characters);
		var n = characters.getNtax();
		var m = characters.getNchar();
		var incomplete = "-?NRYSWKMBDHV" + characters.getGapCharacter() + characters.getMissingCharacter();
		var kept = new ArrayList<Integer>();
		for (var c = 1; c <= m; c++) {
			var complete = true;
			for (var t = 1; t <= n && complete; t++)
				complete = incomplete.indexOf(Character.toUpperCase(characters.get(t, c))) < 0;
			if (complete)
				kept.add(c);
		}
		var labels = new LinkedHashMap<String, String>(); // sequence on the kept columns -> label of its first taxon
		for (var t = 1; t <= n; t++) {
			var buf = new StringBuilder();
			for (var c : kept)
				buf.append(Character.toUpperCase(characters.get(t, c)));
			labels.putIfAbsent(buf.toString(), taxa.getLabel(t));
		}
		var sequences = new ArrayList<>(labels.keySet());
		var k = sequences.size();
		var d = new int[k][k];
		for (var i = 0; i < k; i++)
			for (var j = i + 1; j < k; j++) {
				var differences = 0;
				for (var p = 0; p < sequences.get(i).length(); p++)
					if (sequences.get(i).charAt(p) != sequences.get(j).charAt(p))
						differences++;
				d[i][j] = d[j][i] = differences;
			}
		try (var w = new PrintWriter(output)) {
			w.println("#nexus");
			w.printf("[! Hamming distances of %s: %d of %d columns kept (no gap, missing or ambiguous character), %d sequences, %d distinct]%n",
					input, kept.size(), m, n, k);
			w.printf("BEGIN TAXA;%nDIMENSIONS ntax=%d;%nTAXLABELS%n", k);
			for (var label : labels.values())
				w.printf("'%s'%n", label.replace("'", ""));
			w.printf(";%nEND;%nBEGIN DISTANCES;%nDIMENSIONS ntax=%d;%nFORMAT labels=left diagonal triangle=both;%nMATRIX%n", k);
			var i = 0;
			for (var label : labels.values()) {
				w.printf("'%s'", label.replace("'", ""));
				for (var j = 0; j < k; j++)
					w.print(" " + d[i][j]);
				w.println();
				i++;
			}
			w.printf(";%nEND;%n");
		}
		System.out.printf("%s: %d of %d columns kept, %d sequences, %d distinct, written to %s%n", input, kept.size(), m, n, k, output);
	}
}
