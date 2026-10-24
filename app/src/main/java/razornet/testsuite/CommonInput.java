/*
 *  CommonInput.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.progress.ProgressSilent;
import splitstree6.data.CharactersBlock;
import splitstree6.data.TaxaBlock;
import splitstree6.io.readers.ImportManager;

import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

/**
 * The common input of the haplotype network comparisons: the alignment without every column that contains a
 * gap, a missing or an ambiguous character (complete deletion), with identical sequences merged into one
 * haplotype. Used by {@link CompareHaplotypeMethods} and {@link PopartCheck}.
 *
 * @param sequences  the number of sequences of the alignment
 * @param kept       the columns kept, 1-based
 * @param masked     every sequence on the kept columns, by its label in the alignment
 * @param haplotypes the distinct sequences, by the label of their first member, made safe for file names
 * @param members    the labels of the sequences of each haplotype
 * @author Daniel Huson, 10.2026
 */
public record CommonInput(int sequences, List<Integer> kept, LinkedHashMap<String, String> masked,
						  LinkedHashMap<String, String> haplotypes, LinkedHashMap<String, List<String>> members) {

	public static CommonInput prepare(String input) throws IOException {
		var allTaxa = new TaxaBlock();
		var allCharacters = new CharactersBlock();
		read(input, allTaxa, allCharacters);
		// a gap, a missing character, N or an IUPAC ambiguity code; any other symbol is a state (so are 0 and 1)
		var incomplete = "-?NRYSWKMBDHV" + Character.toUpperCase(allCharacters.getGapCharacter())
						 + Character.toUpperCase(allCharacters.getMissingCharacter());
		var kept = new ArrayList<Integer>();
		for (var j = 1; j <= allCharacters.getNchar(); j++) {
			var complete = true;
			for (var t = 1; t <= allTaxa.getNtax() && complete; t++)
				complete = incomplete.indexOf(Character.toUpperCase(allCharacters.get(t, j))) < 0;
			if (complete)
				kept.add(j);
		}
		var masked = new LinkedHashMap<String, String>();
		var haplotypes = new LinkedHashMap<String, String>(); // label -> sequence
		var sequence2label = new LinkedHashMap<String, String>();
		var members = new LinkedHashMap<String, List<String>>();
		for (var t = 1; t <= allTaxa.getNtax(); t++) {
			var buf = new StringBuilder();
			for (var j : kept)
				buf.append(Character.toUpperCase(allCharacters.get(t, j)));
			masked.put(allTaxa.getLabel(t), buf.toString());
			var label = sequence2label.get(buf.toString());
			if (label == null) {
				label = allTaxa.getLabel(t).replaceAll("[^A-Za-z0-9_.-]", "_");
				while (haplotypes.containsKey(label))
					label += "_";
				sequence2label.put(buf.toString(), label);
				haplotypes.put(label, buf.toString());
				members.put(label, new ArrayList<>());
			}
			members.get(label).add(allTaxa.getLabel(t));
		}
		return new CommonInput(allTaxa.getNtax(), kept, masked, haplotypes, members);
	}

	/**
	 * reads an alignment in any format that SplitsTree imports
	 */
	static void read(String file, TaxaBlock taxa, CharactersBlock characters) throws IOException {
		var reader = ImportManager.getInstance().getReader(file);
		if (reader == null || !CharactersBlock.class.isAssignableFrom(reader.getToClass()))
			throw new IOException("Not an alignment: " + file);
		((splitstree6.io.utils.DataReaderBase<CharactersBlock>) reader).read(new ProgressSilent(), file, taxa, characters);
	}
}
