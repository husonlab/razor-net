/*
 * NewickTree.java Copyright (C) 2026 Daniel H. Huson
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
import java.util.ArrayList;

/**
 * A rooted tree read from a Newick string: for each node its parent, label and the length of the edge to its
 * parent. Node 0 is the root.
 *
 * @author Daniel Huson, 2026
 */
public final class NewickTree {
	private final ArrayList<Integer> parent = new ArrayList<>();
	private final ArrayList<String> label = new ArrayList<>();
	private final ArrayList<Double> length = new ArrayList<>();
	private final ArrayList<Integer> children = new ArrayList<>(); // number of children
	private String text;
	private int position;

	private NewickTree() {
	}

	public static NewickTree parse(String newick) throws IOException {
		var tree = new NewickTree();
		tree.text = newick.trim();
		tree.position = 0;
		tree.parseNode(-1);
		if (tree.position >= tree.text.length() || tree.text.charAt(tree.position) != ';')
			throw new IOException("Newick: ';' expected at position " + tree.position);
		return tree;
	}

	private void parseNode(int parentNode) throws IOException {
		var v = parent.size();
		parent.add(parentNode);
		label.add(null);
		length.add(Double.NaN);
		children.add(0);
		if (parentNode >= 0)
			children.set(parentNode, children.get(parentNode) + 1);
		if (peek() == '(') {
			do {
				position++; // '(' or ','
				parseNode(v);
			} while (peek() == ',');
			if (peek() != ')')
				throw new IOException("Newick: ')' expected at position " + position);
			position++;
		}
		var start = position;
		while (position < text.length() && ",():;".indexOf(text.charAt(position)) < 0)
			position++;
		if (position > start)
			label.set(v, text.substring(start, position).trim());
		if (peek() == ':') {
			start = ++position;
			while (position < text.length() && ",();".indexOf(text.charAt(position)) < 0)
				position++;
			length.set(v, Double.parseDouble(text.substring(start, position).trim()));
		}
	}

	private char peek() {
		return position < text.length() ? text.charAt(position) : 0;
	}

	public int size() {
		return parent.size();
	}

	/** the parent of a node, -1 for the root */
	public int parent(int v) {
		return parent.get(v);
	}

	public String label(int v) {
		return label.get(v);
	}

	/** the length of the edge from a node to its parent, NaN if none was given */
	public double length(int v) {
		return length.get(v);
	}

	public boolean isLeaf(int v) {
		return children.get(v) == 0;
	}
}
