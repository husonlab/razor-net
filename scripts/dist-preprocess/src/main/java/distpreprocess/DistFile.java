package distpreprocess;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Labeled square .dist files used by RazorNet (relaxed Phylip-like format).
 * First non-comment line is n; each of the next n data lines is a taxon label
 * followed by n distances. Lines starting with '#' and blank lines are ignored.
 */
public final class DistFile {
	private DistFile() {
	}

	public record Matrix(List<String> labels, double[][] distances) {
		public int size() {
			return labels.size();
		}

		public double[][] copyDistances() {
			int n = distances.length;
			double[][] copy = new double[n][n];
			for (int i = 0; i < n; i++) {
				System.arraycopy(distances[i], 0, copy[i], 0, n);
			}
			return copy;
		}
	}

	public static Matrix read(Path path) throws IOException {
		List<String> lines = Files.readAllLines(path, StandardCharsets.UTF_8);
		int n = -1;
		int lineNo = 0;
		List<String> labels = new ArrayList<>();
		double[][] distances = null;
		int row = 0;

		for (String raw : lines) {
			lineNo++;
			String line = raw.trim();
			if (line.isEmpty() || line.startsWith("#")) {
				continue;
			}
			if (n < 0) {
				try {
					n = Integer.parseInt(line);
				} catch (NumberFormatException e) {
					throw new DistPreprocessException(
							"Line " + lineNo + ": first non-comment line must be the integer taxon count, found: " + line);
				}
				if (n <= 0) {
					throw new DistPreprocessException("Matrix dimension must be positive, found n=" + n);
				}
				distances = new double[n][n];
				continue;
			}
			if (row >= n) {
				throw new DistPreprocessException("Line " + lineNo + ": more than n=" + n + " data rows");
			}
			String[] tokens = line.split("\\s+");
			if (tokens.length != n + 1) {
				throw new DistPreprocessException(
						"Line " + lineNo + ": expected label plus " + n + " distances (square matrix with labels), found "
								+ tokens.length + " tokens. Lower/upper-triangle-only files are not supported.");
			}
			labels.add(tokens[0]);
			for (int j = 0; j < n; j++) {
				try {
					double value = Double.parseDouble(tokens[j + 1]);
					if (!Double.isFinite(value)) {
						throw new DistPreprocessException(
								"Line " + lineNo + ", column " + (j + 1) + ": non-finite value " + tokens[j + 1]);
					}
					distances[row][j] = value;
				} catch (NumberFormatException e) {
					throw new DistPreprocessException(
							"Line " + lineNo + ", column " + (j + 1) + ": not a number: " + tokens[j + 1]);
				}
			}
			row++;
		}

		if (n < 0) {
			throw new DistPreprocessException("File contains no taxon count");
		}
		if (row < n) {
			throw new DistPreprocessException("Expected " + n + " data rows but found " + row);
		}
		return new Matrix(List.copyOf(labels), distances);
	}

	/**
	 * Writes a labeled square matrix. {@code formatDistance} is used for every entry,
	 * including the diagonal (which the caller should already have set to 0).
	 */
	public static void write(Path path, List<String> headerComments, List<String> labels, double[][] distances,
							 DistanceFormatter formatter) throws IOException {
		int n = labels.size();
		if (distances.length != n) {
			throw new DistPreprocessException("Label count " + n + " does not match matrix size " + distances.length);
		}
		StringBuilder sb = new StringBuilder();
		for (String comment : headerComments) {
			if (!comment.startsWith("#")) {
				sb.append("# ");
			}
			sb.append(comment).append('\n');
		}
		sb.append(n).append('\n');
		for (int i = 0; i < n; i++) {
			if (distances[i].length != n) {
				throw new DistPreprocessException("Matrix is not square at row " + i);
			}
			sb.append(labels.get(i)).append('\t');
			for (int j = 0; j < n; j++) {
				if (j > 0) {
					sb.append(' ');
				}
				sb.append(formatter.format(i, j, distances[i][j]));
			}
			sb.append('\n');
		}
		Path parent = path.getParent();
		if (parent != null) {
			Files.createDirectories(parent);
		}
		Files.writeString(path, sb.toString(), StandardCharsets.UTF_8);
	}

	@FunctionalInterface
	public interface DistanceFormatter {
		String format(int row, int col, double value);
	}
}
