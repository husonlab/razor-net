package distpreprocess;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class DistFileTest {
	@TempDir
	Path tmp;

	@Test
	void roundTripPreservesLabelsAndValues() throws Exception {
		Path input = tmp.resolve("in.dist");
		Files.writeString(input, """
				# comment
				3
				alpha	0 1.5 2
				beta	1.5 0 3
				gamma	2 3 0
				""");
		String original = Files.readString(input);
		DistFile.Matrix m = DistFile.read(input);
		assertEquals(List.of("alpha", "beta", "gamma"), m.labels());
		assertEquals(1.5, m.distances()[0][1]);

		Path output = tmp.resolve("out.dist");
		DistFile.write(output, List.of("# dist-preprocess"), m.labels(), m.distances(),
				(r, c, v) -> r == c ? "0" : Double.toString(v));
		DistFile.Matrix again = DistFile.read(output);
		assertEquals(m.labels(), again.labels());
		assertEquals(1.5, again.distances()[0][1]);
		assertEquals(original, Files.readString(input));
	}

	@Test
	void rejectsTriangleOnlyMatrix() throws Exception {
		Path input = tmp.resolve("lower.dist");
		Files.writeString(input, """
				3
				a	0
				b	1 0
				c	2 3 0
				""");
		assertThrows(DistPreprocessException.class, () -> DistFile.read(input));
	}

	@Test
	void readsLeadingWhitespaceOnN() throws Exception {
		Path input = tmp.resolve("n.dist");
		Files.writeString(input, """
				  2
				x 0 4
				y 4 0
				""");
		DistFile.Matrix m = DistFile.read(input);
		assertEquals(2, m.size());
		assertEquals("x", m.labels().get(0));
		assertEquals(4.0, m.distances()[0][1]);
	}

	@Test
	void writtenFileStartsWithProvenanceAndN() throws Exception {
		Path output = tmp.resolve("w.dist");
		DistFile.write(output,
				List.of("# dist-preprocess", "# source=foo.dist"),
				List.of("a", "b"),
				new double[][]{{0, 1}, {1, 0}},
				(r, c, v) -> r == c ? "0" : "1");
		String text = Files.readString(output);
		assertTrue(text.startsWith("# dist-preprocess\n# source=foo.dist\n2\n"));
		assertTrue(text.contains("a\t0 1"));
	}
}
