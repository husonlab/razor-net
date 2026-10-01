package distpreprocess;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class PreprocessDistancesTest {
	@TempDir
	Path tmp;

	@Test
	void parseRequiresMdsOrRound() {
		assertThrows(DistPreprocessException.class, () -> PreprocessDistances.parse(new String[]{
				"--input", "a.dist", "--output", "b.dist"
		}));
	}

	@Test
	void parseRejectsDimensionWithoutMds() {
		assertThrows(DistPreprocessException.class, () -> PreprocessDistances.parse(new String[]{
				"--input", "a.dist", "--output", "b.dist", "--dimension", "3", "--round", "2"
		}));
	}

	@Test
	void parseDefaultsDimensionToTwo() {
		PreprocessDistances.Options opt = PreprocessDistances.parse(new String[]{
				"-i", "a.dist", "-o", "b.dist", "--mds"
		});
		assertTrue(opt.mds());
		assertEquals(2, opt.dimension());
	}

	@Test
	void roundingOnlyDoesNotChangeInputFile() throws Exception {
		Path input = tmp.resolve("in.dist");
		Path output = tmp.resolve("out.dist");
		String original = """
				3
				a	0 3.333333 1
				b	3.333333 0 2
				c	1 2 0
				""";
		Files.writeString(input, original);

		int code = PreprocessDistances.run(new String[]{
				"--input", input.toString(),
				"--output", output.toString(),
				"--round", "2"
		}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(new ByteArrayOutputStream()));
		assertEquals(0, code);
		assertEquals(original, Files.readString(input));

		DistFile.Matrix out = DistFile.read(output);
		assertEquals(3.3, out.distances()[0][1], 1e-12);
		assertEquals("a", out.labels().get(0));
		String text = Files.readString(output);
		assertTrue(text.contains("3.3"));
		assertTrue(text.contains("# mds=false"));
		assertTrue(text.contains("# round=true significantDigits=2"));
	}

	@Test
	void mdsThenRoundAppliesInThatOrder() throws Exception {
		Path input = tmp.resolve("pts.dist");
		Path output = tmp.resolve("out.dist");
		double s = Math.sqrt(2.0);
		Files.writeString(input, """
				3
				p	0 1 %s
				q	1 0 1
				r	%s 1 0
				""".formatted(s, s));

		int code = PreprocessDistances.run(new String[]{
				"--input", input.toString(),
				"--output", output.toString(),
				"--mds", "--dimension", "2",
				"--round", "3"
		}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(new ByteArrayOutputStream()));
		assertEquals(0, code);

		DistFile.Matrix out = DistFile.read(output);
		// After exact 2D MDS the (p,r) distance is sqrt(2) ≈ 1.41421356, 3 sig figs -> 1.41
		assertEquals(1.41, out.distances()[0][2], 1e-9);
		assertEquals(0.0, out.distances()[0][0]);
		assertEquals(out.distances()[0][2], out.distances()[2][0]);
	}

	@Test
	void roundBeforeMdsDiffersFromDefaultAfterOrder() throws Exception {
		Path input = tmp.resolve("rect.dist");
		Path afterOut = tmp.resolve("after.dist");
		Path beforeOut = tmp.resolve("before.dist");
		// Four corners of a rectangle in the plane (exact Euclidean, 2D embeddable).
		Files.writeString(input, """
				4
				A	0 3.333333 5.555555 4.444444
				B	3.333333 0 4.444444 5.555555
				C	5.555555 4.444444 0 3.333333
				D	4.444444 5.555555 3.333333 0
				""");

		ByteArrayOutputStream errAfter = new ByteArrayOutputStream();
		ByteArrayOutputStream errBefore = new ByteArrayOutputStream();
		int afterCode = PreprocessDistances.run(new String[]{
				"--input", input.toString(),
				"--output", afterOut.toString(),
				"--mds", "--dimension", "2",
				"--round", "2"
		}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(errAfter));
		int beforeCode = PreprocessDistances.run(new String[]{
				"--input", input.toString(),
				"--output", beforeOut.toString(),
				"--mds", "--dimension", "2",
				"--round", "2",
				"--round-when", "before"
		}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(errBefore));
		assertEquals(0, afterCode);
		assertEquals(0, beforeCode);

		DistFile.Matrix after = DistFile.read(afterOut);
		DistFile.Matrix before = DistFile.read(beforeOut);
		assertTrue(Math.abs(after.distances()[0][1] - before.distances()[0][1]) > 1e-3,
				"before and after order should produce different distances");

		String afterText = Files.readString(afterOut);
		String beforeText = Files.readString(beforeOut);
		assertTrue(afterText.contains("when=after"));
		assertTrue(beforeText.contains("when=before"));
		assertTrue(errAfter.toString().contains("comparedTo=sanitizedInput"));
		assertTrue(errBefore.toString().contains("comparedTo=afterRounding"));
	}

	@Test
	void parseRejectsRoundWhenWithoutBothMdsAndRound() {
		assertThrows(DistPreprocessException.class, () -> PreprocessDistances.parse(new String[]{
				"--input", "a.dist", "--output", "b.dist", "--mds", "--round-when", "before"
		}));
		assertThrows(DistPreprocessException.class, () -> PreprocessDistances.parse(new String[]{
				"--input", "a.dist", "--output", "b.dist", "--round", "2", "--round-when", "after"
		}));
		assertThrows(DistPreprocessException.class, () -> PreprocessDistances.parse(new String[]{
				"--input", "a.dist", "--output", "b.dist", "--round-when", "after"
		}));
	}

	@Test
	void parseRejectsInvalidRoundWhenValue() {
		assertThrows(DistPreprocessException.class, () -> PreprocessDistances.parse(new String[]{
				"--input", "a.dist", "--output", "b.dist",
				"--mds", "--round", "2", "--round-when", "first"
		}));
	}

	@Test
	void refusesToOverwriteInput() throws Exception {
		Path input = tmp.resolve("same.dist");
		Files.writeString(input, """
				2
				a	0 1
				b	1 0
				""");
		int code = PreprocessDistances.run(new String[]{
				"--input", input.toString(),
				"--output", input.toString(),
				"--round", "2"
		}, new PrintStream(new ByteArrayOutputStream()), new PrintStream(new ByteArrayOutputStream()));
		assertEquals(2, code);
	}
}
