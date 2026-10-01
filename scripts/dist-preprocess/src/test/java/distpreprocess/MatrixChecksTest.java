package distpreprocess;

import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class MatrixChecksTest {
	private static PrintStream sink() {
		return new PrintStream(new ByteArrayOutputStream());
	}

	@Test
	void tinyDiagonalIsZeroed() {
		double[][] d = {
				{1e-15, 1, 2},
				{1, 0, 3},
				{2, 3, 0}
		};
		double[][] out = MatrixChecks.sanitize(d, List.of("a", "b", "c"), sink());
		assertEquals(0.0, out[0][0]);
		assertEquals(1.0, out[0][1]);
	}

	@Test
	void largeDiagonalIsError() {
		double[][] d = {
				{0.01, 1, 2},
				{1, 0, 3},
				{2, 3, 0}
		};
		assertThrows(DistPreprocessException.class,
				() -> MatrixChecks.sanitize(d, List.of("a", "b", "c"), sink()));
	}

	@Test
	void largeAsymmetryIsError() {
		double[][] d = {
				{0, 1, 2},
				{1.5, 0, 3},
				{2, 3, 0}
		};
		assertThrows(DistPreprocessException.class,
				() -> MatrixChecks.sanitize(d, List.of("a", "b", "c"), sink()));
	}

	@Test
	void clearlyNegativeDistanceIsError() {
		double[][] d = {
				{0, -0.01, 2},
				{-0.01, 0, 3},
				{2, 3, 0}
		};
		assertThrows(DistPreprocessException.class,
				() -> MatrixChecks.sanitize(d, List.of("a", "b", "c"), sink()));
	}

	@Test
	void triangleViolationIsReportedAndNotRepaired() {
		double[][] d = {
				{0, 1, 3},
				{1, 0, 1},
				{3, 1, 0}
		};
		double[][] before = copy(d);
		MatrixChecks.TriangleReport report = MatrixChecks.checkTriangleInequality(d);
		assertTrue(report.triangleInequalityViolations() > 0);
		assertEquals(1.0, report.maxTriangleSlack(), 1e-12);
		assertEquals(3.0, d[0][2]);
		assertTrue(matricesEqual(before, d));
	}

	@Test
	void metricTriangleHasNoViolations() {
		double[][] d = {
				{0, 1, 1},
				{1, 0, 1},
				{1, 1, 0}
		};
		assertEquals(0, MatrixChecks.checkTriangleInequality(d).triangleInequalityViolations());
	}

	private static double[][] copy(double[][] d) {
		double[][] c = new double[d.length][d.length];
		for (int i = 0; i < d.length; i++) {
			System.arraycopy(d[i], 0, c[i], 0, d.length);
		}
		return c;
	}

	private static boolean matricesEqual(double[][] a, double[][] b) {
		for (int i = 0; i < a.length; i++) {
			for (int j = 0; j < a.length; j++) {
				if (a[i][j] != b[i][j]) {
					return false;
				}
			}
		}
		return true;
	}
}
