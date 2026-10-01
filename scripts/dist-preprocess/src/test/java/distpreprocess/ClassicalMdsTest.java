package distpreprocess;

import org.junit.jupiter.api.Test;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ClassicalMdsTest {
	private static PrintStream sink() {
		return new PrintStream(new ByteArrayOutputStream());
	}

	@Test
	void recoversEuclideanDistancesInTwoDimensions() {
		// Points (0,0), (1,0), (0,1)
		double s = Math.sqrt(2.0);
		double[][] din = {
				{0, 1, 1},
				{1, 0, s},
				{1, s, 0}
		};
		ClassicalMds.Result result = ClassicalMds.apply(din, 2, sink());
		assertTrue(result.error().RMSE() < 1e-8, "RMSE=" + result.error().RMSE());
		assertTrue(result.error().relativeRMSE() < 1e-8);
		assertTrue(result.error().maxAbsError() < 1e-8);
		assertTrue(result.eigen().positiveCount() >= 2);
		assertTrue(result.eigen().proportionPositive()[0] > 0);
	}

	@Test
	void collinearPointsFailForDimensionTwo() {
		// Points 0, 1, 3 on a line: intrinsically 1-dimensional
		double[][] din = {
				{0, 1, 3},
				{1, 0, 2},
				{3, 2, 0}
		};
		assertThrows(DistPreprocessException.class, () -> ClassicalMds.apply(din, 2, sink()));
	}

	@Test
	void collinearPointsSucceedForDimensionOne() {
		double[][] din = {
				{0, 1, 3},
				{1, 0, 2},
				{3, 2, 0}
		};
		ClassicalMds.Result result = ClassicalMds.apply(din, 1, sink());
		assertTrue(result.error().RMSE() < 1e-8, "RMSE=" + result.error().RMSE());
	}

	@Test
	void dimensionAboveNMinusOneFails() {
		double[][] din = {
				{0, 1},
				{1, 0}
		};
		assertThrows(DistPreprocessException.class, () -> ClassicalMds.apply(din, 2, sink()));
	}
}
