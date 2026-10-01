package distpreprocess;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

class SignificantRoundingTest {
	@Test
	void roundsThreePointThree() {
		assertEquals("3.3", SignificantRounding.formatRounded(3.333333, 2));
	}

	@Test
	void roundsSmallValue() {
		assertEquals("0.0033", SignificantRounding.formatRounded(0.003333, 2));
	}

	@Test
	void zeroStaysZero() {
		assertEquals("0", SignificantRounding.formatRounded(0.0, 3));
	}

	@Test
	void applyCopiesUpperTriangleAndZeroDiagonal() {
		double[][] in = {
				{0, 3.333333, 1.25},
				{3.333333, 0, 0.003333},
				{1.25, 0.003333, 0}
		};
		double[][] out = SignificantRounding.apply(in, 2);
		assertEquals(0.0, out[0][0]);
		assertEquals(0.0, out[1][1]);
		assertEquals(3.3, out[0][1], 1e-12);
		assertEquals(3.3, out[1][0], 1e-12);
		assertEquals(0.0033, out[1][2], 1e-12);
		assertEquals(0.0033, out[2][1], 1e-12);
	}

	@Test
	void rejectsNonPositiveDigits() {
		assertThrows(DistPreprocessException.class, () -> SignificantRounding.roundValue(1.0, 0));
	}
}
