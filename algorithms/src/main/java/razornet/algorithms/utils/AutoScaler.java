/*
 * AutoScaler.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.algorithms.utils;

/**
 * Scales doubles to integers using up to a given number of digits. The scale factor is the smallest power of ten
 * (at most 10^maxDigits) for which min*scale and max*scale are integers, up to rounding error; min may equal max.
 * <p>
 * Daniel Huson, 10.2025
 */
public class AutoScaler {
	private final double min;
	private final double max;
	private final int maxDigits;
	private final double scale;

	/**
	 * Construct an AutoScaler with given range and digit limit.
	 *
	 * @param min       smallest value
	 * @param max       largest value; may equal {@code min}
	 * @param maxDigits maximum number of decimal digits to use for scaling
	 */
	public AutoScaler(double min, double max, int maxDigits) {
		if (max < min)
			throw new IllegalArgumentException("max must be >= min, got min=" + min + ", max=" + max);
		if (maxDigits < 0)
			throw new IllegalArgumentException("maxDigits must be >= 0");
		this.min = min;
		this.max = max;
		this.maxDigits = maxDigits;
		this.scale = determineScale(min, max, maxDigits);
	}

	/**
	 * Chosen scale factor (power of ten).
	 */
	public double getScale() {
		return scale;
	}

	/**
	 * Converts a double to a scaled integer.
	 */
	public int toInt(double d) {
		return (int) Math.round(d * scale);
	}

	/**
	 * Converts a scaled integer back to a double.
	 */
	public double toDouble(int i) {
		return i / scale;
	}

	/**
	 * the smallest power of ten 10^d (d at most maxDigits) that makes both min and max integers, up to rounding error;
	 * 10^maxDigits if there is none
	 */
	private static double determineScale(double min, double max, int maxDigits) {
		if (isIntegerish(min) && isIntegerish(max))
			return 1.0;

		for (int d = 0; d <= maxDigits; d++) {
			double s = Math.pow(10, d);
			if (isIntegerish(min * s) && isIntegerish(max * s))
				return s;
		}
		return Math.pow(10, maxDigits);
	}

	/**
	 * Returns true if value is very close to an integer.
	 */
	private static boolean isIntegerish(double x) {
		double rounded = Math.rint(x);
		return Math.abs(x - rounded) < 1e-10;
	}

	@Override
	public String toString() {
		return String.format("AutoScaler[min=%g, max=%g, maxDigits=%d, scale=%g]",
				min, max, maxDigits, scale);
	}
}