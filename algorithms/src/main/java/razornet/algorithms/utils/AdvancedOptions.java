/*
 * AdvancedOptions.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.fx.util.ProgramProperties;

/**
 * The switch between the basic and the advanced options of the RazorNet algorithms. The advanced options are
 * listed if the system property {@code razornet.advanced} is true, if the program property
 * {@code RazorNetAdvancedOptions} is true, or after {@link #enable()} has been called, as the command-line tools
 * do.
 */
public final class AdvancedOptions {
	public static final String SYSTEM_PROPERTY = "razornet.advanced";
	public static final String PROGRAM_PROPERTY = "RazorNetAdvancedOptions";

	private static volatile boolean enabledByCode = false;

	private AdvancedOptions() {
	}

	/** whether the advanced options are to be listed */
	public static boolean enabled() {
		return enabledByCode || Boolean.getBoolean(SYSTEM_PROPERTY) || ProgramProperties.get(PROGRAM_PROPERTY, false);
	}

	/** turns advanced mode on for this process; the command-line tools call this so that every option stays reachable */
	public static void enable() {
		enabledByCode = true;
	}
}
