/*
 * Version.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.main;

/**
 * name, version and URLs of the RazorNet program
 */
public class Version {
	static public final String NAME = "RazorNet";
	public static final String SUFFIX = ".stree6";
	static public String HOME_URL = "https://github.com/husonlab/razor-net"; // sources and releases; "Check for Updates" reads the latest release

	public static final String VERSION = resolveVersion();
	public static final String SHORT_DESCRIPTION = NAME + " (version " + VERSION + ") - License GPL v3";

	static public final String WEBSITE_URL = "https://husonlab.github.io/razor-net"; // the user manual

	public static String resolveVersion() {
		// set by the packaged launcher (jpackage --java-options -Dapp.version=...)
		var sys = System.getProperty("app.version");
		if (sys != null && !sys.isBlank())
			return sys;
		// class-path fallback (e.g. running the plain jar): the Implementation-Version from the jar manifest
		var pkg = razornet.main.Version.class.getPackage();
		var v = (pkg != null) ? pkg.getImplementationVersion() : null;
		return (v != null && !v.isBlank()) ? v : "dev";
	}
}
