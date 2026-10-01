/*
 * ExtensionTrace.java Copyright (C) 2026 Daniel H. Huson
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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;

/**
 * Optional diagnostic recorder for one RazorNet run, used by {@link AnalysisLog}: records the calls of the
 * recursive extension for the first restart that reaches it, and the length achieved by every restart. Off
 * unless {@link #start()} has been called.
 * <p>
 * Daniel Huson and Momoko Hayamizu, 8.2026
 */
public final class ExtensionTrace {
	private static volatile boolean enabled = false;
	private static final AtomicBoolean firstRestartClaimed = new AtomicBoolean(false);
	private static final List<String> calls = new ArrayList<>();
	private static final List<Double> restartLengths = new ArrayList<>();
	private static volatile int maxMatrixSize = 0;

	private ExtensionTrace() {
	}

	/** Begins recording, discarding anything held from an earlier run. */
	public static void start() {
		synchronized (calls) {
			calls.clear();
		}
		synchronized (restartLengths) {
			restartLengths.clear();
		}
		firstRestartClaimed.set(false);
		maxMatrixSize = 0;
		enabled = true;
	}

	public static void stop() {
		enabled = false;
	}

	public static boolean isEnabled() {
		return enabled;
	}

	/**
	 * True for exactly one restart -- the first to ask. Callers use it to decide whether to record their
	 * recursion, so the trace covers one restart end to end rather than a mixture of several.
	 */
	public static boolean claimTracedRestart() {
		return enabled && firstRestartClaimed.compareAndSet(false, true);
	}

	/** One {@code extendSubset} call of the traced restart. */
	public static void recordCall(int depth, int subsetSize, int slacklistSize, String slackValues, int pointsAdded,
								  int skipped, int verticesG, int keptFromSubset, int edgesG, int settled,
								  int unconfirmed, String componentSizes, int matrixSize) {
		if (!enabled) return;
		var line = "%5d %8d %10d %8d %8d %8d %9d %8d %8d %10d %8d  %-18s %s"
				.formatted(depth, subsetSize, slacklistSize, pointsAdded, skipped, verticesG, keptFromSubset,
						edgesG, settled, unconfirmed, matrixSize, componentSizes, slackValues);
		synchronized (calls) {
			calls.add(line);
		}
	}

	/** a call that returned at one of the stopping conditions, with the reason */
	public static void recordEarlyReturn(int depth, int subsetSize, int matrixSize, String reason) {
		if (!enabled) return;
		var line = "%5d %8d %10s %8s %8s %8s %9s %8s %8s %10s %8d  %-18s %s"
				.formatted(depth, subsetSize, "-", "-", "-", "-", "-", "-", "-", "-", matrixSize, "-", "STOPPED: " + reason);
		synchronized (calls) {
			calls.add(line);
		}
	}

	public static void recordRestartLength(double length) {
		if (!enabled) return;
		synchronized (restartLengths) {
			restartLengths.add(length);
		}
	}

	public static void noteMatrixSize(int size) {
		if (enabled && size > maxMatrixSize) maxMatrixSize = size;
	}

	public static List<String> callLines() {
		synchronized (calls) {
			return new ArrayList<>(calls);
		}
	}

	public static List<Double> restartLengths() {
		synchronized (restartLengths) {
			return new ArrayList<>(restartLengths);
		}
	}

	public static int maxMatrixSize() {
		return maxMatrixSize;
	}

	/** Column header matching {@link #recordCall}. */
	public static String callHeader() {
		return "%5s %8s %10s %8s %8s %8s %9s %8s %8s %10s %8s  %-18s %s"
				.formatted("depth", "|S|", "|slacklist|", "added", "skipped", "|V(G)|", "keptFromS", "|E(G)|",
						"settled", "unconfirm", "|D|", "components", "c(x) values");
	}
}
