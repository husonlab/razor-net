/*
 * TimeLimit.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.testsuite;

import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Runs a computation with a time limit. A computation that exceeds the limit is abandoned: its thread is
 * interrupted and left to end with the program.
 *
 * @author Daniel Huson, 2026
 */
public final class TimeLimit {
	/**
	 * @param seconds the limit; 0 or less means none
	 * @return the result, or null if the limit was exceeded
	 */
	public static <T> T call(Callable<T> computation, int seconds) throws Exception {
		if (seconds <= 0)
			return computation.call();
		var executor = Executors.newSingleThreadExecutor(runnable -> {
			var thread = new Thread(runnable);
			thread.setDaemon(true);
			return thread;
		});
		var future = executor.submit(computation);
		try {
			return future.get(seconds, TimeUnit.SECONDS);
		} catch (TimeoutException e) {
			future.cancel(true);
			return null;
		} catch (ExecutionException e) {
			if (e.getCause() instanceof Exception cause)
				throw cause;
			throw e;
		} finally {
			executor.shutdown();
		}
	}
}
