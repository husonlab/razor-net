/*
 * TimeoutProgress.java Copyright (C) 2026 Daniel H. Huson
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

import jloda.util.CanceledException;
import jloda.util.progress.ProgressListener;
import jloda.util.progress.ProgressObserver;

/**
 * A {@link ProgressListener} that enforces a time budget: once the deadline has passed,
 * {@link #checkForCancel()} throws a {@link DeadlineReachedException}, which is a {@link CanceledException}, so
 * that the algorithm stops and the caller can return the best result found so far. With a budget, the progress
 * bar shows the fraction of the budget that has been used, and the progress values of the algorithm are ignored.
 * Without a budget ({@code maxSeconds <= 0}) everything is delegated.
 *
 * @author Daniel Huson, 2026
 */
public class TimeoutProgress implements ProgressListener {
	/**
	 * Thrown by {@link #checkForCancel()} once the wall-clock budget is exhausted. Distinct from a genuine
	 * user cancellation (which surfaces as a plain {@link CanceledException}).
	 */
	public static class DeadlineReachedException extends CanceledException {
		public DeadlineReachedException(String message) {
			super(message);
		}
	}

	private static final long BAR_UPDATE_INTERVAL_MILLIS = 100L; // throttle UI bar updates to ~10/s

	private final ProgressListener delegate;
	private final long startMillis;
	private final long budgetMillis;   // total budget in ms; <= 0 means no limit (delegate verbatim)
	private final long deadlineMillis; // absolute deadline; <= 0 means no limit
	private volatile boolean deadlineReached = false;
	private volatile long lastBarUpdateMillis = 0L;

	/**
	 * @param delegate   the underlying progress listener that drives the UI and user cancellation
	 * @param maxSeconds wall-clock budget in seconds counted from now; {@code <= 0} disables the deadline
	 */
	public TimeoutProgress(ProgressListener delegate, double maxSeconds) {
		this.delegate = delegate;
		this.startMillis = System.currentTimeMillis();
		this.budgetMillis = maxSeconds > 0 ? Math.round(maxSeconds * 1000.0) : 0L;
		this.deadlineMillis = budgetMillis > 0 ? startMillis + budgetMillis : 0L;
	}

	/**
	 * with a budget, sets the progress bar of the delegate to elapsed/budget (throttled; never throws)
	 */
	private void updateTimeBar(long now) {
		if (budgetMillis <= 0)
			return;
		if (now - lastBarUpdateMillis < BAR_UPDATE_INTERVAL_MILLIS)
			return;
		lastBarUpdateMillis = now;
		var elapsed = Math.max(0L, Math.min(budgetMillis, now - startMillis));
		delegate.setMaximum(budgetMillis);
		delegate.setProgressIgnoreCancel(elapsed);
	}

	/**
	 * @return the seconds left until the deadline, 0 once it has passed, or {@code Double.POSITIVE_INFINITY}
	 * when there is no budget; used to hand a phase a capped share of what remains
	 */
	public double secondsRemaining() {
		if (deadlineMillis <= 0)
			return Double.POSITIVE_INFINITY;
		return Math.max(0.0, (deadlineMillis - System.currentTimeMillis()) / 1000.0);
	}

	/**
	 * @return true once the deadline has been observed to have passed (set the first time it is enforced)
	 */
	public boolean isDeadlineReached() {
		if (!deadlineReached && deadlineMillis > 0 && System.currentTimeMillis() >= deadlineMillis)
			deadlineReached = true;
		return deadlineReached;
	}

	@Override
	public void checkForCancel() throws CanceledException {
		var now = System.currentTimeMillis();
		if (deadlineMillis > 0 && now >= deadlineMillis) {
			deadlineReached = true;
			throw new DeadlineReachedException("RazorNet: time budget exhausted");
		}
		updateTimeBar(now);
		delegate.checkForCancel();
	}

	// ---- delegation; with a budget, the progress values of the algorithm are ignored ----

	@Override
	public void setMaximum(long total) {
		if (budgetMillis > 0)
			return; // time-driven: we own the scale
		delegate.setMaximum(total);
	}

	@Override
	public void setProgress(long current) throws CanceledException {
		if (budgetMillis > 0) { // time-driven: refresh the time bar, ignore the algorithm's value, never throw
			updateTimeBar(System.currentTimeMillis());
			return;
		}
		delegate.setProgress(current);
	}

	@Override
	public void setProgressIgnoreCancel(long current) {
		if (budgetMillis > 0) {
			updateTimeBar(System.currentTimeMillis());
			return;
		}
		delegate.setProgressIgnoreCancel(current);
	}

	@Override
	public long getProgress() {
		return delegate.getProgress();
	}

	@Override
	public void addObserver(ProgressObserver observer) {
		delegate.addObserver(observer);
	}

	@Override
	public void setTasks(String taskName, String subtaskName) {
		delegate.setTasks(taskName, subtaskName);
	}

	@Override
	public void setSubtask(String subtaskName) {
		delegate.setSubtask(subtaskName);
	}

	@Override
	public void setCancelable(boolean enabled) {
		delegate.setCancelable(enabled);
	}

	@Override
	public boolean isUserCancelled() {
		return delegate.isUserCancelled();
	}

	@Override
	public void setUserCancelled(boolean userCancelled) {
		delegate.setUserCancelled(userCancelled);
	}

	@Override
	public void incrementProgress() throws CanceledException {
		if (budgetMillis > 0) {
			updateTimeBar(System.currentTimeMillis());
			return;
		}
		delegate.incrementProgress();
	}

	@Override
	public void incrementProgressIgnoreCancel() {
		if (budgetMillis > 0) {
			updateTimeBar(System.currentTimeMillis());
			return;
		}
		delegate.incrementProgressIgnoreCancel();
	}

	@Override
	public void close() {
		delegate.close();
	}

	@Override
	public boolean isCancelable() {
		return delegate.isCancelable();
	}

	@Override
	public void reportTaskCompleted() {
		delegate.reportTaskCompleted();
	}

	@Override
	public void setDebug(boolean debug) {
		delegate.setDebug(debug);
	}

	@Override
	public void setPause(boolean pause) {
		delegate.setPause(pause);
	}

	@Override
	public boolean getPause() {
		return delegate.getPause();
	}
}
