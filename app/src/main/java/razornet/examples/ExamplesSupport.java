/*
 * ExamplesSupport.java Copyright (C) 2026 Daniel H. Huson
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

package razornet.examples;

import javafx.collections.ListChangeListener;
import javafx.scene.control.MenuItem;
import jloda.fx.examples.ExamplesManager;
import jloda.fx.window.IMainWindow;
import jloda.fx.window.MainWindowManager;
import jloda.fx.window.NotificationManager;
import splitstree6.io.FileLoader;
import splitstree6.window.MainWindow;

import java.io.IOException;
import java.util.ArrayList;

/**
 * the example documents bundled with RazorNet: they replace SplitsTree's own in the File -> Examples menu of every window
 * <p>
 * This class sits in the same package as the example resources, so the resource directory is derived from its own
 * package name. Every resource ends in .dat, including index.dat, which lists the examples in menu order with their
 * display names and real file names; scripts/make-examples.sh writes them from examples/bundled/. The examples are
 * complete SplitsTree6 documents (.stree6) with their networks computed; like SplitsTree's, they are loaded from memory
 * and never written to disk.
 */
public class ExamplesSupport {
	public static final String RESOURCE_DIR = "/" + ExamplesSupport.class.getPackageName().replace('.', '/');

	private static ExamplesManager manager;

	public static ExamplesManager getManager() {
		if (manager == null)
			manager = ExamplesManager.createOrNull(ExamplesSupport.class, RESOURCE_DIR);
		return manager;
	}

	/**
	 * puts RazorNet's examples into the Examples menu of every main window, those that exist and those to come
	 */
	public static void installInAllWindows() {
		var windows = MainWindowManager.getInstance().getMainWindows();
		windows.forEach(ExamplesSupport::install);
		windows.addListener((ListChangeListener<IMainWindow>) change -> {
			while (change.next()) {
				if (change.wasAdded())
					change.getAddedSubList().forEach(ExamplesSupport::install);
			}
		});
	}

	/**
	 * replaces the entries of the window's Examples menu, which SplitsTree's presenter has filled with its own, by
	 * RazorNet's. The items are set directly: the menu's disabled state is bound by the presenter to the running workflow
	 */
	public static void install(IMainWindow window) {
		if (!(window instanceof MainWindow mainWindow))
			return;
		var manager = getManager();
		if (manager == null)
			return;
		var items = new ArrayList<MenuItem>();
		for (var entry : manager.entries()) {
			var item = new MenuItem(entry.displayName());
			item.setMnemonicParsing(false);
			item.setOnAction(e -> {
				try {
					FileLoader.applyText(mainWindow, entry.fileName(), manager.readText(entry), ex -> NotificationManager.showError("Open example failed: " + ex));
				} catch (IOException ex) {
					NotificationManager.showError("Open example failed: " + ex.getMessage());
				}
			});
			items.add(item);
		}
		mainWindow.getController().getExamplesMenu().getItems().setAll(items);
	}

	private ExamplesSupport() {
	}
}
