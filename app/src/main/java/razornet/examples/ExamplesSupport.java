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

import javafx.application.Platform;
import javafx.collections.ListChangeListener;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuItem;
import javafx.scene.control.SeparatorMenuItem;
import jloda.fx.examples.ExampleEntry;
import jloda.fx.examples.ExamplesManager;
import jloda.fx.window.IMainWindow;
import jloda.fx.window.MainWindowManager;
import jloda.fx.window.NotificationManager;
import splitstree6.io.FileLoader;
import splitstree6.window.MainWindow;

import java.io.IOException;
import java.nio.file.Path;

/**
 * the bundled example files and the File -> Examples menu that opens them
 * <p>
 * This class sits in the same package as the example resources, so the resource directory is derived from its own
 * package name. Every resource ends in .dat, including index.dat, which lists the examples in menu order with their
 * display names and real file names; scripts/make-examples.sh writes them from examples/bundled/.
 */
public class ExamplesSupport {
	public static final String RESOURCE_DIR = "/" + ExamplesSupport.class.getPackageName().replace('.', '/');
	private static final String MENU_ID = "examplesMenu";

	private static ExamplesManager manager;

	public static ExamplesManager getManager() {
		if (manager == null)
			manager = ExamplesManager.createOrNull(ExamplesSupport.class, RESOURCE_DIR);
		return manager;
	}

	/**
	 * installs the Examples menu in every main window, those that exist and those to come
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
	 * adds the Examples menu below Open Recent in the window's File menu, unless it is there already
	 */
	public static void install(IMainWindow window) {
		if (!(window instanceof MainWindow mainWindow))
			return;
		var manager = getManager();
		if (manager == null)
			return;
		var controller = mainWindow.getController();
		var items = controller.getFileMenu().getItems();
		if (items.stream().anyMatch(item -> MENU_ID.equals(item.getId())))
			return;
		var menu = new Menu("Examples");
		menu.setId(MENU_ID);
		String previous = null;
		for (var entry : manager.entries()) {
			// a separator wherever the method named before the colon changes ("RazorNet: ...", "RazorHaploNet: ...")
			var name = entry.displayName();
			var group = name.contains(":") ? name.substring(0, name.indexOf(':')) : "";
			if (previous != null && !group.equals(previous))
				menu.getItems().add(new SeparatorMenuItem());
			previous = group;
			var item = new MenuItem(name);
			item.setMnemonicParsing(false); // an underscore in a name (F_ST) is shown as it is
			item.setOnAction(e -> open(mainWindow, entry));
			// SplitsTree's window presenter drops the action of every menu item and disables it whenever the focused
			// tab changes, and its tabs then re-enable the items they handle; this item belongs to no tab, so it
			// puts itself back on the next pulse
			item.disableProperty().addListener((v, o, n) -> {
				if (n)
					Platform.runLater(() -> {
						item.setOnAction(e -> open(mainWindow, entry));
						item.setDisable(false);
					});
			});
			menu.getItems().add(item);
		}
		var pos = items.indexOf(controller.getOpenRecentMenu());
		items.add(pos < 0 ? items.size() : pos + 1, menu);
	}

	/**
	 * writes the example to the temporary directory and opens it like any file: in this window, if its input is
	 * empty, otherwise in a new one
	 */
	private static void open(MainWindow mainWindow, ExampleEntry entry) {
		try {
			var dir = Path.of(System.getProperty("java.io.tmpdir"), "RazorNet-examples");
			var file = getManager().materialize(entry, dir);
			FileLoader.apply(mainWindow, file.toString(), ex -> NotificationManager.showError("Open example failed: " + ex));
		} catch (IOException ex) {
			NotificationManager.showError("Open example failed: " + ex.getMessage());
		}
	}

	private ExamplesSupport() {
	}
}
