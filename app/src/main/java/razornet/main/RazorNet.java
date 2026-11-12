/*
 * RazorNet.java Copyright (C) 2026 Daniel H. Huson
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

import javafx.application.Application;
import javafx.geometry.Point2D;
import javafx.stage.Stage;
import jloda.fx.message.MessageWindow;
import jloda.fx.util.ArgsOptions;
import jloda.fx.util.Icebergs;
import jloda.fx.util.ProgramProperties;
import jloda.fx.util.ResourceManagerFX;
import jloda.fx.window.MainWindowManager;
import jloda.fx.window.NotificationManager;
import jloda.fx.window.SplashScreen;
import jloda.fx.window.WindowGeometry;
import jloda.util.Basic;
import jloda.util.ProgramExecutorService;
import jloda.util.UsageException;
import razornet.examples.ExamplesSupport;
import splitstree6.cite.ExtractCitations;
import splitstree6.main.AppProfile;
import splitstree6.main.SplitsTree6;
import splitstree6.window.MainWindow;

import java.io.File;
import java.time.Duration;

/**
 * The RazorNet program: a JavaFX application built on SplitsTree6 that offers the RazorNet algorithms and related
 * methods (see {@link RazorNetProfile}).
 */
public class RazorNet extends Application {
	private static boolean showSplash = true;

	@Override
	public void init() {
		Runtime.getRuntime().addShutdownHook(new Thread(() -> {
			ProgramProperties.store();
		}));
		ProgramProperties.setUseGUI(true);
		ExtractCitations.addKeyCitation(razornet.algorithms.RazorNet.CITATION_KEY, razornet.algorithms.RazorNet.CITATION);
	}

	/**
	 * main
	 */
	public static void main(String[] args) {
		Basic.restoreSystemOut(System.err); // send system out to system err
		Basic.startCollectionStdErr();
		ResourceManagerFX.addResourceRoot(SplitsTree6.class, "splitstree6.resources");
		AppProfile.setProfile(new RazorNetProfile());

		// as in SplitsTree6: icebergs (enlarged hit areas behind shapes), which the network view requires
		Icebergs.setEnabled(true);

		ProgramProperties.getProgramIconsFX().setAll(Utils.getImage(RazorNet.class, "RazorNet-512x512.png"));
		ProgramProperties.setProgramName(Version.NAME);
		ProgramProperties.setProgramVersion(Version.SHORT_DESCRIPTION);
		// the inherited menu item Help -> "Open User Manual..." opens the RazorNet manual
		splitstree6.main.Version.WEBSITE_URL = Version.WEBSITE_URL;
		// the inherited menu item Help -> "Check for Updates..." checks the RazorNet repository and version
		splitstree6.main.Version.HOME_URL = Version.HOME_URL;
		splitstree6.main.Version.UPDATE_VERSION = Version.VERSION;
		ProgramProperties.setProgramLicence("""
				Copyright (C) 2026. This program comes with ABSOLUTELY NO WARRANTY.
				This is free software, licensed under the terms of the GNU General Public License, Version 3.
				Sources available at: https://github.com/husonlab/razor-net
				""");
		SplashScreen.setVersionString(ProgramProperties.getProgramVersion());
		SplashScreen.setImage(Utils.getImage(RazorNet.class, "RazorNet-Splash.png"));
		SplashScreen.setLabelAnchor(new Point2D(190, 10));

		try {
			parseArguments(args);
		} catch (Throwable th) {
			//catch any exceptions and the like that propagate up to the top level
			if (!th.getMessage().startsWith("Help")) {
				System.err.println("Fatal error:" + "\n" + th);
				Basic.caught(th);
			}
			System.exit(1);
		}

		launch(args);
	}

	protected static void parseArguments(String[] args) throws UsageException {
		var options = new ArgsOptions(args, RazorNet.class, "Run the RazorNet on an MSA or distances");
		options.setAuthors("Momoko Hayamizu and Daniel H. Huson");
		options.setLicense(ProgramProperties.getProgramLicence());
		options.setVersion(ProgramProperties.getProgramVersion());

		options.comment(ArgsOptions.OTHER);

		final var propertiesFile = options.getOption("-p", "propertiesFile", "Properties file", getDefaultPropertiesFile());
		final var showVersion = options.getOption("-V", "version", "Show version string", false);
		final var silentMode = options.getOption("-S", "silentMode", "Silent mode", false);
		showSplash = !options.getOption("-s", "hideSplash", "Hide startup splash screen", false);
		ProgramExecutorService.setNumberOfCoresToUse(options.getOption("-t", "threads", "Maximum number of threads to use in a parallel algorithm (0=all available)", 0));
		ProgramProperties.setConfirmQuit(options.getOption("-q", "confirmQuit", "Confirm quit on exit", ProgramProperties.isConfirmQuit()));
		ProgramProperties.put("MaxNumberRecentFiles", 100);
		options.done();

		ProgramProperties.load(propertiesFile);

		// default for RazorNet: split networks are drawn as outlines, unless the user has chosen a diagram
		final var splitsDiagramKey = "splitstree6.view.splits.viewer.SplitsView" + "optionDiagram";
		if (ProgramProperties.get(splitsDiagramKey, "").isBlank())
			ProgramProperties.put(splitsDiagramKey, "Outline");

		if (silentMode) {
			Basic.stopCollectingStdErr();
			Basic.hideSystemErr();
			Basic.hideSystemOut();
		}

		if (showVersion) {
			System.err.println(ProgramProperties.getProgramVersion());
			System.err.println(jloda.util.Version.getVersion(RazorNet.class, ProgramProperties.getProgramName()));
			System.err.println("Java version: " + System.getProperty("java.version"));
		}
	}

	@Override
	public void start(Stage stage) throws Exception {
		if (showSplash)
			SplashScreen.showSplash(Duration.ofSeconds(3));
		try {
			stage.setTitle("Untitled - " + ProgramProperties.getProgramName());
			NotificationManager.setShowNotifications(true);
			ExamplesSupport.installInAllWindows(); // RazorNet's examples in File -> Examples, in this window and every later one
			final var window = new MainWindow();
			WindowGeometry.setToStage(stage);
			window.show(stage, stage.getX(), stage.getY(), stage.getWidth(), stage.getHeight());
			WindowGeometry.listenToStage(stage);
			MainWindowManager.getInstance().addMainWindow(window);
			MessageWindow.getInstance().setVisible(true);
			MessageWindow.getInstance().setVisible(false);
		} catch (Exception ex) {
			Basic.caught(ex);
			throw ex;
		}
	}

	@Override
	public void stop() {
		ProgramProperties.store();
		System.exit(0);
	}

	public static String getDefaultPropertiesFile() {
		if (ProgramProperties.isMacOS())
			return System.getProperty("user.home") + "/Library/Preferences/RazorNet.def";
		else
			return System.getProperty("user.home") + File.separator + ".RazorNet.def";
	}
}
