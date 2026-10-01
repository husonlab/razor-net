/**
 * The RazorNet program, its command-line tools and its test suite.
 */
module razornet {
	requires javafx.base;
	requires javafx.controls;
	requires javafx.fxml;
	requires javafx.graphics;
	requires transitive jloda_core;
	requires jloda_fx;
	requires jloda_phylogeny;
	requires transitive splitstreesix;
	requires razornet.algorithms;

	requires java.xml;
	requires java.sql;
	requires jdk.compiler;

	exports razornet.main;
	opens razornet.main;

	exports razornet.testsuite;
	exports razornet.tools;
}