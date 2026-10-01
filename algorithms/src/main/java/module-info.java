/**
 * The RazorNet algorithms for SplitsTree6 and the RazorNet app.
 * <p>
 * Registered through {@code splitstree6.algorithms.AlgorithmProvider}: on the module path via the
 * {@code provides} clause below, on the class path via {@code META-INF/services}.
 */
module razornet.algorithms {
	requires splitstreesix;

	provides splitstree6.algorithms.AlgorithmProvider with razornet.RazorNetProvider;

	exports razornet;
	exports razornet.algorithms;
	exports razornet.algorithms.strict;
	exports razornet.algorithms.relaxed;
	exports razornet.algorithms.utils;
	exports razornet.algorithms.preprocessing;
	exports razornet.algorithms.tcs;
	exports razornet.algorithms.tightspan;
}
