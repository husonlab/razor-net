package distpreprocess;

import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

/**
 * Command-line entry point: optional classical MDS and/or significant-digit rounding.
 * When both are enabled, {@link RoundWhen} selects whether rounding runs before or after MDS.
 */
public final class PreprocessDistances {
	private PreprocessDistances() {
	}

	public enum RoundWhen {
		AFTER,
		BEFORE;

		static RoundWhen parse(String value) {
			return switch (value) {
				case "after" -> AFTER;
				case "before" -> BEFORE;
				default -> throw new DistPreprocessException(
						"--round-when must be 'after' or 'before', found: " + value);
			};
		}

		String cliValue() {
			return name().toLowerCase();
		}
	}

	public record Options(
			Path input,
			Path output,
			boolean mds,
			int dimension,
			boolean round,
			int significantDigits,
			RoundWhen roundWhen,
			boolean roundWhenExplicit
	) {
	}

	public static void main(String[] args) {
		int code = run(args, System.out, System.err);
		if (code != 0) {
			System.exit(code);
		}
	}

	public static int run(String[] args, PrintStream out, PrintStream err) {
		try {
			Options options = parse(args);
			execute(options, out, err);
			return 0;
		} catch (HelpRequested e) {
			return 0;
		} catch (DistPreprocessException e) {
			err.println("error: " + e.getMessage());
			return 2;
		} catch (IOException e) {
			err.println("error: " + e.getMessage());
			return 1;
		}
	}

	public static Options parse(String[] args) {
		Path input = null;
		Path output = null;
		boolean mds = false;
		boolean dimensionSet = false;
		int dimension = 2;
		boolean round = false;
		int significantDigits = -1;
		RoundWhen roundWhen = RoundWhen.AFTER;
		boolean roundWhenExplicit = false;

		for (int i = 0; i < args.length; i++) {
			String arg = args[i];
			try {
				switch (arg) {
					case "--help", "-h" -> {
						printUsage(System.err);
						throw new HelpRequested();
					}
					case "--input", "-i" -> input = Path.of(requireValue(args, ++i, arg));
					case "--output", "-o" -> output = Path.of(requireValue(args, ++i, arg));
					case "--mds" -> mds = true;
					case "--dimension", "-d" -> {
						dimension = Integer.parseInt(requireValue(args, ++i, arg));
						dimensionSet = true;
					}
					case "--round", "-r" -> {
						round = true;
						significantDigits = Integer.parseInt(requireValue(args, ++i, arg));
					}
					case "--round-when" -> {
						roundWhen = RoundWhen.parse(requireValue(args, ++i, arg));
						roundWhenExplicit = true;
					}
					default -> throw new DistPreprocessException("Unknown argument: " + arg);
				}
			} catch (NumberFormatException e) {
				throw new DistPreprocessException("Bad numeric value near '" + arg + "'");
			} catch (ArrayIndexOutOfBoundsException e) {
				throw new DistPreprocessException("Missing value after " + arg);
			}
		}

		if (input == null) {
			throw new DistPreprocessException("Missing --input");
		}
		if (output == null) {
			throw new DistPreprocessException("Missing --output");
		}
		if (!mds && !round) {
			throw new DistPreprocessException("Specify --mds and/or --round; refusing a no-op copy");
		}
		if (dimensionSet && !mds) {
			throw new DistPreprocessException("--dimension requires --mds");
		}
		if (mds && dimension < 1) {
			throw new DistPreprocessException("dimension must be >= 1, found " + dimension);
		}
		if (round && significantDigits < 1) {
			throw new DistPreprocessException("significantDigits must be >= 1, found " + significantDigits);
		}
		if (roundWhenExplicit && !(mds && round)) {
			throw new DistPreprocessException("--round-when requires both --mds and --round");
		}
		return new Options(input, output, mds, dimension, round, significantDigits, roundWhen, roundWhenExplicit);
	}

	public static void execute(Options options, PrintStream out, PrintStream err) throws IOException {
		Path input = options.input().toAbsolutePath().normalize();
		Path output = options.output().toAbsolutePath().normalize();
		if (!Files.isRegularFile(input)) {
			throw new DistPreprocessException("Input file does not exist: " + input);
		}
		if (input.equals(output)) {
			throw new DistPreprocessException("Output path must differ from the input file");
		}

		DistFile.Matrix matrix = DistFile.read(input);
		double[][] current = MatrixChecks.sanitize(matrix.distances(), matrix.labels(), err);
		List<String> labels = matrix.labels();

		boolean both = options.mds() && options.round();
		if (both && options.roundWhen() == RoundWhen.BEFORE) {
			current = applyRounding(current, labels, options.significantDigits(), err);
			current = applyMds(current, labels, options.dimension(), "afterRounding", err);
		} else {
			if (options.mds()) {
				current = applyMds(current, labels, options.dimension(), "sanitizedInput", err);
			}
			if (options.round()) {
				current = applyRounding(current, labels, options.significantDigits(), err);
			}
		}

		List<String> header = buildHeader(options);
		boolean useRoundedText = options.round()
				&& !(both && options.roundWhen() == RoundWhen.BEFORE);
		writeOutput(output, header, labels, current, options.significantDigits(), useRoundedText);
		out.println("Wrote: " + output);
	}

	private static double[][] applyMds(double[][] input, List<String> labels, int dimension,
									   String comparedTo, PrintStream err) {
		err.printf("classicalMDS: dimension=%d n=%d comparedTo=%s%n", dimension, labels.size(), comparedTo);
		ClassicalMds.Result mds = ClassicalMds.apply(input, dimension, err);
		double[][] out = mds.distances();
		MatrixChecks.checkTriangleInequality(out).print(err, "afterMDS", labels);
		return out;
	}

	private static double[][] applyRounding(double[][] input, List<String> labels, int significantDigits,
											PrintStream err) {
		err.printf("rounding: significantDigits=%d%n", significantDigits);
		double[][] out = SignificantRounding.apply(input, significantDigits);
		MatrixChecks.checkTriangleInequality(out).print(err, "afterRounding", labels);
		return out;
	}

	private static List<String> buildHeader(Options options) {
		List<String> header = new ArrayList<>();
		header.add("# dist-preprocess");
		header.add("# source=" + options.input().getFileName());
		header.add("# mds=" + options.mds() + (options.mds() ? " dimension=" + options.dimension() : ""));
		if (options.round()) {
			String roundLine = "# round=true significantDigits=" + options.significantDigits();
			if (options.mds()) {
				roundLine += " when=" + options.roundWhen().cliValue();
			}
			header.add(roundLine);
		} else {
			header.add("# round=false");
		}
		return header;
	}

	private static void writeOutput(Path output, List<String> header, List<String> labels, double[][] matrix,
									int significantDigits, boolean useRoundedText) throws IOException {
		DistFile.write(output, header, labels, matrix, (row, col, value) -> {
			if (row == col) {
				return "0";
			}
			if (useRoundedText) {
				return SignificantRounding.formatRounded(value, significantDigits);
			}
			return Double.toString(value);
		});
	}

	private static String requireValue(String[] args, int index, String flag) {
		if (index >= args.length) {
			throw new DistPreprocessException("Missing value after " + flag);
		}
		return args[index];
	}

	static void printUsage(PrintStream err) {
		err.print("""
				PreprocessDistances - apply optional classical MDS and/or significant-digit rounding to a .dist file

				  -i, --input FILE         input .dist (required)
				  -o, --output FILE        output .dist (required; input is never overwritten)
				      --mds                enable classical MDS
				  -d, --dimension N        MDS dimension (default: 2; requires --mds)
				  -r, --round N            round each distance to N significant digits
				      --round-when after|before
				                           when both --mds and --round are set: run rounding after MDS (default)
				                           or before MDS; requires both --mds and --round
				  -h, --help               show this help

				Default order with both enabled: input -> MDS -> round -> output.
				Fails if --mds is set and there are fewer than N positive eigenvalues.
				""");
	}

	/** Used only so --help can leave run() without being treated as a failure by tests that catch it. */
	static final class HelpRequested extends DistPreprocessException {
		HelpRequested() {
			super("help");
		}
	}
}
