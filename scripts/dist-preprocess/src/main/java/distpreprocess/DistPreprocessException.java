package distpreprocess;

/**
 * User-facing failure: bad input, insufficient MDS rank, or invalid arguments.
 * These errors must stop the program; they are never repaired silently.
 */
public class DistPreprocessException extends RuntimeException {
	public DistPreprocessException(String message) {
		super(message);
	}
}
