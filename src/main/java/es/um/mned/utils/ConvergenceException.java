package es.um.mned.utils;

public class ConvergenceException extends Exception {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public ConvergenceException(String message) {
		super(message);
	}

	public ConvergenceException(Throwable cause) {
		super(cause);
	}

	public ConvergenceException(String message, Throwable cause) {
		super(message, cause);
	}

	public ConvergenceException(String message, Throwable cause, boolean enableSuppression,
			boolean writableStackTrace) {
		super(message, cause, enableSuppression, writableStackTrace);
	}

}
