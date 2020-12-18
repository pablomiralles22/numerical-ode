package es.um.mned.utils;

import es.um.mned.interpolation.ExtendedStateFunction;

/*
 * This method assumes f is a one-dimensional state function
 */

public class Newton1D {
	
	private static final double DEFAULT_TOL = 1e-8;
	private static final int MAX_IT = 30;

	
	public static double solve(
			ExtendedStateFunction f,
			double start,
			double tol
			) throws Exception {
		double x = start;
		
		for(int i = 0; i < MAX_IT; ++i) {
			double next = x - f.getState(x, 0) / f.getDerivative(x, 0);
			if(Math.abs(x - next) < tol)
				return next;
			x = next;
		}
		
		throw new Exception("Newton method did not converge.");
	}
	
	public static double solve(
			ExtendedStateFunction f,
			double start
			) throws Exception {
		return solve(f, start, DEFAULT_TOL);
	}
	
}
