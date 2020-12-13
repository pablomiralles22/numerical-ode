package es.um.mned.ode;

public abstract class ExtendedInitialValueProblem extends InitialValueProblem {
	
	public ExtendedInitialValueProblem(double t0, double[] x0) {
		super(t0, x0);
	}

	public abstract double[] getDerivativeDY(double time, double[] state);
	
}
