package es.um.mned.ode;

public interface ExtendedInitialValueProblem extends InitialValueProblem {
	
	public double[] getDerivativeDY(double time, double[] state);
	
}
