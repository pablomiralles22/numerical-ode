package es.um.mned.methods;

import java.util.Arrays;

import es.um.mned.auxiliar.ExtendedFunction1D;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.tools.Newton1D;

public class ImplicitNewton1DFixedStepMethod extends FixedStepMethod {
	
	public static ImplicitNewton1DFixedStepMethod backwardsEulerNewton1DFixedStepMethod(
			ExtendedInitialValueProblem problem,
			double step,
			double tolerance
			) {
		return new ImplicitNewton1DFixedStepMethod(
				problem,
				new double[] {1},
				new double[] {1},
				step,
				tolerance
				);
	}
	
	private static class 
	Implicit1DMethodExtendedEquation implements ExtendedFunction1D {

		protected double t;
		protected double[][] points;
		protected double[][] derivatives;
		private double step;
		private double[] a;
		private double[] b;
		
		ExtendedInitialValueProblem ivp;
		
		public Implicit1DMethodExtendedEquation(ExtendedInitialValueProblem ivp, double[] a, double b[], double step) {
			this.ivp = ivp;
			this.a = Arrays.copyOf(a, a.length);
			this.b = Arrays.copyOf(b, b.length);
			this.step = step;
			points = new double[a.length][];
			points = new double[b.length - 1][];
		}

		@Override
		public double getValue(double w) {
			double value = w;
			
			for(int i=0; i<a.length; ++i)
				value -= a[i] * points[i][0];
			
			for(int i=0; i<b.length-1; ++i)
				value -= step * b[i] * derivatives[i][0];
			
			value -= step * b[b.length - 1] * (ivp.getDerivative(t+step, new double[]{w} ))[0];
			
			return value;
		}

		@Override
		public double getDerivative(double w) {
			return 1 - step * b[b.length - 1] * (ivp.getDerivativeDY(t+step, new double[]{w} ))[0];
		}
		
	}
	
	Implicit1DMethodExtendedEquation mEquation;
	double tolerance;

	private ImplicitNewton1DFixedStepMethod(ExtendedInitialValueProblem problem, double[] a, double[] b, double step, double tolerance) {
		super(problem, step);
		mEquation = new Implicit1DMethodExtendedEquation(problem, a, b, step);
		this.tolerance = tolerance;
	}

	@Override
	public double doStep(double deltaTime, double time, double[] state) {
        try {
			state[0] = Newton1D.solve(mEquation, state[0], tolerance);
		} catch (Exception e) {
			System.out.println("Newton didn't converge in BackwardsEuler");
			return Double.NaN;
		}
        return time+deltaTime;
	}

}
