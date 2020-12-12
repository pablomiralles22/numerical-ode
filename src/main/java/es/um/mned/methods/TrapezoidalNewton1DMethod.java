package es.um.mned.methods;

import es.um.mned.interpolation.ExtendedFunction1D;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.tools.Newton1D;

public class TrapezoidalNewton1DMethod extends FixedStepMethod {
	
	private static class 
	Trapezoidal1DMethodExtendedEquation implements ExtendedFunction1D {

		double t, x, h;
		double derivative;
		ExtendedInitialValueProblem ivp;
		
		public Trapezoidal1DMethodExtendedEquation(ExtendedInitialValueProblem ivp) {
			this.ivp = ivp;
		}

		@Override
		public double getValue(double w) {
			return w - x - h / 2. * (
					(ivp.getDerivative(t+h, new double[]{w} ))[0] + derivative
					);
		}

		@Override
		public double getDerivative(double w) {
			return 1 - h / 2. * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0];
		}
		
		public void changeParameters(double t, double x, double deltaTime, double derivative) {
			this.t = t;
			this.x = x;
			this.h = deltaTime;
			this.derivative = derivative;
		}
		
	}
	
	Trapezoidal1DMethodExtendedEquation mEquation;
	double mTolerance;
	double[] mDerivative;

	public TrapezoidalNewton1DMethod(ExtendedInitialValueProblem problem, double step, double tolerance) {
		super(problem, step);
		mEquation = new Trapezoidal1DMethodExtendedEquation(problem);
		mTolerance = tolerance;
	}

	@Override
	public double doStep(double deltaTime, double time, double[] state) {
		mDerivative = mProblem.getDerivative(time, state);
		mEquation.changeParameters(time, state[0], deltaTime, mDerivative[0]);
        try {
			state[0] = Newton1D.solve(mEquation, state[0], mTolerance);
		} catch (Exception e) {
			System.out.println("Newton didn't converge in BackwardsEuler");
			return Double.NaN;
		}
        return time+deltaTime;
	}

}
