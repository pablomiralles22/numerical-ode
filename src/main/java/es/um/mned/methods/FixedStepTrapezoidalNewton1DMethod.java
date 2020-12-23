package es.um.mned.methods;

import es.um.mned.interpolation.ExtendedStateFunction;
import es.um.mned.ode.Event;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.Newton1D;

public class FixedStepTrapezoidalNewton1DMethod extends FixedStepMethod {
	
	/**
	 * Auxiliar class to use in Newton method
	 */
	protected static class 
	Trapezoidal1DMethodExtendedEquation implements ExtendedStateFunction {

		protected double t, x, h;
		protected double derivative;
		ExtendedInitialValueProblem ivp;
		
		public Trapezoidal1DMethodExtendedEquation(ExtendedInitialValueProblem ivp, double step) {
			this.ivp = ivp;
			h = step;
		}

		@Override
		public double[] getState(double w) {
			return new double[] {
				w - x - h / 2. * ((ivp.getDerivative(t+h, new double[]{w} ))[0] + derivative)
			};
		}
		
		@Override
		public double getState(double w, int index) {
			return w - x - h / 2. * ((ivp.getDerivative(t+h, new double[]{w} ))[0] + derivative);
		}

		@Override
		public double[] getDerivative(double w) {
			return new double[] {
					1 - h / 2. * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0]
			};
		}
		
		@Override
		public double getDerivative(double w, int index) {
			return 1 - h / 2. * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0];
		}
		
	}
	
	/*
	 * ========================================
	 * Attributes
	 * ========================================
	 */
	
	Trapezoidal1DMethodExtendedEquation mEquation;
	double mTolerance;
	
	/*
	 * ========================================
	 * Constructors
	 * ========================================
	 */

	public FixedStepTrapezoidalNewton1DMethod(ExtendedInitialValueProblem problem, double step, double tolerance) {
		super(problem, step);
		mEquation = new Trapezoidal1DMethodExtendedEquation(problem, step);
		mTolerance = tolerance;
	}
	
	public FixedStepTrapezoidalNewton1DMethod(ExtendedInitialValueProblem problem, double step, double tolerance, Event event) {
		this(problem, step, tolerance);
		super.setEvent(event);
	}
	
	/*
	 * ========================================
	 * Step and order
	 * ========================================
	 */
	
	@Override
    public int getOrder() {
    	return 2;
    }

	@Override
	public double doStep(double deltaTime, double time, double[] state) throws ConvergenceException {
		mEquation.derivative = mProblem.getDerivative(time, state)[0];
		mEquation.t = time;
		mEquation.x = state[0];
		
		state[0] = Newton1D.solve(mEquation, state[0], mTolerance);
        return time+deltaTime;
	}

}
