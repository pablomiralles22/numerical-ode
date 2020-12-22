package es.um.mned.methods;

import es.um.mned.interpolation.ExtendedStateFunction;
import es.um.mned.ode.Event;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.Newton1D;

public class FixedStepBackwardsEulerNewton1DMethod extends FixedStepMethod {
	
	protected static class 
	BackwardsEuler1DMethodExtendedEquation implements ExtendedStateFunction {

		double t; // current time
		double x; // current state
		double h;
		ExtendedInitialValueProblem ivp;
		
		public BackwardsEuler1DMethodExtendedEquation(ExtendedInitialValueProblem ivp, double step) {
			this.ivp = ivp;
			h = step;
		}
		
		@Override
		public double[] getState(double w) {
			return new double[] {
				w - x - h * (ivp.getDerivative(t+h, new double[]{w} ))[0]
			};
		}

		@Override
		public double getState(double w, int index) {
			return w - x - h * (ivp.getDerivative(t+h, new double[]{w} ))[0];
		}

		@Override
		public double[] getDerivative(double w) {
			return new double[] {
				1 - h * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0]
			};
		}
		
		@Override
		public double getDerivative(double w, int index) {
			return 1 - h * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0];
		}
		
	}
	
	BackwardsEuler1DMethodExtendedEquation mEquation;
	double mTolerance;

	public FixedStepBackwardsEulerNewton1DMethod(ExtendedInitialValueProblem problem, double step, double tolerance) {
		super(problem, step);
		mEquation = new BackwardsEuler1DMethodExtendedEquation(problem, step);
		mTolerance = tolerance;
	}
	
	public FixedStepBackwardsEulerNewton1DMethod(ExtendedInitialValueProblem problem, double step, double tolerance, Event event) {
		this(problem, step, tolerance);
		super.setEvent(event);
	}
	
    
    @Override
    public int getOrder() {
    	return 1;
    }

	@Override
	public double doStep(double deltaTime, double time, double[] state) throws ConvergenceException {
		mEquation.x = state[0];
		mEquation.t = time;
		state[0] = Newton1D.solve(mEquation, state[0]);
        return time+deltaTime;
	}

}
