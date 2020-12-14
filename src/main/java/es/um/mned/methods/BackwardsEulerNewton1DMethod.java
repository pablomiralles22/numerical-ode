package es.um.mned.methods;

import es.um.mned.interpolation.ExtendedStateFunction;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.tools.Newton1D;

public class BackwardsEulerNewton1DMethod extends FixedStepMethod {
	
	private static class 
	BackwardsEuler1DMethodExtendedEquation implements ExtendedStateFunction {

		double t; // current time
		double x; // current state
		double h;
		ExtendedInitialValueProblem ivp;
		
		public BackwardsEuler1DMethodExtendedEquation(ExtendedInitialValueProblem ivp) {
			this.ivp = ivp;
		}
		
		@Override
		public double[] getState(double w) {
			return new double[] {
					w - x - h * (ivp.getDerivative(t+h, new double[]{w} ))[0]
			};
		}

		@Override
		public double getState(double w, int index) {
			return w - x - 
					h * (ivp.getDerivative(t+h, new double[]{w} ))[0];
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
		
		public void changeParameters(double t, double x, double deltaTime) {
			this.t = t;
			this.x = x;
			this.h = deltaTime;
		}
		
	}
	
	BackwardsEuler1DMethodExtendedEquation mEquation;
	double mTolerance;

	public BackwardsEulerNewton1DMethod(ExtendedInitialValueProblem problem, double step, double tolerance) {
		super(problem, step);
		mEquation = new BackwardsEuler1DMethodExtendedEquation(problem);
		mTolerance = tolerance;
	}
	
    
    @Override
    public int getOrder() {
    	return 1;
    }

	@Override
	public double doStep(double deltaTime, double time, double[] state) {
		mEquation.changeParameters(time, state[0], deltaTime);
        try {
			state[0] = Newton1D.solve(mEquation, state[0]);
		} catch (Exception e) {
			System.out.println("Newton didn't converge in BackwardsEuler");
			return Double.NaN;
		}
        return time+deltaTime;
	}

}
