package es.um.mned.methods;

import java.util.Optional;

import es.um.mned.interpolation.ExtendedStateFunction;
import es.um.mned.ode.ConvergenceException;
import es.um.mned.ode.Event;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.utils.Newton1D;

public class FixedStepBackwardsEulerNewton1DMethod extends FixedStepMethod {
	
	/**
	 * Auxiliar class to use in Newton method
	 */
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
	
	/*
	 * ========================================
	 * Attributes
	 * ========================================
	 */
	
	BackwardsEuler1DMethodExtendedEquation equation;
	double tolerance = 0.0; // This will use Newton's default tol
	
	/*
	 * ========================================
	 * Constructors
	 * ========================================
	 */

	public FixedStepBackwardsEulerNewton1DMethod(ExtendedInitialValueProblem problem, double step, Optional<Event> event) {
		super(problem, step, event);
		equation = new BackwardsEuler1DMethodExtendedEquation(problem, step);
	}
	
	/*
	 * ========================================
	 * Step, order and tol
	 * ========================================
	 */

	/**
	 * Set tolerance for Newton method (advanced users)
	 * @param tolerance
	 */
	public void setNewtonTolerance(double tolerance) {
		this.tolerance = tolerance;
	}
    
    @Override
    public int getOrder() {
    	return 1;
    }

	@Override
	public double doStep(double deltaTime, double time, double[] state) throws ConvergenceException {
		equation.x = state[0];
		equation.t = time;
		state[0] = Newton1D.solve(equation, state[0], tolerance);
        return time+deltaTime;
	}

}
