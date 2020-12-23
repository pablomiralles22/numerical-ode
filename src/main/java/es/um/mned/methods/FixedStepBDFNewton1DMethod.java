package es.um.mned.methods;

import es.um.mned.interpolation.ExtendedStateFunction;
import es.um.mned.ode.Event;
import es.um.mned.ode.ExtendedInitialValueProblem;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.Newton1D;
import es.um.mned.methods.FixedStepTrapezoidalNewton1DMethod.Trapezoidal1DMethodExtendedEquation;

public class FixedStepBDFNewton1DMethod extends FixedStepMethod {
	
	/**
	 * Auxiliar class to use in Newton method
	 */
	protected static class 
	BDF1DMethodExtendedEquation implements ExtendedStateFunction {

		protected double t, h;
        protected double[] states; // w_{i}, w_{i-1}, ...
        protected double[] a; // a[i]*w_{i}, a[i-1]* w_{i-1}, ...
        protected double b;
		ExtendedInitialValueProblem ivp;
		
		public BDF1DMethodExtendedEquation (ExtendedInitialValueProblem ivp, double step, double[] a, double b) {
			this.ivp = ivp;
            this.a = a;
            this.b = b;
			h = step;
            states = new double[a.length];
		}

		@Override
		public double[] getState(double w) {
            double result = w;
            for(int i=0; i < a.length; ++i)
                result += a[i] * states[i];
            result -= b * h * (ivp.getDerivative(t+h, new double[]{w} ))[0];

			return new double[] {
				result
			};
		}
		
		@Override
		public double getState(double w, int index) {
            double result = w;
            for(int i=0; i < a.length; ++i)
                result += a[i] * states[i];
            result -= b * h * (ivp.getDerivative(t+h, new double[]{w} ))[0];

			return result;
		}

		@Override
		public double[] getDerivative(double w) {
			return new double[] {
				1 - b * h * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0]
			};
		}
		
		@Override
		public double getDerivative(double w, int index) {
			return 1 - b * h * (ivp.getDerivativeDY(t+h, new double[]{w} ))[0];
		}
		
	}
	
	/*
	 * ========================================
	 * Attributes
	 * ========================================
	 */
	BDF1DMethodExtendedEquation equation;
    Trapezoidal1DMethodExtendedEquation startEquation;
    int order;
    int startSteps;
	double tolerance = 0.0; // this will use Newton's default tol
	
	/*
	 * ========================================
	 * Constructors
	 * ========================================
	 */

	/**
	 * @param problem InitialValueProblem that implements partial derivative of y
	 * @param order Order of the BDF method, between 2 and 6. Order 1 is 
	 * 	implemented elsewhere (Backwards Euler)
	 * @param step Size of the steps to take
	 * @param tolerance Tolerance for Newton method to solve the equation
	 */
	public FixedStepBDFNewton1DMethod(
            ExtendedInitialValueProblem problem,
            int order,
            double step
            ) {

		super(problem, step);

        this.order = order;
        startSteps = order - 1;

        double[] a;
        double b;
        switch(order) {
            case 2:
                a = new double[] {-4./3., 1./3.};
                b = 2./3.;
                break;
            case 3:
                a = new double[] {-18./11., 9./11., -2./11.};
                b = 6./11.;
                break;
            case 4:
                a = new double[] {-48./25., 36./25., -16./25., 3./25.};
                b = 12./25.;
                break;
            case 5:
                a = new double[] {-300./137., 300./137., -200./137., 75./137., -12./137.};
                b = 60./137.;
                break;
            case 6:
                a = new double[] {-360./147., 450./147., -400./147., 225./147., -72./147., 10./147.};
                b = 60./147.;
                break;
            default:
            	throw new IllegalArgumentException("Order not available.");
        }
		equation = new BDF1DMethodExtendedEquation(problem, step, a, b);
		/*
		 * this kills BDFs of order > 3, but I don't have better implicit methods.
		 * if I had them I would initialize the first points and startEquation inside the switch.
		 * Notice that order 3 still works well because error in each trapezoidal step has O(h^3) error.
		 * I don't think it's a good idea to use RK methods, as it would lead to issues with rigid problems.
		 */
        startEquation = new Trapezoidal1DMethodExtendedEquation(problem, step);
	}
	

	/**
	 * @param problem InitialValueProblem that implements partial derivative of y
	 * @param order Order of the BDF method, between 1 and 6.
	 * @param step Size of the steps to take
	 * @param tolerance Tolerance for Newton method to solve the equation
	 * @param event Event object
	 */
	public FixedStepBDFNewton1DMethod(
            ExtendedInitialValueProblem problem,
            int order,
            double step,
            Event event
            ) {
		this(problem, order, step);
		super.setEvent(event);
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
    	return order;
    }

	@Override
	public double doStep(double deltaTime, double time, double[] state) throws ConvergenceException {
        if(startSteps > 0) { // method start
        	equation.states[--startSteps] = state[0];
            startEquation.derivative = mProblem.getDerivative(time, state)[0];
            startEquation.t = time;
            startEquation.x = state[0];
            state[0] = Newton1D.solve(startEquation, state[0], tolerance);
            return time+deltaTime;
        }

		equation.t = time;
		
        for(int i = equation.states.length - 1; i>0; --i)
            equation.states[i] = equation.states[i-1];
        equation.states[0] = state[0];
		
		state[0] = Newton1D.solve(equation, state[0], tolerance);
        return time+deltaTime;
	}

}
