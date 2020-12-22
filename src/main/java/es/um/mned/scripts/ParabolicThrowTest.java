package es.um.mned.scripts;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.Event;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.problems.ParabolicThrowWithFriction;
import es.um.mned.utils.ConvergenceException;

public class ParabolicThrowTest {
	
	static double t0 = 0.;
	static double x0 = 0, vx0 = 100;
	static double y0 = 300, vy0 = 0.;
	static double mGravity = 9.8;
	
	static final double analyticalZeroAt = Math.sqrt(2. * y0 / 9.8);
	
    static public class TrueSol implements StateFunction {
        public double[] getState(double time) {
            return new double[] { 
            		vx0 * time,
            		vx0,
            		y0 - 0.5*mGravity*time*time,
            		-mGravity
            		};
        }
        public double getState(double time, int index) {
            switch (index) {
                case 0 : return vx0 * time;
                case 1 : return vx0;
                case 2 : return y0 - 0.5*mGravity*time*time;
                case 3 : return -mGravity;
                default : return Double.NaN;
            }
        }
            
        public double yZeroAt() {
            return Math.sqrt(2*y0/mGravity);
        }

}
	
    private static class YCross extends Event {

		public YCross(boolean blocking, double tolerance) {
			super(blocking, tolerance);
		}

		@Override
		public double crossFunction(double time, double[] state) {
			return state[2];
		}

		@Override
		public void crossAction(double time, double[] state) {
			System.out.println("=============================");
			System.out.println("Zero at: " + time);
			System.out.println("Analytical zero at: " + analyticalZeroAt);
			double relError = (time - analyticalZeroAt) /  analyticalZeroAt;
			System.out.println("Relative error: " + relError);
			System.out.println("=============================");
		}
    	
    }
	
    public static void main(String[] args) throws ConvergenceException {
    	// Params
    	double frictionCoefficient = 0.;
    	double maxTime = 15.;
    	double hStep = 1e-4;
    	double zeroTol = 1e-8;
    	// Problem
        InitialValueProblem problem = 
        		new ParabolicThrowWithFriction(t0, new double[] { x0, vx0, y0, vy0 }, frictionCoefficient);
        // Event
        Event yCross = new YCross(true, zeroTol);
        // Method
        FixedStepMethod method = new FixedStepEulerMethod(problem, hStep, yCross);
        // Analytical sol
        StateFunction sol = new TrueSol();

        
        method.solve(maxTime);
        System.out.println ("Evaluations = " + problem.getEvaluationCounter());
        System.out.println("Max error: " + method.getSolution().getMaxError(sol, new int[] {0,2}));
        
        
//        DisplaySolution.listError(method.getSolution(), sol, new int[]{0, 2});
//        DisplaySolution.statePlot(method.getSolution(), 0, 2);
//        DisplaySolution.timePlot(method.getSolution(), new int[]{1,3});
//        DisplaySolution.timePlot(method.getSolution());
    }

}
