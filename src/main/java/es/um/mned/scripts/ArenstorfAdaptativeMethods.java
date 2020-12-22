package es.um.mned.scripts;

import es.um.mned.methods.*;
import es.um.mned.ode.*;
import es.um.mned.problems.ArenstorfOrbits;
import es.um.mned.utils.*;

public class ArenstorfAdaptativeMethods {


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
			System.out.println("X coordenate: " + state[0]);
			System.out.println("Time: " + time);
			System.out.println("=============================");
		}
    	
    }
    
    public static void main(String[] args) throws ConvergenceException {
    	// Parameters
        double hStep = -1.0e-2;
        double tolerance = 1.0e-8;
        double maxT = -ArenstorfOrbits.PERIOD * 2.3;
        // Problem
        InitialValueProblem problem = new ArenstorfOrbits(
        		0.,
        		new double[] { 0.994, 0.0 , 0.0, -2.00158510637908252240537862224 }
        		);
        
        // Event
        Event yCross = new YCross(false, tolerance);
        
        // Method
//        FixedStepMethod method = new FixedStepPredictorCorrector4Method(problem,hStep, yCross);
//        FixedStepMethod method = new AdaptiveStepPredictorCorrector4Method(problem,hStep, tolerance, yCross);
//        FixedStepMethod method = new AdaptiveStepRKFehlbergMethod(problem,hStep, tolerance, yCross);
        FixedStepMethod method = new AdaptiveStepRK4Method(problem,hStep, tolerance, yCross);
//        FixedStepMethod method = new FixedStepModifiedEulerMethod(problem,hStep, yCross);
        
        
        method.solve(maxT);
        
        System.out.println ("Evaluations = " + problem.getEvaluationCounter());

        DisplaySolution.statePlot(method.getSolution(), 0, 2,(int) Math.floor(1.0e-2/hStep));
        if (method instanceof AdaptiveStepMethod) 
            DisplaySequence.plot(((AdaptiveStepMethod) method).getSolution().getStepList());
    }
}
