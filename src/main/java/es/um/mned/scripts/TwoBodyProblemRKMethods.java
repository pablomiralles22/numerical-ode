package es.um.mned.scripts;

import java.util.Arrays;

import es.um.mned.methods.AdaptiveStepMethod;
import es.um.mned.methods.AdaptiveStepRKFehlbergMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.Event;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.problems.TwoBodyProblem;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.DisplaySequence;
import es.um.mned.utils.DisplaySolution;

public class TwoBodyProblemRKMethods {

	private static double[] initState = new double[] { 152.100533, 0.0 , 0.0, 0.105444 }; // x,vx,y,vy
    
    private static class YCross extends Event {
    	
    	private int loopCount;

		public YCross(double tolerance) {
			super(tolerance);
			loopCount = 0;
		}

		@Override
		public double crossFunction(double time, double[] state) {
			return state[2];
		}

		@Override
		public void crossAction(double time, double[] state) {
			System.out.println("=============================");
			System.out.println("X coordenate: " + state[0]);
			System.out.println("Year: " + (int)(time/(24*365)));
			System.out.println("Day: " + (int)(time/(24) % 365 + 1));
			System.out.println("Hour: " + (int)(time) % 24);
			System.out.println("Time: " + time);
			System.out.println("=============================");
			loopCount++;
		}
		
		@Override
		public boolean stopCondition() {
			return loopCount >= 20;
		}
    	
    }
	
    public static void main(String[] args) {
    	// Parameters
        double hStep = -10;
        double tolerance = 1.0e-8;
        int maxYears = 12;
        double tMax = -maxYears * 365 * 24;
        
        // Problem
        InitialValueProblem problem = new TwoBodyProblem(0., Arrays.copyOf(initState, initState.length));
        // Event
        Event yCross = new YCross(tolerance);

        // Methods
//        FixedStepMethod method = new FixedStepModifiedEulerMethod(problem,hStep, yCross);
//        FixedStepMethod method = new FixedStepPredictorCorrector4Method(problem,hStep,yCross);
//        FixedStepMethod method = new AdaptiveStepEulerMethod(problem,hStep, tolerance, yCross);
//        FixedStepMethod method = new AdaptiveStepPredictorCorrector4Method(problem,hStep, tolerance, yCross);
        FixedStepMethod method = new AdaptiveStepRKFehlbergMethod(problem,hStep, tolerance, yCross);
        
        
        
        try {
			method.solve(tMax);
		} catch (ConvergenceException e) {
			e.printStackTrace();
		}
        
        
        System.out.println ("Evaluations = "+problem.getEvaluationCounter());
        	

        DisplaySolution.statePlot(method.getSolution(), 0, 2);
        if (method instanceof AdaptiveStepMethod) 
            DisplaySequence.plot(((AdaptiveStepMethod) method).getSolution().getStepList());
    }

}
