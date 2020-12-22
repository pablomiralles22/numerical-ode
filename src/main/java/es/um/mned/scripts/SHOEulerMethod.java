package es.um.mned.scripts;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.problems.SimpleHarmonicOscillator;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.DisplaySolution;

public class SHOEulerMethod {

    
    static private double Xo = 1.5;
    static private double Vo = 0;
    
    public static void main(String[] args) throws ConvergenceException {
    	// Params
        double maxTime = 40;
        double tolerance = 1.0e-2;
        double initialStep = 0.1;
//      int[] indexes = new int[]{0};

    	
    	// Problem
        SimpleHarmonicOscillator shoProblem = new SimpleHarmonicOscillator(0., new double[] {Xo, Vo});
        
        // Analytical sol
        StateFunction sol = new SimpleHarmonicOscillator.TrueSol();

    	// -----------------------------------


        // Euler Extrapolation
        NumericalSolution solution = FixedStepEulerMethod.extrapolateToTolerance(shoProblem, maxTime, 
                                                                        tolerance, initialStep, 1.0e-6);
        
        double hStep = solution.get(1).getTime() - solution.get(0).getTime();
        System.out.println ("Solution accepted for h = "+hStep);
        System.out.println ("Evaluations = "+shoProblem.getEvaluationCounter()+"\n");
        shoProblem.resetEvaluationCounter();
        
        // Euler method h

        FixedStepMethod method = new FixedStepEulerMethod(shoProblem,hStep);
        method.solve(maxTime);
        System.out.println ("Max Error for h ("+hStep+") is "+ method.getSolution().getMaxError(sol));
        // Euler method h/2
        FixedStepMethod method2 = new FixedStepEulerMethod(shoProblem,hStep/2);
        method2.solve(maxTime);
        
        System.out.println ("Max Error for h/2 (" + hStep/2 + ")is " + solution.getMaxError(sol));
        System.out.println ("Max Error for extrapolation is " + solution.getMaxError(sol));
        System.out.println ("Evaluations = "+shoProblem.getEvaluationCounter()+"\n");
        
        
//        DisplaySolution.listError(solution, sol, indexes,10);
//        DisplaySolution.statePlot(solution, 0, 1, 10);
//        DisplaySolution.timePlot(solution, indexes, 10);
//        DisplaySolution.timePlot(solution);
    }
}
