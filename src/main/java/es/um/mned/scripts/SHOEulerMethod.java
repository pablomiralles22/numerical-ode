package es.um.mned.scripts;

import java.util.Optional;
import java.util.Vector;

import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.ConvergenceException;
import es.um.mned.ode.Event;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.problems.SimpleHarmonicOscillator;

public class SHOEulerMethod {
	
    static private double l = 0.7;
    static private double m = 1.0;
    static private double k = 1.5;

    static private double b = 0.; // 0.3
    static private double amp = 0.; // 0.4
    static private double freq = 1.3; // 2.4
    
    static private double Xo = 1.5;
    static private double Vo = 0;
    
    private static class XCross extends Event {
    	
    	private Vector<Double> zeros;

		public XCross() {
			super();
			zeros = new Vector<>();
		}

		@Override
		public double crossFunction(double time, double[] state) {
			return state[0];
		}

		@Override
		public void crossAction(double time, double[] state) {
			System.out.println("=============================");
			System.out.println("X coordenate: " + state[0]);
			System.out.println("Time: " + time);
			System.out.println("=============================");
			zeros.addElement(time);
		}
		
		@Override
		public boolean stopCondition() {
			return false;
		}
		
		public double getFrequency() {
			return zeros.elementAt(2) - zeros.elementAt(0);
		}
    	
    }
    
    public static void main(String[] args) throws ConvergenceException {
    	// Params
        double maxTime = 40;
        double tolerance = 1.0e-2;
        double initialStep = 0.1;
//      int[] indexes = new int[]{0};

    	
    	// Problem
        SimpleHarmonicOscillator shoProblem = new SimpleHarmonicOscillator(
                0., new double[] {Xo, Vo},
                l, m, k,
                b, amp, freq
                );
        
        // Analytical sol
        StateFunction sol = shoProblem.getTrueSol();

    	// -----------------------------------


        // Euler Extrapolation
        NumericalSolution solution = FixedStepEulerMethod.extrapolateToTolerance(shoProblem, maxTime, 
                                                                        tolerance, initialStep, 1.0e-6);
        
        double hStep = solution.get(1).getTime() - solution.get(0).getTime();
        System.out.println ("Solution accepted for h = "+hStep);
        System.out.println ("Evaluations = "+shoProblem.getEvaluationCounter()+"\n");
        shoProblem.resetEvaluationCounter();
        
        // Euler method h
        XCross eulerCrossH = new XCross();
        FixedStepMethod method = new FixedStepEulerMethod(shoProblem, hStep, Optional.of(eulerCrossH));
        method.solve(maxTime);
        System.out.println ("Max Error for h ("+hStep+") is "+ method.getSolution().getMaxError(sol));
        System.out.println("Frequency with h: " + eulerCrossH.getFrequency());
        
        
        // Euler method h/2
        XCross eulerCrossH2 = new XCross();
        FixedStepMethod method2 = new FixedStepEulerMethod(shoProblem,hStep/2, Optional.of(eulerCrossH2));
        method2.solve(maxTime);
        System.out.println("Frequency with h/2: " + eulerCrossH2.getFrequency());
        
        
        // Error comparison
        System.out.println ("Max Error for h/2 (" + hStep/2 + ")is " + solution.getMaxError(sol));
        System.out.println ("Max Error for extrapolation is " + solution.getMaxError(sol));
        System.out.println ("Evaluations = "+shoProblem.getEvaluationCounter()+"\n");
        
        
//        DisplaySolution.listError(solution, sol, indexes,10);
//        DisplaySolution.statePlot(solution, 0, 1, 10);
//        DisplaySolution.timePlot(solution, indexes, 10);
//        DisplaySolution.timePlot(solution);
    }
}
