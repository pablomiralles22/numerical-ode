package es.um.mned.scripts;

import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.problems.SimpleHarmonicOscillator;
import es.um.mned.problems.SimpleHarmonicOscillator.TrueSol;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.DisplaySolution;

public class SHOEulerMethod {

    
    static private double Xo = 1.5;
    static private double Vo = 0;
    
    public static void main(String[] args) throws ConvergenceException {
        SimpleHarmonicOscillator shoProblem = new SimpleHarmonicOscillator(0., new double[] {Xo, Vo});
        //FixedStepMethod method = new FixedStepEulerMethod(shoProblem,1.0e-4);

        double maxTime = 40;
        double tolerance = 1.0e-2;
        double initialStep = 0.1;
        
        NumericalSolution solution = FixedStepEulerMethod.extrapolateToTolerance(shoProblem, maxTime, 
                                                                        tolerance, initialStep, 1.0e-6);
        TrueSol sol = new TrueSol();
        
        int[] indexes = new int[]{0};
        double hStep = solution.get(1).getTime() - solution.get(0).getTime();
        System.out.println ("Solution accepted for h = "+hStep);
        System.out.println ("Evaluations = "+shoProblem.getEvaluationCounter()+"\n");
        shoProblem.resetEvaluationCounter();

        FixedStepMethod method = new FixedStepEulerMethod(shoProblem,hStep);
        method.solve(maxTime);
        System.out.println ("Max Error for h ("+hStep+") is "+DisplaySolution.maximumError(method.getSolution(), sol, indexes));
        FixedStepMethod method2 = new FixedStepEulerMethod(shoProblem,hStep/2);
        method2.solve(maxTime);
        System.out.println ("Max Error for h/2 ("+hStep/2+")is "+DisplaySolution.maximumError(method2.getSolution(), sol, indexes));
        System.out.println ("Max Error for extrapolation is "+DisplaySolution.maximumError(solution, sol, indexes));
        System.out.println ("Evaluations = "+shoProblem.getEvaluationCounter()+"\n");
        DisplaySolution.listError(solution, new TrueSol(), indexes,10);

        
        if (false) {

            DisplaySolution.statePlot(solution, 0, 1, 10);
            DisplaySolution.timePlot(solution, indexes, 10);
//            DisplaySolution.timePlot(solution);
        }
    }
}
