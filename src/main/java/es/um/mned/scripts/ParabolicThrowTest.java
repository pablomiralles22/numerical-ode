package es.um.mned.scripts;

import es.um.mned.interpolation.EulerMethodInterpolator;
import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.FixedStepEulerMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolutionPoint;
import es.um.mned.problems.ParabolicThrowWithFriction;
import es.um.mned.problems.ParabolicThrowWithFriction.TrueSol;
import es.um.mned.utils.BisectionMethod;
import es.um.mned.utils.ConvergenceException;
import es.um.mned.utils.DisplaySolution;

public class ParabolicThrowTest {
	
    public static void main(String[] args) throws ConvergenceException {
        InitialValueProblem problem = new ParabolicThrowWithFriction(0., new double[] { 0, 100, 300, 0 }, 0.);
        FixedStepMethod method = new FixedStepEulerMethod(problem,1.0e-2);
        
        NumericalSolutionPoint previousPoint, currentPoint;
        
        if (false) method.solve(15);
        else { 
            previousPoint = currentPoint = method.getSolution().getLastPoint();
            while (currentPoint.getState(2)>0) {
                previousPoint = currentPoint;
                currentPoint = method.step();
            }
            previousPoint.println();
            currentPoint.println();
        }
        
        TrueSol sol = new TrueSol();
        if (false) { // find zero
            StateFunction interpolator = new EulerMethodInterpolator(problem, previousPoint);
            double zeroYAt;
			try {
				zeroYAt = BisectionMethod.findZero (interpolator, previousPoint.getTime(), currentPoint.getTime(), 1.0e-8, 2);
			} catch (ConvergenceException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
            if (Double.isNaN(zeroYAt)) {
                System.out.print ("Zero not found!!!");
            }
            else{
                double yZero = interpolator.getState(zeroYAt, 2);
                System.out.println ("Zero at t="+zeroYAt+", y = "+yZero);
                System.out.println ("Analytical zero at t="+sol.yZeroAt());
            }
        }
        System.out.println ("Evaluations = "+problem.getEvaluationCounter());

        DisplaySolution.listError(method.getSolution(), new TrueSol(), new int[]{0, 2});

        if (false) {

            DisplaySolution.statePlot(method.getSolution(), 0, 2);
            DisplaySolution.timePlot(method.getSolution(), new int[]{1,3});
            DisplaySolution.timePlot(method.getSolution());
        }
    }

}
