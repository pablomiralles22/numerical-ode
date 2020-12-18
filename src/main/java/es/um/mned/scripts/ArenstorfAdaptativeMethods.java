package es.um.mned.scripts;

import es.um.mned.interpolation.HermiteInterpolator;
import es.um.mned.interpolation.StateFunction;
import es.um.mned.methods.AdaptiveStepMethod;
import es.um.mned.methods.AdaptiveStepPredictorCorrector4Method;
import es.um.mned.methods.AdaptiveStepRKFehlbergMethod;
import es.um.mned.methods.FixedStepMethod;
import es.um.mned.methods.FixedStepPredictorCorrector4Method;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolutionPoint;
import es.um.mned.problems.ArenstorfOrbits;
import es.um.mned.utils.BisectionMethod;
import es.um.mned.utils.DisplaySequence;
import es.um.mned.utils.DisplaySolution;

public class ArenstorfAdaptativeMethods {


    static private double computeCrossing (InitialValueProblem problem, FixedStepMethod method,
            NumericalSolutionPoint fromPoint, NumericalSolutionPoint toPoint, 
            double tolerance, int index) {
//        StateFunction interpolator = new EulerMethodInterpolator(problem, fromPoint);
//        StateFunction interpolator = new FixedStepMethodInterpolator(method, fromPoint);
        StateFunction interpolator = new HermiteInterpolator(problem, fromPoint, toPoint);
        double zeroAt = BisectionMethod.findZero (interpolator, fromPoint.getTime(), toPoint.getTime(), tolerance, index);
        if (Double.isNaN(zeroAt)) {
            System.out.print ("Zero not found!!!");
        }
        else {
            double[] zeroState = interpolator.getState(zeroAt);
            System.out.println ("Zero at t="+zeroAt+", x="+zeroState[0]+", vx="+zeroState[1]+", y="+zeroState[2]+", vy="+zeroState[3]);
        }        
        return zeroAt;
    }
    
    public static void main(String[] args) {
        double hStep = 1.0e-2;
        double tolerance = 1.0e-8;
        InitialValueProblem problem = new ArenstorfOrbits(
        		0.,
        		new double[] { 0.994, 0.0 , 0.0, -2.00158510637908252240537862224 }
        		);
        FixedStepMethod method = new FixedStepPredictorCorrector4Method(problem,hStep);
        method = new AdaptiveStepPredictorCorrector4Method(problem,hStep, tolerance);
        method = new AdaptiveStepRKFehlbergMethod(problem,hStep, tolerance);
        
//        FixedStepMethod method = new FixedStepModifiedEulerMethod(problem,hStep);
        
        NumericalSolutionPoint previousPoint, currentPoint;
        
        previousPoint = currentPoint = method.getSolution().getLastPoint();
        double time = problem.getInitialTime();
        while (time<ArenstorfOrbits.PERIOD*2.3) {
            previousPoint = currentPoint;
            currentPoint = method.step();
            if (currentPoint==null) {
                System.out.println ("Method failed at t="+previousPoint.getTime()+" !!!");
                System.exit(2);

            }
            //if ((time>0.9*sPeriod && time<1.1*sPeriod) || (time>1.8*sPeriod && time<2.1*sPeriod)) 
            {
                if (currentPoint.getState(2)*previousPoint.getState(2)<0) {
                    computeCrossing(problem, method, previousPoint, currentPoint, 1.0e-6, 2);
                    // Crossed axis
                    //if (currentPoint.getTime()>sPeriod) break;
                }
            }
            time = currentPoint.getTime();
        }
        previousPoint.println();
        currentPoint.println();
        
        
        System.out.println ("Evaluations = " + problem.getEvaluationCounter());

        DisplaySolution.statePlot(method.getSolution(), 0, 2,(int) Math.floor(1.0e-2/hStep));
        if (method instanceof AdaptiveStepMethod) 
            DisplaySequence.plot(((AdaptiveStepMethod) method).getSolution().getStepList());
    }
}
