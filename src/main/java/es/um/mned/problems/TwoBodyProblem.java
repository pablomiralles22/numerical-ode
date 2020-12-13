/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.problems;

import java.util.Arrays;
import es.um.mned.interpolation.StateFunction;
import es.um.mned.interpolation.*;
import es.um.mned.methods.*;
import es.um.mned.ode.*;
import es.um.mned.tools.*;

/**
 *
 * @author paco
 */
public class TwoBodyProblem extends InitialValueProblem {

	private static double[] initState = new double[] { 152.100533, 0.0 , 0.0, 0.105444 }; // x,vx,y,vy

	static private double sG = 8.6498928e-4;

    private double mM1 = 1988.5, mM2 = 0.0059724;
    private double mConstant = sG * (mM1 + mM2);
    
    // ------------------
    // Implementation of InitialValueProblem
    // ------------------

    public TwoBodyProblem(double t0, double[] x0) {
		super(t0, x0);
	}
    
    public double[] getDerivative(double t, double[] x) {
        double div  = Math.pow(x[0]*x[0]+x[2]*x[2],1.5);
        return new double[] { x[1], -mConstant * x[0] / div, 
                              x[3], -mConstant * x[2] / div };
    }

    // ------------------
    // End of implementation of InitialValueProblem
    // ------------------

    static private double computeCrossing (String label, InitialValueProblem problem, FixedStepMethod method,
            NumericalSolutionPoint fromPoint, NumericalSolutionPoint toPoint, 
            double tolerance, int index) {
//        StateFunction interpolator = new EulerMethodInterpolator(problem, fromPoint);
//      StateFunction interpolator = new FixedStepMethodInterpolator(method, fromPoint);
        StateFunction interpolator = new HermiteInterpolator(problem, fromPoint, toPoint);
        double zeroAt = BisectionMethod.findZero (interpolator, fromPoint.getTime(), toPoint.getTime(), tolerance, index);
        if (Double.isNaN(zeroAt)) {
            System.out.print (label+" not found!!!");
        }
        else {
            double xZero = interpolator.getState(zeroAt, 0);
            System.out.println (label+" at t="+zeroAt+", x = "+xZero);
        }        
        return zeroAt;
    }
    
    public static void main(String[] args) {
        double hStep = 10;
        double tolerance = 1.0e-8;

        InitialValueProblem problem = new TwoBodyProblem(0., Arrays.copyOf(initState, initState.length));
        FixedStepMethod method = new FixedStepModifiedEulerMethod(problem,10);
        method = new FixedStepPredictorCorrector4Method(problem,10);
        method = new AdaptiveStepPredictorCorrector4Method(problem,hStep, tolerance);
        method = new AdaptiveStepRKFehlbergMethod(problem,hStep, tolerance);

        NumericalSolutionPoint previousPoint, currentPoint;
        
        previousPoint = currentPoint = method.getSolution().getLastPoint();
        int day = 0;
        int maxYears = 10, currentYear = 0;
        double lastYearStart = 0;
        boolean done = false;
        while (day<370*(maxYears+1) && currentYear<=maxYears) {
            previousPoint = currentPoint;
            currentPoint = method.step();
            if (currentPoint.getState(2)*previousPoint.getState(2)<0) {
                String label = currentPoint.getState(3)>0 ? "Aphelium  " : "Perihelium";
                
                double yearTime = computeCrossing(label,problem, method, previousPoint, currentPoint, 1.0e-6, 2);
                // Crossed axis
                if (currentPoint.getState(3)>0) {
                    double currentYearHours = yearTime-lastYearStart;
                    currentYear++;
                    lastYearStart = yearTime;
                    int days  = (int) Math.floor(currentYearHours/24.);
                    double hours = (currentYearHours - days*24.);
                    System.out.println ("====================================");
                    System.out.println ("YEAR: "+currentYear);
                    System.out.println ("Days : "+days+ " : Hours = " +hours);
                    System.out.println ("Days : "+currentYearHours/24.);
                    System.out.println ("Hours : "+currentYearHours);
                    System.out.println ("Total Hours : "+currentPoint.getTime());
                    System.out.println ("====================================");
                }
            }
            double time = currentPoint.getTime();
            if (time>day*24) {
                day++;
            }
        }
        previousPoint.println();
        currentPoint.println();
        
        
        System.out.println ("Evaluations ="+method.getEvaluationCounter());

        DisplaySolution.statePlot(method.getSolution(), 0, 2);
        if (method instanceof AdaptiveStepMethod) 
            DisplaySequence.plot(((AdaptiveStepMethod) method).getStepList());
    }
}
