/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.methods;

import java.util.Iterator;
import java.util.Optional;

import es.um.mned.ode.ConvergenceException;
import es.um.mned.ode.Event;
import es.um.mned.ode.InitialValueProblem;
import es.um.mned.ode.NumericalSolution;
import es.um.mned.ode.NumericalSolutionPoint;

public class FixedStepEulerMethod extends FixedStepMethod {
    
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public FixedStepEulerMethod(InitialValueProblem problem, double step, Optional<Event> event) {
        super(problem,step, event);
    }
    
    @Override
    public int getOrder() {
    	return 1;
    }

    /**
     * Euler method implementation
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     */
    public double doStep(double deltaTime, double time, double[] state) {
        double[] derivative = mProblem.getDerivative(time, state);
        for (int i=0; i<state.length; i++) {
            state[i] = state[i] + deltaTime * derivative[i];
        }
        return time+deltaTime;
    }
        
    static public NumericalSolution extrapolate (NumericalSolution fullStep, NumericalSolution halfStep) {
        NumericalSolution extrapolatedSolution = NumericalSolution.createExtrapolationSol(fullStep);
        Iterator<NumericalSolutionPoint> iteratorFull = fullStep.iterator();
        Iterator<NumericalSolutionPoint> iteratorHalf = halfStep.iterator();
        while (iteratorFull.hasNext() && iteratorHalf.hasNext()) {
            NumericalSolutionPoint pointFull = iteratorFull.next();
            double[] stateFull = pointFull.getState();
            double[] stateHalf = iteratorHalf.next().getState();
            for (int i=0; i<stateFull.length; i++) {
                stateFull[i] = (2*stateHalf[i]-stateFull[i]); // Reuse stateFull array
            }
            extrapolatedSolution.add(pointFull.getTime(), stateFull);
            if (!iteratorHalf.hasNext()) return extrapolatedSolution;
            iteratorHalf.next();
        }
        return extrapolatedSolution;
    }

    
    /**
     * Uses Richardson extrapolation for Euler method
     * @param problem
     * @param maxTime
     * @param tolerance
     * @param initialStep
     * @param minStepAllowed
     * @return 
     * @throws ConvergenceException 
     */
    static public NumericalSolution extrapolateToTolerance(InitialValueProblem problem, 
            double maxTime, double tolerance, 
            double initialStep, double minStepAllowed) throws ConvergenceException {
        
        double h = initialStep;
        FixedStepEulerMethod methodFull = new FixedStepEulerMethod(problem,h, Optional.empty());
        methodFull.solve(maxTime);
        NumericalSolution solutionFull = methodFull.getSolution();
        while (Math.abs(h)>Math.abs(minStepAllowed)) {
            System.out.println ("Trying for h = "+h+"...");
            FixedStepEulerMethod methodHalf = new FixedStepEulerMethod(problem,h/2, Optional.empty());
            methodHalf.solve(maxTime);
            NumericalSolution solutionHalf = methodHalf.getSolution();
            double maxError = maxHalfStepError(solutionFull,solutionHalf);
            System.out.println ("- Error (for h/2) ~= "+maxError);
            if (maxError<tolerance) {
                System.out.println ("Tolerance reached for h = "+h+"\n");
                return extrapolate(solutionFull,solutionHalf);
            }
            h /= 2;
            solutionFull = solutionHalf;
        }
        return null;
    }

}
