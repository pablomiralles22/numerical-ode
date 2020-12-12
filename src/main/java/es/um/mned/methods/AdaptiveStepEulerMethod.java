/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package es.um.mned.methods;

import es.um.mned.methods.AdaptiveStepMethod;
import es.um.mned.ode.InitialValueProblem;

/**
 * Fixed Step Euler Method
 * 
 * @author F. Esquembre
 * @version September 2020
 */
public class AdaptiveStepEulerMethod extends AdaptiveStepMethod {
    private double mCurrentStep;
    private double mMinimumStepAllowed; // Non-convergence minimum
    private double[] mHalfStepState;
    private double[] mFullStepState; 
    
    /**
     * Initializes the method for a given InitialValueProblem
     * @param InitialValueProblem problem 
     * @param step the fixed step to take. If negative, we'd solve backwards in time
     */
    public AdaptiveStepEulerMethod(InitialValueProblem problem, double step, double tolerance) {
        super(problem,step);
        super.setTolerance(tolerance);
        mCurrentStep = step;
        mHalfStepState = problem.getInitialState();
        mFullStepState = problem.getInitialState();
        mMinimumStepAllowed = step/1.0e6;
    }

    /**
     * Extrapolated Euler method implementation
     * @param deltaTime the step to take
     * @param time the current time
     * @param state the current state
     * @return the value of time of the step taken, state will contain the updated state
     */
    public double doStep(double deltaTime, double time, double[] state) {
        while (mCurrentStep>=mMinimumStepAllowed) {
            double[] derivative = mProblem.getDerivative(time, state);
            super.addToEvaluationCounter(2);
            double halfStep = mCurrentStep/2;
            for (int i=0; i<state.length; i++) {
                mHalfStepState[i] = state[i] + halfStep     * derivative[i];
                mFullStepState[i] = state[i] + mCurrentStep * derivative[i];
            }
            derivative = mProblem.getDerivative(time+halfStep, mHalfStepState);
            double error = 0;
            for (int i=0; i<state.length; i++) {
                mHalfStepState[i] += halfStep * derivative[i];
                double errorInIndex = mHalfStepState[i]-mFullStepState[i];
                error += errorInIndex*errorInIndex;
            }
            error = Math.sqrt(error);
            if (error<mTolerance*mCurrentStep) {
                for (int i=0; i<state.length; i++) {
                    //state[i] = mHalfStepState[i]; 
                    state[i] = 2*mHalfStepState[i] - mFullStepState[i];
                }
                time += mCurrentStep;
                // Adapt step
                if (error<1.0e-10) mCurrentStep = 2*mCurrentStep;
                else {
                    double q = 0.84*(mTolerance*mCurrentStep)/error;
                    mCurrentStep *= q;
                }
                //System.out.println ("ACCEPTED: t = "+time+ " New step is "+mCurrentStep+ " error = "+error);
                return time;
            }
            // Try a new smaller step
            double q = 0.84*(mTolerance*mCurrentStep)/error;
            mCurrentStep *= q;
            System.out.println ("REJECTED: t = "+time+ " New step is "+mCurrentStep+ " error = "+error);
        }
        // Was not able to reach tolerance before going below mMinimumStepAllowed
        return Double.NaN; 
    }
    
}
